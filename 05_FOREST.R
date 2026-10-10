#!/usr/bin/env Rscript
# FOREST.R
# FOREST: Feature-frequency Out-of-fold Resampling with Elastic-net Selection and Training.
# 单文件运行：input.csv -> 训练集内差异/频率选特征 -> 重复训练与验证。
# 第 1 列 ID，第 2 列 Tumor/Control，第 3 列起为数值特征。
# 默认 source() 后立即运行；参数可在当前页面先定义 frequency_config 或 n_repeats。
# 直接使用固定默认参数运行，无需先进行参数预筛选或读取已有分析结果。
# 用 options(frequency_oof.functions_only = TRUE) 仅加载函数。
# 详细使用方式、输出定义和失败处理见随附 README.md 的第 5 节。

frequency_default_config <- function() {
  list(
    input_file = "input.csv",
    output_dir = ".",
    positive_class = "Tumor",
    train_fraction = 0.67,
    TOP_n = 10L,
    alpha_fixed = 0.5,
    lambda_ratio_min = 0.001,
    n_repeats = 100L,               # 外层拆分和每次训练集内频率子抽样共用此次数
    n_lambda = 50L,
    pvalue_cutoff = 0.05,            # NULL 表示不做 p 值筛选
    inner_fraction = 0.8,
    diff_method = "wilcox",
    frequency_alpha_mode = "pooled",
    cv_nfolds = 5L,
    plot_iter = 1L,
    seed = 123L,
    coef_tol = 1e-8,
    fit_maxit = 1000000L,
    fpr_grid = seq(0, 1, length.out = 101L),
    ensure_oof_coverage = TRUE,
    require_complete_oof = TRUE,
    install_missing_packages = TRUE,
    save_models = TRUE,
    overwrite = FALSE,
    verbose = TRUE
  )
}

run_frequency_holdouts <- function(
    data,
    positive_class,
    train_fraction_grid,
    alpha_fixed_grid,
    lambda_ratio_min_grid,
    TOP_n,
    n_repeats,
    n_lambda,
    pvalue_cutoff,
    inner_fraction = 0.8,
    diff_method = c("wilcox", "welch"),
    frequency_alpha_mode = c("pooled", "by_alpha"),
    cv_nfolds = 5L,
    plot_iter = 1L,
    seed = 123L,
    coef_tol = 1e-8,
    fit_maxit = 1000000L,
    keep_differential = TRUE,
    verbose = TRUE,
    keep_all_repeats = FALSE,
    ensure_oof_coverage = FALSE) {

  for (pkg in c("glmnet", "pROC")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("请先安装 R 包：", pkg, call. = FALSE)
    }
  }
  diff_method <- match.arg(diff_method)
  frequency_alpha_mode <- match.arg(frequency_alpha_mode)
  is_count <- function(z, lower) {
    is.numeric(z) && length(z) > 0L && all(is.finite(z)) &&
      all(z >= lower) && all(z == floor(z)) &&
      all(z <= .Machine$integer.max)
  }
  stopifnot(
    is.data.frame(data), ncol(data) >= 4L, "Group" %in% names(data),
    !anyNA(data$Group), !anyDuplicated(names(data)),
    length(positive_class) == 1L, !anyNA(positive_class),
    is_count(TOP_n, 2L),
    is_count(n_repeats, 2L), length(n_repeats) == 1L,
    is_count(n_lambda, 2L), length(n_lambda) == 1L,
    is_count(cv_nfolds, 3L), length(cv_nfolds) == 1L,
    is_count(plot_iter, 1L), length(plot_iter) == 1L, plot_iter <= n_repeats,
    is.logical(keep_all_repeats), length(keep_all_repeats) == 1L,
    !anyNA(keep_all_repeats),
    is.logical(ensure_oof_coverage), length(ensure_oof_coverage) == 1L,
    !anyNA(ensure_oof_coverage),
    is_count(seed, 0L), length(seed) == 1L,
    is_count(fit_maxit, 1L), length(fit_maxit) == 1L,
    is.numeric(train_fraction_grid), length(train_fraction_grid) > 0L,
    all(is.finite(train_fraction_grid)),
    all(train_fraction_grid > 0 & train_fraction_grid < 1),
    is.numeric(alpha_fixed_grid), length(alpha_fixed_grid) > 0L,
    all(is.finite(alpha_fixed_grid)),
    all(alpha_fixed_grid > 0 & alpha_fixed_grid <= 1),
    is.numeric(lambda_ratio_min_grid), length(lambda_ratio_min_grid) > 0L,
    all(is.finite(lambda_ratio_min_grid)),
    all(lambda_ratio_min_grid > 0 & lambda_ratio_min_grid < 1),
    length(inner_fraction) == 1L, is.finite(inner_fraction),
    inner_fraction > 0, inner_fraction < 1,
    is.null(pvalue_cutoff) || (is.numeric(pvalue_cutoff) &&
      length(pvalue_cutoff) == 1L && is.finite(pvalue_cutoff) &&
      pvalue_cutoff > 0 && pvalue_cutoff < 1),
    length(coef_tol) == 1L, is.finite(coef_tol), coef_tol >= 0,
    is.logical(keep_differential), length(keep_differential) == 1L,
    !is.na(keep_differential),
    is.logical(verbose), length(verbose) == 1L, !is.na(verbose)
  )
  data <- as.data.frame(data)
  feat_cols <- names(data)[seq.int(3L, ncol(data))]
  stopifnot(
    !any(c("Group", "ID") %in% feat_cols),
    !anyNA(feat_cols), all(nzchar(feat_cols)),
    all(vapply(data[, feat_cols, drop = FALSE], is.numeric, logical(1)))
  )
  x <- as.matrix(data[, feat_cols, drop = FALSE])
  positive_class <- as.character(positive_class)
  lev <- levels(factor(data$Group))
  stopifnot(length(lev) == 2L, positive_class %in% lev,
            !any(is.infinite(x)))
  y <- factor(data$Group, levels = c(setdiff(lev, positive_class), positive_class))
  ID <- if ("ID" %in% names(data)) as.character(data$ID) else rownames(data)
  stopifnot(!anyNA(ID), all(nzchar(ID)), !anyDuplicated(ID), min(table(y)) >= 5L)

  fractions <- unique(train_fraction_grid)
  alphas <- unique(alpha_fixed_grid)
  minima <- unique(lambda_ratio_min_grid)
  tops <- unique(as.integer(TOP_n))
  parameter_grid <- expand.grid(
    train_fraction = fractions, TOP_n = tops,
    alpha_fixed = alphas, lambda_ratio_min = minima,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  parameter_grid <- cbind(Parameter = seq_len(nrow(parameter_grid)), parameter_grid)
  n_holdouts <- length(fractions) * n_repeats
  stopifnot(seed + 4000 + n_holdouts <= .Machine$integer.max)
  n_evaluations <- nrow(parameter_grid) * n_repeats
  if (verbose) message(sprintf(
    paste0("收到参数：train_fraction={%s}, TOP_n={%s}, alpha={%s}, ",
           "lambda_ratio_min={%s}；%d 组参数，%d 次外层拆分，计划 %d 次最终评估。"),
    paste(fractions, collapse = ","), paste(tops, collapse = ","),
    paste(alphas, collapse = ","), paste(minima, collapse = ","),
    nrow(parameter_grid), n_holdouts, n_evaluations
  ))
  if (verbose) message(sprintf(
    "统一重复次数 n_repeats = %d；每次外层拆分的训练集内进行 %d 次特征频率子抽样。",
    n_repeats, n_repeats
  ))

  stratified_sample <- function(z, fraction) {
    sort(unlist(lapply(levels(z), function(g) {
      ii <- which(z == g)
      ii[sample.int(length(ii), floor(length(ii) * fraction))]
    }), use.names = FALSE))
  }
  # 先生成外层拆分，避免内部拟合或 CV 的随机数消耗改变下一次拆分。
  split_plans <- lapply(seq_along(fractions), function(f) {
    set.seed(as.integer(seed + f - 1L))
    lapply(seq_len(n_repeats), function(i) stratified_sample(y, fractions[f]))
  })

  # 随机拆分后仅根据 ID、分组和已计划的验证次数修补遗漏。
  # 保持每次每类的训练/验证样本数不变；不参考表达量或拟合结果。
  split_adjustments <- data.frame(
    Holdout = character(), train_fraction = numeric(), Repeat = integer(),
    Removed_Validation_ID = character(), Added_Validation_ID = character(),
    stringsAsFactors = FALSE
  )
  if (ensure_oof_coverage) {
    for (f in seq_along(fractions)) {
      va_plans <- lapply(split_plans[[f]], function(tr) setdiff(seq_along(y), tr))
      set.seed(as.integer(seed + 3000 + f))
      for (g in levels(y)) {
        ii <- which(y == g)
        n_val_g <- length(ii) - floor(length(ii) * fractions[f])
        if (n_val_g < 2L || n_repeats * n_val_g < length(ii)) {
          stop(sprintf(
            "分组 %s 无法覆盖全部样本：每次验证 %d 个，至少需要 %d 次外层重复；验证每类还需至少 2 个样本。",
            g, n_val_g, ceiling(length(ii) / max(1L, n_val_g))
          ), call. = FALSE)
        }
        counts <- tabulate(unlist(va_plans, use.names = FALSE), nbins = length(y))
        uncovered <- ii[counts[ii] == 0L]
        for (u in uncovered) {
          donor_repeats <- which(vapply(va_plans, function(va) {
            any(y[va] == g & counts[va] > 1L)
          }, logical(1)))
          stopifnot(length(donor_repeats) > 0L)
          r <- donor_repeats[sample.int(length(donor_repeats), 1L)]
          donors <- va_plans[[r]][y[va_plans[[r]]] == g &
                                   counts[va_plans[[r]]] > 1L]
          donor <- donors[sample.int(length(donors), 1L)]
          va_plans[[r]][match(donor, va_plans[[r]])] <- u
          va_plans[[r]] <- sort(va_plans[[r]])
          counts[donor] <- counts[donor] - 1L
          counts[u] <- counts[u] + 1L
          split_adjustments <- rbind(split_adjustments, data.frame(
            Holdout = sprintf("H%04d", (f - 1L) * n_repeats + r),
            train_fraction = fractions[f], Repeat = r,
            Removed_Validation_ID = ID[donor], Added_Validation_ID = ID[u],
            stringsAsFactors = FALSE
          ))
        }
      }
      stopifnot(all(tabulate(unlist(va_plans, use.names = FALSE),
                            nbins = length(y)) > 0L))
      split_plans[[f]] <- lapply(va_plans, function(va) setdiff(seq_along(y), va))
    }
    if (verbose && nrow(split_adjustments)) {
      message(sprintf("为覆盖全部样本，修补 %d 个验证名额；训练比例和外层重复次数保持不变。",
                      nrow(split_adjustments)))
    }
  }

  capture_conditions <- function(expr) {
    notes <- character()
    value <- tryCatch(withCallingHandlers(
      force(expr), warning = function(w) {
        notes <<- c(notes, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ), error = function(e) e)
    list(value = value, warning = paste(unique(notes), collapse = " | "))
  }
  bind_list <- function(z) {
    if (!length(z)) return(data.frame())
    out <- do.call(rbind, z)
    rownames(out) <- NULL
    out
  }

  # 不预先按 MAD、缺失率、方差或 topN 缩小检验集合。
  # 检验基于训练集原始观测值，先检验，再填补。
  differential_train <- function(raw, yt) {
    neg <- yt != positive_class
    pos <- yt == positive_class
    if (is.null(pvalue_cutoff)) {
      return(data.frame(
        Feature = colnames(raw),
        N_Negative = colSums(is.finite(raw[neg, , drop = FALSE])),
        N_Positive = colSums(is.finite(raw[pos, , drop = FALSE])),
        P = NA_real_, Warning = "", Reason = "", Testable = FALSE,
        Significant = TRUE, stringsAsFactors = FALSE
      ))
    }
    rows <- lapply(seq_len(ncol(raw)), function(j) {
      a <- raw[neg, j]
      b <- raw[pos, j]
      a <- a[is.finite(a)]
      b <- b[is.finite(b)]
      reason <- note <- ""
      p <- NA_real_
      if (length(a) < 2L || length(b) < 2L) {
        reason <- "至少一组的非缺失观测少于 2 个。"
      } else if (length(unique(c(a, b))) < 2L) {
        reason <- "所有观测值相同，无法检验。"
      } else {
        z <- capture_conditions(if (diff_method == "wilcox") {
          stats::wilcox.test(a, b, alternative = "two.sided",
                             paired = FALSE, exact = FALSE)$p.value
        } else {
          stats::t.test(a, b, alternative = "two.sided",
                        paired = FALSE, var.equal = FALSE)$p.value
        })
        note <- z$warning
        if (inherits(z$value, "error")) {
          reason <- conditionMessage(z$value)
        } else if (length(z$value) == 1L && is.finite(z$value) &&
                   z$value >= 0 && z$value <= 1) {
          p <- as.numeric(z$value)
        } else {
          reason <- "检验未返回有效 p 值。"
        }
      }
      data.frame(Feature = colnames(raw)[j], N_Negative = length(a),
                 N_Positive = length(b), P = p, Warning = note,
                 Reason = reason, stringsAsFactors = FALSE)
    })
    tab <- bind_list(rows)
    tab$Testable <- is.finite(tab$P)
    # 直接使用原始 p 值筛选；不可检验项不入选。
    tab$Significant <- tab$Testable & tab$P < pvalue_cutoff
    tab
  }

  # 差异分析之后的数值处理：仅去掉无法填补/标准化的列，不做特征排序。
  preprocess_train <- function(raw) {
    med <- apply(raw, 2, stats::median, na.rm = TRUE)
    z <- raw[, is.finite(med), drop = FALSE]
    med <- med[colnames(z)]
    if (ncol(z) < 2L) stop("可填补的差异变量少于 2 个。", call. = FALSE)
    for (j in seq_len(ncol(z))) z[is.na(z[, j]), j] <- med[j]
    sigma <- apply(z, 2, stats::sd)
    mu <- colMeans(z)
    keep <- is.finite(sigma) & sigma > 0 & is.finite(mu)
    z <- z[, keep, drop = FALSE]
    if (ncol(z) < 2L) stop("可标准化的差异变量少于 2 个。", call. = FALSE)
    med <- med[colnames(z)]
    mu <- mu[colnames(z)]
    sigma <- sigma[colnames(z)]
    z <- sweep(sweep(z, 2, mu, "-"), 2, sigma, "/")
    stopifnot(all(is.finite(z)))
    list(x = z, median = med, mu = mu, sigma = sigma)
  }
  transform_new <- function(raw, prep, feats) {
    z <- raw[, feats, drop = FALSE]
    for (j in seq_along(feats)) z[is.na(z[, j]), j] <- prep$median[feats[j]]
    z <- sweep(sweep(z, 2, prep$mu[feats], "-"), 2, prep$sigma[feats], "/")
    stopifnot(all(is.finite(z)))
    z
  }
  make_lambda <- function(xs, ys, alpha, ratio_min) {
    lambda_max <- max(abs(drop(crossprod(xs, ys - mean(ys))))) /
      (nrow(xs) * alpha)
    if (!is.finite(lambda_max) || lambda_max <= 0) {
      stop("lambda_max 无效，无法建立该路径。", call. = FALSE)
    }
    ratio <- exp(seq(0, log(ratio_min), length.out = n_lambda))
    list(lambda_max = lambda_max, lambda = lambda_max * ratio)
  }
  validate_path <- function(fit) {
    if (length(fit$jerr) != 1L || !is.finite(fit$jerr) || fit$jerr > 0L) {
      stop("glmnet 返回致命或无效错误码。", call. = FALSE)
    }
    beta <- as.matrix(fit$beta)
    if (ncol(beta) != length(fit$lambda) || length(fit$a0) != ncol(beta)) {
      stop("系数矩阵与返回的 lambda 不匹配。", call. = FALSE)
    }
    valid <- is.finite(fit$lambda) & fit$lambda > 0 &
      is.finite(fit$a0) & colSums(!is.finite(beta)) == 0L
    if (!any(valid)) stop("路径没有有效模型。", call. = FALSE)
    list(beta = beta, valid = valid)
  }

  diagnostic_rows <- final_diagnostic_rows <- list()
  # 每个下限单独统计；所有 alpha 先分别计数，之后可池化或分开排序。
  screen_frequency <- function(raw_train, yt, candidates, inner_plans,
                               ratio_min, key, fraction, repeat_id) {
    counts <- matrix(0, nrow = length(candidates), ncol = length(alphas),
                     dimnames = list(candidates, as.character(alphas)))
    totals <- partial <- failed <- numeric(length(alphas))
    for (i in seq_len(n_repeats)) {
      sub <- inner_plans[[i]]
      # 同一子抽样的所有 alpha 共用相同预处理，全部候选变量进入该步骤。
      prep_result <- capture_conditions(preprocess_train(
        raw_train[sub, candidates, drop = FALSE]
      ))
      ys <- as.integer(yt[sub] == positive_class)
      for (a in seq_along(alphas)) {
        diag <- data.frame(
          Holdout = key, train_fraction = fraction, Repeat = repeat_id,
          Inner_Repeat = i, alpha_fixed = alphas[a], lambda_ratio_min = ratio_min,
          N_Candidates = length(candidates), N_Usable = NA_integer_,
          Requested = n_lambda, Returned = 0L, Used = 0L,
          Jerr = NA_integer_, Status = "failed", Warning = prep_result$warning,
          Error = "", stringsAsFactors = FALSE
        )
        z <- capture_conditions({
          if (inherits(prep_result$value, "error")) stop(prep_result$value)
          xs <- prep_result$value$x
          grid <- make_lambda(xs, ys, alphas[a], ratio_min)
          fit <- glmnet::glmnet(
            xs, ys, family = "binomial", alpha = alphas[a],
            lambda = grid$lambda, standardize = FALSE,
            thresh = 1e-8, maxit = as.integer(fit_maxit)
          )
          checked <- validate_path(fit)
          stopifnot(identical(rownames(checked$beta), colnames(xs)))
          list(fit = fit, checked = checked, feats = colnames(xs))
        })
        diag$Warning <- paste(c(prep_result$warning, z$warning)[
          nzchar(c(prep_result$warning, z$warning))], collapse = " | ")
        if (!inherits(prep_result$value, "error")) {
          diag$N_Usable <- ncol(prep_result$value$x)
        }
        if (inherits(z$value, "error")) {
          diag$Error <- conditionMessage(z$value)
          failed[a] <- failed[a] + 1L
        } else {
          ans <- z$value
          fit <- ans$fit
          valid <- ans$checked$valid
          beta <- ans$checked$beta[, valid, drop = FALSE]
          diag$Returned <- length(fit$lambda)
          diag$Used <- ncol(beta)
          diag$Jerr <- fit$jerr
          diag$Status <- if (fit$jerr == 0L && all(valid) &&
                             length(fit$lambda) == n_lambda) "complete" else "partial"
          if (diag$Status == "partial") partial[a] <- partial[a] + 1L
          # 子抽样中无法使用的候选变量，在这些模型中计为未入选。
          counts[ans$feats, a] <- counts[ans$feats, a] +
            rowSums(abs(beta) > coef_tol)
          totals[a] <- totals[a] + ncol(beta)
        }
        diagnostic_rows[[length(diagnostic_rows) + 1L]] <<- diag
      }
    }
    list(counts = counts, totals = totals, partial = partial, failed = failed)
  }

  auc_roc <- function(response, probability) {
    stopifnot(length(response) == length(probability), all(is.finite(probability)))
    roc <- pROC::roc(response, as.numeric(probability), levels = levels(y),
                     direction = "<", quiet = TRUE)
    list(roc = roc, auc = as.numeric(pROC::auc(roc)))
  }
  # TOP_n 个输入变量在完整外层训练集上拟合；CV 不接触外层验证集。
  fit_final <- function(prep, feats, yt, raw_val, yv, alpha, ratio_min, foldid,
                        expected_top_n) {
    stopifnot(length(feats) == expected_top_n, !anyDuplicated(feats))
    xt <- prep$x[, feats, drop = FALSE]
    stopifnot(ncol(xt) == expected_top_n, identical(colnames(xt), feats))
    yt01 <- as.integer(yt == positive_class)
    grid <- make_lambda(xt, yt01, alpha, ratio_min)
    fit <- glmnet::cv.glmnet(
      xt, yt01, family = "binomial", alpha = alpha,
      lambda = grid$lambda, foldid = foldid,
      type.measure = "deviance", grouped = FALSE, keep = TRUE,
      standardize = FALSE, thresh = 1e-8, maxit = as.integer(fit_maxit)
    )
    checked <- validate_path(fit$glmnet.fit)
    master_index <- match(fit$lambda, fit$glmnet.fit$lambda)
    matched <- !is.na(master_index)
    master_valid <- rep(FALSE, length(fit$lambda))
    master_valid[matched] <- checked$valid[master_index[matched]]
    # 仅在所有 CV 留出预测均有效的 lambda 上选参，避免比较不完整折结果。
    stopifnot(is.matrix(fit$fit.preval), ncol(fit$fit.preval) == length(fit$lambda))
    valid <- is.finite(fit$cvm) & master_valid &
      colSums(!is.finite(fit$fit.preval)) == 0L
    if (!any(valid)) stop("最终模型没有完整有效的训练集 CV 点。", call. = FALSE)
    best_deviance <- min(fit$cvm[valid])
    chosen <- max(fit$lambda[valid & fit$cvm == best_deviance])
    chosen_index <- match(chosen, fit$glmnet.fit$lambda)
    final_nz <- sum(abs(checked$beta[, chosen_index]) > coef_tol)
    # 只有模型和 lambda 确定之后，才处理外层验证集并计算 AUC。
    xv <- transform_new(raw_val, prep, feats)
    stopifnot(ncol(xv) == expected_top_n, identical(colnames(xv), colnames(xt)))
    pt <- as.numeric(stats::predict(fit, xt, s = chosen, type = "response"))
    pv <- as.numeric(stats::predict(fit, xv, s = chosen, type = "response"))
    trn <- auc_roc(yt, pt)
    val <- auc_roc(yv, pv)
    fit$fit.preval <- NULL
    list(fit = fit, final_lambda = chosen,
         final_ratio = chosen / grid$lambda_max, cv_deviance = best_deviance,
         n_nonzero = final_nz, pt = pt, pv = pv,
         trn = trn, val = val, jerr = fit$glmnet.fit$jerr)
  }

  result_rows <- frequency_tables <- selected_tables <- differential_tables <- list()
  models <- roc_examples <- predictions <- splits <- list()
  h_id <- 0L
  for (f in seq_along(fractions)) {
    fraction <- fractions[f]
    params_f <- parameter_grid[parameter_grid$train_fraction == fraction, , drop = FALSE]
    for (i in seq_len(n_repeats)) {
      h_id <- h_id + 1L
      key <- sprintf("H%04d", h_id)
      tr <- split_plans[[f]][[i]]
      va <- setdiff(seq_along(y), tr)
      yt <- y[tr]
      raw_train <- x[tr, , drop = FALSE]
      splits[[key]] <- list(train_fraction = fraction, Repeat = i, tr = tr, va = va,
                           train_ID = ID[tr], validation_ID = ID[va])
      if (verbose) message(sprintf("外层拆分 [%d/%d] train_fraction=%.2f, repeat=%d",
                                   h_id, n_holdouts, fraction, i))

      outer <- capture_conditions({
        if (min(table(yt)) < 3L || min(table(y[va])) < 2L) {
          stop("该比例不能保证训练每类至少 3 个、验证每类至少 2 个样本。", call. = FALSE)
        }
        if (min(floor(as.numeric(table(yt)) * inner_fraction)) < 2L) {
          stop("训练子抽样每类至少需要 2 个样本，请增大 inner_fraction。", call. = FALSE)
        }
        differential <- differential_train(raw_train, yt)
        sig <- differential$Feature[differential$Significant]
        if (keep_differential) differential_tables[[length(differential_tables) + 1L]] <-
          cbind(Holdout = key, train_fraction = fraction, Repeat = i, differential)
        if (length(sig) < 2L) stop(
          if (is.null(pvalue_cutoff)) "训练集候选变量少于 2 个。" else
            sprintf("训练集 p value < %.3g 的差异变量少于 2 个。", pvalue_cutoff),
          call. = FALSE
        )
        prep <- preprocess_train(raw_train[, sig, drop = FALSE])
        candidates <- colnames(prep$x)
        set.seed(as.integer(seed + 1000 + h_id))
        inner_plans <- lapply(seq_len(n_repeats),
                              function(j) stratified_sample(yt, inner_fraction))
        # 所有 TOP_n/alpha/lambda 下限使用相同 CV 折划分。
        set.seed(as.integer(seed + 2000 + h_id))
        folds <- as.integer(min(cv_nfolds, min(table(yt))))
        foldid <- integer(length(yt))
        for (g in levels(yt)) {
          ii <- which(yt == g)
          foldid[ii[sample.int(length(ii))]] <- rep(seq_len(folds), length.out = length(ii))
        }
        list(differential = differential, prep = prep, candidates = candidates,
             inner_plans = inner_plans, foldid = foldid, n_diff = length(sig))
      })
      outer_error <- if (inherits(outer$value, "error")) conditionMessage(outer$value) else ""

      selection <- list()
      if (!nzchar(outer_error)) {
        o <- outer$value
        for (m in seq_along(minima)) {
          screened <- screen_frequency(raw_train, yt, o$candidates, o$inner_plans,
                                        minima[m], key, fraction, i)
          scopes <- if (frequency_alpha_mode == "pooled") list(seq_along(alphas)) else
            lapply(seq_along(alphas), function(a) a)
          selection[[m]] <- vector("list", length(scopes))
          for (s in seq_along(scopes)) {
            aa <- scopes[[s]]
            n_models <- sum(screened$totals[aa])
            count <- rowSums(screened$counts[, aa, drop = FALSE])
            stopifnot(all(count <= n_models))
            tab <- data.frame(
              Feature = names(count), Count = as.numeric(count),
              Frequency = if (n_models > 0) as.numeric(count) / n_models else NA_real_,
              Pvalue = o$differential$P[match(names(count), o$differential$Feature)],
              stringsAsFactors = FALSE
            )
            tab <- tab[order(-tab$Count, tab$Feature, method = "radix"), , drop = FALSE]
            tab$Rank <- seq_len(nrow(tab))
            rownames(tab) <- NULL
            selection[[m]][[s]] <- list(
              table = tab, n_models = n_models,
              n_requested = n_repeats * length(aa) * n_lambda,
              n_partial = sum(screened$partial[aa]), n_failed = sum(screened$failed[aa])
            )
            frequency_tables[[length(frequency_tables) + 1L]] <- cbind(
              Holdout = key, train_fraction = fraction, Repeat = i,
              lambda_ratio_min = minima[m],
              Selection_Alpha = paste(alphas[aa], collapse = ";"),
              N_Models = n_models, tab
            )
          }
        }
      }

      for (p in seq_len(nrow(params_f))) {
        pars <- params_f[p, , drop = FALSE]
        row <- cbind(pars, data.frame(
          Holdout = key, Repeat = i, N_TRN = length(tr), N_VAL = length(va),
          N_Diff = NA_integer_, N_Usable_Diff = NA_integer_,
          N_PositiveFrequency = NA_integer_, N_Used = NA_integer_,
          N_Screen_Models = NA_real_, N_Requested_Models = NA_real_,
          N_PartialPaths = NA_integer_, N_FailedPaths = NA_integer_,
          Final_Lambda = NA_real_, Final_Lambda_Ratio = NA_real_,
          CV_Deviance = NA_real_, N_Final_Nonzero = NA_integer_,
          `TRN-AUC` = NA_real_, `VAL-AUC` = NA_real_,
          Top_Features = "", Warning = outer$warning, Error = outer_error,
          check.names = FALSE, stringsAsFactors = FALSE
        ))
        if (!nzchar(outer_error)) {
          m <- match(pars$lambda_ratio_min, minima)
          a <- match(pars$alpha_fixed, alphas)
          s <- if (frequency_alpha_mode == "pooled") 1L else a
          sel <- selection[[m]][[s]]
          tab <- sel$table
          eligible <- tab$Feature[tab$Count > 0]
          row$N_Diff <- o$n_diff
          row$N_Usable_Diff <- length(o$candidates)
          row$N_PositiveFrequency <- length(eligible)
          row$N_Screen_Models <- sel$n_models
          row$N_Requested_Models <- sel$n_requested
          row$N_PartialPaths <- sel$n_partial
          row$N_FailedPaths <- sel$n_failed
          if (sel$n_models == 0) {
            row$Error <- "该频率统计范围内没有有效模型。"
          } else if (length(eligible) < pars$TOP_n) {
            row$Error <- sprintf("正频率差异变量仅 %d 个，少于 TOP_n=%d。",
                                 length(eligible), pars$TOP_n)
          } else {
            feats <- head(eligible, pars$TOP_n)
            row$N_Used <- length(feats)
            row$Top_Features <- paste(feats, collapse = ";")
            selected_tables[[length(selected_tables) + 1L]] <- cbind(
              pars, Holdout = key, Repeat = i,
              tab[match(feats, tab$Feature), , drop = FALSE]
            )
            z <- capture_conditions(fit_final(
              o$prep, feats, yt, x[va, , drop = FALSE], y[va],
              pars$alpha_fixed, pars$lambda_ratio_min, o$foldid,
              expected_top_n = pars$TOP_n
            ))
            row$Warning <- paste(c(row$Warning, z$warning)[
              nzchar(c(row$Warning, z$warning))], collapse = " | ")
            final_diagnostic_rows[[length(final_diagnostic_rows) + 1L]] <- data.frame(
              Parameter = pars$Parameter, Holdout = key, Repeat = i,
              Warning = z$warning,
              Jerr = if (inherits(z$value, "error")) NA_integer_ else z$value$jerr,
              Error = if (inherits(z$value, "error")) conditionMessage(z$value) else "",
              stringsAsFactors = FALSE
            )
            if (inherits(z$value, "error")) {
              row$Error <- conditionMessage(z$value)
            } else {
              ans <- z$value
              row$Final_Lambda <- ans$final_lambda
              row$Final_Lambda_Ratio <- ans$final_ratio
              row$CV_Deviance <- ans$cv_deviance
              row$N_Final_Nonzero <- ans$n_nonzero
              row[["TRN-AUC"]] <- ans$trn$auc
              row[["VAL-AUC"]] <- ans$val$auc
              # 正式分析可保存全部重复；默认仍只保存 plot_iter，以控制内存占用。
              if (keep_all_repeats || i == plot_iter) {
                model_key <- sprintf("P%03d_%s", pars$Parameter, key)
                models[[model_key]] <- list(
                  Parameter = pars$Parameter, Holdout = key, Repeat = i,
                  fit = ans$fit, feats = feats, final_lambda = ans$final_lambda,
                  preprocess = list(median = o$prep$median[feats],
                                    mu = o$prep$mu[feats], sigma = o$prep$sigma[feats]),
                  levels = levels(y), positive_class = positive_class,
                  train_ID = ID[tr], validation_ID = ID[va], foldid = o$foldid
                )
                roc_examples[[model_key]] <- list(TRN = ans$trn$roc, VAL = ans$val$roc)
                predictions[[model_key]] <- rbind(
                  data.frame(Parameter = pars$Parameter, Holdout = key, Repeat = i,
                             ID = ID[tr], Set = "TRN", Group = as.character(yt),
                             Probability = ans$pt),
                  data.frame(Parameter = pars$Parameter, Holdout = key, Repeat = i,
                             ID = ID[va], Set = "VAL", Group = as.character(y[va]),
                             Probability = ans$pv)
                )
              }
            }
          }
        }
        evaluation_id <- length(result_rows) + 1L
        if (verbose) {
          prefix <- sprintf(
            paste0("  最终训练/验证 [%d/%d] train_fraction=%.2f, repeat=%d, ",
                   "TOP_n=%d, alpha=%.3g, lambda_ratio_min=%.3g"),
            evaluation_id, n_evaluations, fraction, i,
            pars$TOP_n, pars$alpha_fixed, pars$lambda_ratio_min
          )
          if (nzchar(row$Error)) {
            message(prefix, "；失败：", row$Error)
          } else {
            message(prefix, sprintf(
              "；实际输入=%d，最终非零=%d，TRN-AUC=%.4f，VAL-AUC=%.4f",
              row$N_Used, row$N_Final_Nonzero, row[["TRN-AUC"]], row[["VAL-AUC"]]
            ))
          }
        }
        result_rows[[evaluation_id]] <- row
      }
      if (verbose && nzchar(outer_error)) message("  本次拆分的全部组合记为失败：", outer_error)
    }
  }
  holdout_results <- bind_list(result_rows)
  # 成功记录必须分别使用各自指定的 TOP_n 个输入变量；非零系数可以更少。
  successful <- holdout_results$Error == "" &
    is.finite(holdout_results[["TRN-AUC"]]) & is.finite(holdout_results[["VAL-AUC"]])
  stopifnot(
    all(holdout_results$N_Used[successful] == holdout_results$TOP_n[successful]),
    all(holdout_results$N_Final_Nonzero[successful] <= holdout_results$N_Used[successful])
  )
  summarize_auc <- function(z) {
    if (!length(z)) return(rep(NA_real_, 3L))
    c(mean(z), unname(stats::quantile(z, 0.025)), unname(stats::quantile(z, 0.975)))
  }
  indexfilterresult <- bind_list(lapply(seq_len(nrow(parameter_grid)), function(p) {
    z <- holdout_results[holdout_results$Parameter == p, , drop = FALSE]
    ok <- z$Error == "" & is.finite(z[["TRN-AUC"]]) & is.finite(z[["VAL-AUC"]])
    trn <- summarize_auc(z[["TRN-AUC"]][ok])
    val <- summarize_auc(z[["VAL-AUC"]][ok])
    cbind(parameter_grid[p, , drop = FALSE], data.frame(
      `TRN-AUC` = trn[1], `TRN-AUClower` = trn[2], `TRN-AUCupper` = trn[3],
      `VAL-AUC` = val[1], `VAL-AUClower` = val[2], `VAL-AUCupper` = val[3],
      `VAL-AUCwidth` = val[3] - val[2],
      n_success = sum(ok), n_failed = nrow(z) - sum(ok),
      check.names = FALSE
    ))
  }))
  indexfilterresult <- indexfilterresult[
    order(-indexfilterresult[["VAL-AUC"]], -indexfilterresult[["VAL-AUClower"]],
          indexfilterresult$TOP_n, indexfilterresult$Parameter, na.last = TRUE), , drop = FALSE
  ]
  rownames(indexfilterresult) <- NULL
  stopifnot(nrow(holdout_results) == nrow(parameter_grid) * n_repeats,
            all(indexfilterresult$n_success + indexfilterresult$n_failed == n_repeats))
  screen_diagnostics <- bind_list(diagnostic_rows)
  if (!nrow(screen_diagnostics)) screen_diagnostics <- data.frame(
    Holdout = character(), train_fraction = numeric(), Repeat = integer(),
    Inner_Repeat = integer(), alpha_fixed = numeric(), lambda_ratio_min = numeric(),
    N_Candidates = integer(), N_Usable = integer(), Requested = integer(),
    Returned = integer(), Used = integer(), Jerr = integer(), Status = character(),
    Warning = character(), Error = character(), stringsAsFactors = FALSE
  )
  final_diagnostics <- bind_list(final_diagnostic_rows)
  if (!nrow(final_diagnostics)) final_diagnostics <- data.frame(
    Parameter = integer(), Holdout = character(), Repeat = integer(),
    Warning = character(), Jerr = integer(), Error = character(), stringsAsFactors = FALSE
  )
  list(
    indexfilterresult = indexfilterresult,
    holdout_results = holdout_results,
    feature_frequency = bind_list(frequency_tables),
    selected_features = bind_list(selected_tables),
    differential_results = bind_list(differential_tables),
    screen_diagnostics = screen_diagnostics,
    final_diagnostics = final_diagnostics,
    indexfilter_errors = holdout_results[holdout_results$Error != "", , drop = FALSE],
    models = models, roc_examples = roc_examples, predictions = bind_list(predictions),
    splits = splits, split_adjustments = split_adjustments,
    parameter_grid = parameter_grid,
    config = list(positive_class = positive_class,
                  train_fraction_grid = fractions, TOP_n = tops,
                  alpha_fixed_grid = alphas, lambda_ratio_min_grid = minima,
                  n_repeats = n_repeats,
                  n_lambda = n_lambda,
                  inner_fraction = inner_fraction, pvalue_cutoff = pvalue_cutoff,
                  diff_method = diff_method, frequency_alpha_mode = frequency_alpha_mode,
                  selection_alpha_grid = alphas, cv_nfolds = cv_nfolds,
                  plot_iter = plot_iter, seed = seed, coef_tol = coef_tol,
                  fit_maxit = fit_maxit, keep_differential = keep_differential,
                  keep_all_repeats = keep_all_repeats,
                  ensure_oof_coverage = ensure_oof_coverage)
  )
}

# 对保存的模型预测新样本。概率始终对应 positive_class。
predict_frequency_model <- function(model, newdata) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("请先安装 R 包 glmnet。")
  stopifnot(all(model$feats %in% colnames(newdata)),
            all(vapply(newdata[, model$feats, drop = FALSE], is.numeric, logical(1))))
  z <- as.matrix(newdata[, model$feats, drop = FALSE])
  stopifnot(!any(is.infinite(z)))
  for (j in seq_along(model$feats)) {
    z[is.na(z[, j]), j] <- model$preprocess$median[model$feats[j]]
  }
  z <- sweep(sweep(z, 2, model$preprocess$mu, "-"), 2, model$preprocess$sigma, "/")
  stopifnot(all(is.finite(z)))
  as.numeric(stats::predict(model$fit, z, s = model$final_lambda, type = "response"))
}

# 绘制 plot_iter 拆分中的训练/验证 ROC；Parameter 对应结果表中的稳定编号。
plot_frequency_roc <- function(result, parameter_id = 1L) {
  eligible <- vapply(result$models, function(m) {
    m$Parameter == parameter_id && m$Repeat == result$config$plot_iter
  }, logical(1))
  keys <- names(result$models)[eligible]
  if (!length(keys)) stop("指定参数在 plot_iter 拆分中没有成功且已保存的模型。", call. = FALSE)
  r <- result$roc_examples[[keys[1L]]]
  graphics::plot(r$TRN, col = "#2878B5", lwd = 2, legacy.axes = TRUE,
                 main = sprintf("Parameter %d, repeat %d", parameter_id, result$config$plot_iter))
  graphics::plot(r$VAL, col = "#D95319", lwd = 2, add = TRUE)
  graphics::legend("bottomright", bty = "n", lwd = 2, col = c("#2878B5", "#D95319"),
                   legend = c(sprintf("TRN AUC = %.3f", as.numeric(pROC::auc(r$TRN))),
                              sprintf("VAL AUC = %.3f", as.numeric(pROC::auc(r$VAL)))))
  invisible(r)
}


frequency_result_lists <- function(result, fpr_grid = seq(0, 1, length.out = 101L)) {
  if (!isTRUE(result$config$keep_all_repeats)) {
    stop("请以 keep_all_repeats=TRUE 重新运行正式分析，再生成全部重复的明细。",
         call. = FALSE)
  }
  stopifnot(is.numeric(fpr_grid), length(fpr_grid) >= 2L,
            all(is.finite(fpr_grid)), all(fpr_grid >= 0 & fpr_grid <= 1),
            all(diff(fpr_grid) > 0))
  add_parameters <- function(z) {
    cols <- c("train_fraction", "TOP_n", "alpha_fixed", "lambda_ratio_min")
    ii <- match(z$Parameter, result$parameter_grid$Parameter)
    stopifnot(!anyNA(ii))
    out <- cbind(z, result$parameter_grid[ii, cols, drop = FALSE])
    rownames(out) <- NULL
    out
  }

  pred <- result$predictions
  if (!nrow(pred)) {
    pred <- data.frame(Parameter = integer(), Holdout = character(),
                       Repeat = integer(), ID = character(), Set = character(),
                       Group = character(), Probability = numeric())
  }
  results_list <- list(add_parameters(pred))
  # 包含成功和失败拆分；失败的 AUC 为 NA，原因在 Error 列。
  auc_list <- list(result$holdout_results)

  model_list <- list(add_parameters(data.frame(
    Parameter = integer(), Holdout = character(), Repeat = integer(),
    Feature = character(), Coefficient = numeric(), Final_Lambda = numeric()
  )))
  curve_list <- list(add_parameters(data.frame(
    Parameter = integer(), Holdout = character(), Repeat = integer(),
    Set = character(), FPR = numeric(), TPR = numeric()
  )))
  for (key in names(result$models)) {
    m <- result$models[[key]]
    b <- as.matrix(stats::coef(m$fit, s = m$final_lambda))
    # 包含截距和零系数；变量系数对应训练集标准化后的尺度。
    model_list[[length(model_list) + 1L]] <- add_parameters(data.frame(
      Parameter = m$Parameter, Holdout = m$Holdout, Repeat = m$Repeat,
      Feature = rownames(b), Coefficient = as.numeric(b[, 1L]),
      Final_Lambda = m$final_lambda, stringsAsFactors = FALSE
    ))
    rr <- result$roc_examples[[key]]
    stopifnot(!is.null(rr$TRN), !is.null(rr$VAL))
    for (set_name in c("TRN", "VAL")) {
      r <- rr[[set_name]]
      # 统一 FPR 网格只用于曲线展示；AUC 仍来自原始预测的 pROC 计算。
      tpr <- stats::approx(
        x = 1 - r$specificities, y = r$sensitivities,
        xout = fpr_grid, ties = max, rule = 2
      )$y
      curve_list[[length(curve_list) + 1L]] <- add_parameters(data.frame(
        Parameter = m$Parameter, Holdout = m$Holdout, Repeat = m$Repeat,
        Set = set_name, FPR = fpr_grid, TPR = tpr, stringsAsFactors = FALSE
      ))
    }
  }
  list(results_list = results_list, model_list = model_list,
       auc_list = auc_list, curve_list = curve_list)
}



# ============================ 输入与运行入口 ============================
.frequency_bind_rows <- function(z) {
  out <- do.call(rbind, z)
  if (is.null(out)) return(data.frame())
  rownames(out) <- NULL
  out
}

read_frequency_input <- function(input_file) {
  if (!is.character(input_file) || length(input_file) != 1L ||
      is.na(input_file) || !file.exists(input_file)) {
    stop("找不到 input.csv。请把文件放在 getwd() 所示目录，或设置 input_file 为完整路径。",
         call. = FALSE)
  }
  d <- utils::read.csv(
    input_file, header = TRUE, row.names = NULL, check.names = FALSE,
    stringsAsFactors = FALSE, colClasses = "character", fill = FALSE,
    na.strings = c("", "NA", "NaN"), fileEncoding = "UTF-8",
    comment.char = ""
  )
  if (nrow(d) < 1L || ncol(d) < 4L) {
    stop("CSV 需要至少 1 行样本和 4 列：ID、Group、至少 2 个特征。",
         call. = FALSE)
  }
  names(d)[1L] <- sub(intToUtf8(65279L), "", names(d)[1L], fixed = TRUE)
  if (anyNA(names(d)) || any(!nzchar(trimws(names(d)))) ||
      anyDuplicated(names(d))) {
    stop("CSV 列名必须非空且不重复。", call. = FALSE)
  }
  names(d)[1:2] <- c("ID", "Group")
  if (any(c("ID", "Group") %in% names(d)[-(1:2)])) {
    stop("特征列不能命名为 ID 或 Group。", call. = FALSE)
  }
  d$ID <- trimws(d$ID)
  d$Group <- trimws(d$Group)
  if (anyNA(d$ID) || any(!nzchar(d$ID)) || anyDuplicated(d$ID)) {
    stop("第 1 列样本 ID 必须非空且唯一；同一样本不能出现多行。",
         call. = FALSE)
  }
  if (anyNA(d$Group) || !setequal(unique(d$Group), c("Tumor", "Control"))) {
    stop("第 2 列分组必须同时包含 Tumor 和 Control，区分大小写且不能缺失。",
         call. = FALSE)
  }
  for (j in seq.int(3L, ncol(d))) {
    raw <- trimws(d[[j]])
    missing <- is.na(raw) | raw %in% c("", "NA", "NaN")
    value <- suppressWarnings(as.numeric(raw))
    bad <- !missing & (is.na(value) | !is.finite(value))
    if (any(bad)) {
      k <- which(bad)[1L]
      stop(sprintf("特征 %s 在样本 %s 中含非数值或无穷大：%s",
                   names(d)[j], d$ID[k], raw[k]), call. = FALSE)
    }
    value[missing] <- NA_real_
    d[[j]] <- value
  }
  d
}

.frequency_require_packages <- function(install_missing = TRUE) {
  pkgs <- c("glmnet", "pROC")
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) && install_missing) {
    message("正在安装缺少的 R 包：", paste(missing, collapse = ", "))
    attempt <- tryCatch({
      utils::install.packages(missing, repos = "https://cloud.r-project.org",
                              dependencies = NA)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(attempt)) message("自动安装未完成：", attempt)
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) {
    stop("缺少 R 包：", paste(missing, collapse = ", "),
         "。请先 install.packages(c(\"glmnet\", \"pROC\"))，再运行脚本。",
         call. = FALSE)
  }
}

.frequency_validate_config <- function(cfg, d) {
  count_fields <- c("TOP_n", "n_repeats", "n_lambda",
                    "cv_nfolds", "plot_iter", "seed", "fit_maxit")
  for (k in count_fields) {
    v <- cfg[[k]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
        v != floor(v) || v > .Machine$integer.max || v < 0) {
      stop("配置 ", k, " 必须是单个非负整数。", call. = FALSE)
    }
  }
  for (k in c("ensure_oof_coverage", "require_complete_oof",
               "install_missing_packages", "save_models", "overwrite", "verbose")) {
    v <- cfg[[k]]
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("配置 ", k, " 必须是 TRUE 或 FALSE。", call. = FALSE)
    }
  }
  for (k in c("train_fraction", "alpha_fixed", "lambda_ratio_min", "inner_fraction")) {
    v <- cfg[[k]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0 ||
        (if (k == "alpha_fixed") v > 1 else v >= 1)) {
      stop("配置 ", k, " 超出允许范围。每次正式运行使用一个固定参数组合。",
           call. = FALSE)
    }
  }
  if (!identical(cfg$positive_class, "Tumor")) {
    stop("本入口的阳性类固定为 Tumor，以保证 Tumor_Probability 的含义一致。",
         call. = FALSE)
  }
  if (cfg$TOP_n < 2L || cfg$TOP_n > ncol(d) - 2L ||
      cfg$n_repeats < 2L ||
      cfg$n_lambda < 2L || cfg$cv_nfolds < 3L || cfg$fit_maxit < 1L ||
      cfg$plot_iter < 1L || cfg$plot_iter > cfg$n_repeats) {
    stop("检查整数配置：TOP_n 不超过特征数，重复及 lambda 点数至少为 2，CV 折数至少为 3。",
         call. = FALSE)
  }
  if (!is.null(cfg$pvalue_cutoff) &&
      (!is.numeric(cfg$pvalue_cutoff) || length(cfg$pvalue_cutoff) != 1L ||
       !is.finite(cfg$pvalue_cutoff) || cfg$pvalue_cutoff <= 0 || cfg$pvalue_cutoff >= 1)) {
    stop("pvalue_cutoff 应为 (0,1) 内的数值或 NULL。", call. = FALSE)
  }
  if (!is.numeric(cfg$coef_tol) || length(cfg$coef_tol) != 1L ||
      !is.finite(cfg$coef_tol) || cfg$coef_tol < 0) {
    stop("coef_tol 必须是单个非负有限数。", call. = FALSE)
  }
  if (!is.character(cfg$output_dir) || length(cfg$output_dir) != 1L ||
      is.na(cfg$output_dir) || !nzchar(cfg$output_dir)) {
    stop("output_dir 必须是一个非空目录路径。", call. = FALSE)
  }
  if (!is.character(cfg$diff_method) || length(cfg$diff_method) != 1L ||
      !cfg$diff_method %in% c("wilcox", "welch") ||
      !is.character(cfg$frequency_alpha_mode) ||
      length(cfg$frequency_alpha_mode) != 1L ||
      !cfg$frequency_alpha_mode %in% c("pooled", "by_alpha")) {
    stop("diff_method 使用 wilcox/welch；frequency_alpha_mode 使用 pooled/by_alpha。",
         call. = FALSE)
  }
  if (!is.numeric(cfg$fpr_grid) || length(cfg$fpr_grid) < 2L ||
      any(!is.finite(cfg$fpr_grid)) || any(cfg$fpr_grid < 0 | cfg$fpr_grid > 1) ||
      any(diff(cfg$fpr_grid) <= 0)) {
    stop("fpr_grid 必须是在 [0,1] 内严格递增的数值网格。", call. = FALSE)
  }
  n <- as.numeric(table(factor(d$Group, levels = c("Control", "Tumor"))))
  tr <- floor(n * cfg$train_fraction)
  va <- n - tr
  if (min(n) < 5L || min(tr) < 3L || min(va) < 2L ||
      min(floor(tr * cfg$inner_fraction)) < 2L) {
    stop("样本量/比例不足：每类总数至少 5，外层训练至少 3，验证至少 2，内层子抽样至少 2。",
         call. = FALSE)
  }
  if (cfg$ensure_oof_coverage && any(cfg$n_repeats * va < n)) {
    stop(sprintf("覆盖全部样本至少需要 %d 次外层重复，请增大 n_repeats。",
                 max(ceiling(n / va))), call. = FALSE)
  }
}

# config 只需指定要覆盖的项；NULL 会保留，故可关闭 p 值筛选。
run_frequency_pipeline <- function(config = list()) {
  defaults <- frequency_default_config()
  if (!is.list(config) || is.null(names(config)) && length(config) ||
      anyDuplicated(names(config))) {
    stop("config 必须是命名 list，参数名不能重复。", call. = FALSE)
  }
  unknown <- setdiff(names(config), names(defaults))
  if (length(unknown)) {
    stop("未知配置名：", paste(unknown, collapse = ", "), call. = FALSE)
  }
  cfg <- utils::modifyList(defaults, config, keep.null = TRUE)
  d <- read_frequency_input(cfg$input_file)
  .frequency_validate_config(cfg, d)

  files <- c(
    "oof_predictions.csv", "oof_predictions_summary.csv",
    "results_df.csv", "model_df.csv", "auc_df.csv", "curve_df.csv",
    "indexfilterresult.csv", "selected_features.csv", "feature_frequency.csv",
    "differential_results.csv", "screen_diagnostics.csv", "final_diagnostics.csv",
    "indexfilter_errors.csv", "split_membership.csv", "split_adjustments.csv",
    "run_config.R", "session_info.txt", "run_status.txt", "run_log.txt"
  )
  if (cfg$save_models) files <- c(files, "frequency_analysis.rds")
  if (!dir.exists(cfg$output_dir) &&
      !dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("无法创建输出目录：", cfg$output_dir, call. = FALSE)
  }
  cfg$input_file <- normalizePath(cfg$input_file, winslash = "/", mustWork = TRUE)
  cfg$output_dir <- normalizePath(cfg$output_dir, winslash = "/", mustWork = TRUE)
  paths <- file.path(cfg$output_dir, files)
  if (cfg$input_file %in% paths) stop("输入文件不能与输出文件使用同一路径。")
  existing <- paths[file.exists(paths)]
  if (length(existing) && !cfg$overwrite) {
    stop("输出目录已有同名结果。请更换 output_dir，或显式设置 overwrite=TRUE：",
         paste(basename(existing), collapse = ", "), call. = FALSE)
  }
  .frequency_require_packages(cfg$install_missing_packages)
  old_digits <- options(digits = 17)
  on.exit(options(old_digits), add = TRUE)
  log_path <- file.path(cfg$output_dir, "run_log.txt")
  writeLines(character(), log_path, useBytes = TRUE)
  utils::capture.output(dput(cfg), file = file.path(cfg$output_dir, "run_config.R"))
  write_status <- function(lines) {
    writeLines(enc2utf8(lines), file.path(cfg$output_dir, "run_status.txt"), useBytes = TRUE)
  }
  write_status(c("Status: RUNNING", paste("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
                 paste("Input_MD5:", unname(tools::md5sum(cfg$input_file)))))
  status_finalized <- FALSE
  on.exit({
    if (!status_finalized) try(write_status(c(
      "Status: ERROR", "Error: 分析或结果导出中断，请查看 R 控制台错误和 run_log.txt。"
    )), silent = TRUE)
  }, add = TRUE)
  on.exit(utils::capture.output(utils::sessionInfo(),
                               file = file.path(cfg$output_dir, "session_info.txt")), add = TRUE)
  analysis <- tryCatch(withCallingHandlers(
    run_frequency_holdouts(
      data = d, positive_class = cfg$positive_class,
      train_fraction_grid = cfg$train_fraction, TOP_n = cfg$TOP_n,
      alpha_fixed_grid = cfg$alpha_fixed, lambda_ratio_min_grid = cfg$lambda_ratio_min,
      n_repeats = cfg$n_repeats,
      n_lambda = cfg$n_lambda, pvalue_cutoff = cfg$pvalue_cutoff,
      inner_fraction = cfg$inner_fraction, diff_method = cfg$diff_method,
      frequency_alpha_mode = cfg$frequency_alpha_mode, cv_nfolds = cfg$cv_nfolds,
      plot_iter = cfg$plot_iter, seed = cfg$seed, coef_tol = cfg$coef_tol,
      fit_maxit = cfg$fit_maxit, keep_differential = TRUE,
      verbose = cfg$verbose, keep_all_repeats = TRUE,
      ensure_oof_coverage = cfg$ensure_oof_coverage
    ),
    message = function(m) cat(conditionMessage(m), file = log_path, append = TRUE)
  ), error = function(e) {
    write_status(c("Status: ERROR", paste("Error:", conditionMessage(e))))
    status_finalized <<- TRUE
    stop(e)
  })

  output_lists <- frequency_result_lists(analysis, fpr_grid = cfg$fpr_grid)
  bind_rows <- .frequency_bind_rows
  results_df <- bind_rows(output_lists$results_list)
  model_df <- bind_rows(output_lists$model_list)
  auc_df <- bind_rows(output_lists$auc_list)
  curve_df <- bind_rows(output_lists$curve_list)

  # 折外概率只来自外层 VAL，不使用 TRN 或内部 CV 的 fit.preval。
  oof <- results_df[results_df$Set == "VAL", , drop = FALSE]
  names(oof)[names(oof) == "Probability"] <- "Tumor_Probability"
  stopifnot(all(is.finite(oof$Tumor_Probability)),
            all(oof$Tumor_Probability >= 0 & oof$Tumor_Probability <= 1),
            !anyDuplicated(oof[, c("Parameter", "Holdout", "ID"), drop = FALSE]))
  for (key in names(analysis$models)) {
    m <- analysis$models[[key]]
    ids <- oof$ID[oof$Holdout == m$Holdout & oof$Parameter == m$Parameter]
    stopifnot(setequal(ids, m$validation_ID), !any(ids %in% m$train_ID))
  }
  by_id <- split(oof$Tumor_Probability, oof$ID)
  planned_counts <- table(factor(
    unlist(lapply(analysis$splits, function(s) s$validation_ID), use.names = FALSE),
    levels = d$ID
  ))
  oof_summary <- bind_rows(lapply(seq_len(nrow(d)), function(i) {
    p <- by_id[[d$ID[i]]]
    n <- length(p)
    data.frame(
      ID = d$ID[i], Group = d$Group[i],
      N_Planned_Validation = as.integer(planned_counts[d$ID[i]]),
      N_Validation_Predictions = n,
      Tumor_Probability_Mean = if (n) mean(p) else NA_real_,
      Tumor_Probability_Median = if (n) stats::median(p) else NA_real_,
      Tumor_Probability_SD = if (n > 1L) stats::sd(p) else NA_real_,
      Tumor_Probability_Min = if (n) min(p) else NA_real_,
      Tumor_Probability_Max = if (n) max(p) else NA_real_,
      Has_OOF = n > 0L, stringsAsFactors = FALSE
    )
  }))
  stopifnot(identical(oof_summary$ID, d$ID),
            all(oof_summary$N_Validation_Predictions <= oof_summary$N_Planned_Validation))
  split_membership <- bind_rows(lapply(names(analysis$splits), function(key) {
    s <- analysis$splits[[key]]
    rbind(
      data.frame(Holdout = key, Repeat = s$Repeat, Set = "TRN",
                 ID = s$train_ID, Group = d$Group[match(s$train_ID, d$ID)],
                 train_fraction = s$train_fraction),
      data.frame(Holdout = key, Repeat = s$Repeat, Set = "VAL",
                 ID = s$validation_ID, Group = d$Group[match(s$validation_ID, d$ID)],
                 train_fraction = s$train_fraction)
    )
  }))
  tables <- list(
    oof_predictions = oof, oof_predictions_summary = oof_summary,
    results_df = results_df, model_df = model_df, auc_df = auc_df, curve_df = curve_df,
    indexfilterresult = analysis$indexfilterresult,
    selected_features = analysis$selected_features,
    feature_frequency = analysis$feature_frequency,
    differential_results = analysis$differential_results,
    screen_diagnostics = analysis$screen_diagnostics,
    final_diagnostics = analysis$final_diagnostics,
    indexfilter_errors = analysis$indexfilter_errors,
    split_membership = split_membership, split_adjustments = analysis$split_adjustments
  )
  for (nm in names(tables)) {
    tab <- tables[[nm]]
    # 全部拆分失败时，频率/特征表可能没有列；仍保存可识别的空表。
    if (!ncol(tab)) tab <- data.frame(Note = character())
    utils::write.csv(tab, file.path(cfg$output_dir, paste0(nm, ".csv")),
                     row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  }
  missing_ids <- oof_summary$ID[!oof_summary$Has_OOF]
  n_success <- sum(!nzchar(auc_df$Error))
  bundle <- list(config = cfg, analysis = analysis,
                 results_df = results_df, model_df = model_df,
                 auc_df = auc_df, curve_df = curve_df,
                 oof_predictions = oof, oof_predictions_summary = oof_summary,
                 output_files = paths)
  if (cfg$save_models) saveRDS(bundle, file.path(cfg$output_dir, "frequency_analysis.rds"))
  write_status(c(
    if (length(missing_ids)) "Status: INCOMPLETE_OOF" else "Status: COMPLETE",
    paste("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Input_MD5:", unname(tools::md5sum(cfg$input_file))),
    paste("Samples:", nrow(d)), paste("Outer_Repeats:", cfg$n_repeats),
    paste("Frequency_Subsamples_Per_Holdout:", cfg$n_repeats),
    paste("Successful_Models:", n_success),
    paste("Failed_Models:", nrow(auc_df) - n_success),
    paste("Samples_With_OOF:", sum(oof_summary$Has_OOF)),
    paste("Samples_Without_OOF:", length(missing_ids)),
    paste("Missing_IDs:", paste(missing_ids, collapse = ";"))
  ))
  status_finalized <- TRUE
  if (cfg$verbose) {
    message(sprintf("已保存结果至 %s；成功模型 %d/%d，折外预测覆盖 %d/%d 个样本。",
                    cfg$output_dir, n_success, nrow(auc_df),
                    sum(oof_summary$Has_OOF), nrow(d)))
  }
  if (length(missing_ids) && cfg$require_complete_oof) {
    stop(sprintf(
      "%d 个样本尚无有效折外概率。结果已保存；请查看 oof_predictions_summary.csv 和 indexfilter_errors.csv。不会用训练集概率填补。",
      length(missing_ids)
    ), call. = FALSE)
  }
  invisible(bundle)
}

# source() / Rscript 均可直接启动。源文件不会调用其他项目脚本。
.frequency_is_cli <- sys.nframe() == 0L
.frequency_args <- if (.frequency_is_cli) commandArgs(trailingOnly = TRUE) else character()
if (.frequency_is_cli && any(.frequency_args %in% c("--help", "-h"))) {
  cat("Usage: Rscript FOREST.R [input.csv] [output_directory]\n")
  cat("默认输入 input.csv，默认保存到当前工作目录。参数可通过 frequency_config 或脚本顶部默认配置设置。\n")
} else if (!isTRUE(getOption("frequency_oof.functions_only", FALSE))) {
  .frequency_cfg <- if (exists("frequency_config", inherits = FALSE)) frequency_config else list()
  # 当前页面可直接定义 n_repeats；frequency_config 中的显式设置优先。
  if (is.list(.frequency_cfg) && !("n_repeats" %in% names(.frequency_cfg)) &&
      exists("n_repeats", inherits = FALSE)) {
    .frequency_cfg$n_repeats <- n_repeats
  }
  if (length(.frequency_args) > 2L) stop("最多接受两个位置参数：输入 CSV、输出目录。")
  if (length(.frequency_args) >= 1L) .frequency_cfg$input_file <- .frequency_args[1L]
  if (length(.frequency_args) >= 2L) .frequency_cfg$output_dir <- .frequency_args[2L]
  frequency_run <- run_frequency_pipeline(.frequency_cfg)
  results_df <- frequency_run$results_df
  model_df <- frequency_run$model_df
  auc_df <- frequency_run$auc_df
  curve_df <- frequency_run$curve_df
  oof_predictions <- frequency_run$oof_predictions
  oof_predictions_summary <- frequency_run$oof_predictions_summary
  indexfilterresult <- frequency_run$analysis$indexfilterresult
}

