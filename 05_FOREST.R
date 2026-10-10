#!/usr/bin/env Rscript
# FOREST.R
# FOREST: Feature-frequency Out-of-fold Resampling with Elastic-net Selection and Training.
# Standalone workflow: input.csv -> training-only differential/frequency selection -> repeated training and validation.
# Column 1: ID; column 2: Tumor/Control group; columns 3 onward: numeric features.
# Runs immediately when sourced; define frequency_config or n_repeats in the calling session to override defaults.
# Runs directly with the default settings; no preliminary parameter screening or prior result object is required.
# Set options(frequency_oof.functions_only = TRUE) to load functions without running the analysis.
# See section 5 of README.md for usage, output definitions, and failure handling.

frequency_default_config <- function() {
  list(
    input_file = "input.csv",
    output_dir = "outfile",         # All exports go into an outfile subdirectory
    positive_class = "Tumor",
    train_fraction = 0.67,
    TOP_n = 1000L,
    alpha_grid = seq(0.5, 1, by = 0.05), # 11 alpha values shared by frequency estimation and final training CV
    lambda_ratio_min = 0.001,
    n_repeats = 100L,               # Outer training/validation splits only
    n_lambda = 50L,
    pvalue_cutoff = 0.05,            # NULL disables p-value filtering
    diff_method = "wilcox",
    cv_nfolds = 5L,
    plot_iter = 1L,
    seed = 123L,
    coef_tol = 1e-8,
    fit_maxit = 1000000L,
    fpr_grid = seq(0, 1, length.out = 101L),
    roc_show_iterations = TRUE,     # Faint repeat curves, mean curves, and repeat-quantile ribbons
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
    alpha_grid,
    lambda_ratio_min_grid,
    TOP_n,
    n_repeats,
    n_lambda,
    pvalue_cutoff,
    diff_method = c("wilcox", "welch"),
    cv_nfolds = 5L,
    plot_iter = 1L,
    seed = 123L,
    coef_tol = 1e-8,
    fit_maxit = 1000000L,
    keep_differential = TRUE,
    verbose = TRUE,
    keep_all_repeats = TRUE,
    ensure_oof_coverage = FALSE) {

  for (pkg in c("glmnet", "pROC")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Please install R package: ", pkg, call. = FALSE)
    }
  }
  diff_method <- match.arg(diff_method)
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
    is.numeric(alpha_grid), length(alpha_grid) > 0L,
    all(is.finite(alpha_grid)),
    all(alpha_grid >= 0.5 & alpha_grid <= 1),
    !anyDuplicated(alpha_grid),
    is.numeric(lambda_ratio_min_grid), length(lambda_ratio_min_grid) > 0L,
    all(is.finite(lambda_ratio_min_grid)),
    all(lambda_ratio_min_grid > 0 & lambda_ratio_min_grid < 1),
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
  alphas <- sort(alpha_grid)
  minima <- unique(lambda_ratio_min_grid)
  tops <- unique(as.integer(TOP_n))
  parameter_grid <- expand.grid(
    train_fraction = fractions, TOP_n = tops,
    lambda_ratio_min = minima,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  parameter_grid <- cbind(Parameter = seq_len(nrow(parameter_grid)), parameter_grid)
  parameter_grid$alpha_grid <- paste(alphas, collapse = ";")
  n_holdouts <- length(fractions) * n_repeats
  stopifnot(seed + 4000 + n_holdouts <= .Machine$integer.max)
  n_evaluations <- nrow(parameter_grid) * n_repeats
  if (verbose) message(sprintf(
    paste0("Parameters: train_fraction={%s}, TOP_n={%s}, alpha={%s}, ",
           "lambda_ratio_min={%s}; %d parameter combinations, %d outer splits, %d planned final evaluations."),
    paste(fractions, collapse = ","), paste(tops, collapse = ","),
    paste(alphas, collapse = ","), paste(minima, collapse = ","),
    nrow(parameter_grid), n_holdouts, n_evaluations
  ))
  if (verbose) message(sprintf(
    "Each complete training set uses %d alpha values x %d lambda values; outer repeats n_repeats = %d.",
    length(alphas), n_lambda, n_repeats
  ))

  stratified_sample <- function(z, fraction) {
    sort(unlist(lapply(levels(z), function(g) {
      ii <- which(z == g)
      ii[sample.int(length(ii), floor(length(ii) * fraction))]
    }), use.names = FALSE))
  }
  # Generate outer splits first so random draws during training CV cannot change subsequent splits.
  split_plans <- lapply(seq_along(fractions), function(f) {
    set.seed(as.integer(seed + f - 1L))
    lapply(seq_len(n_repeats), function(i) stratified_sample(y, fractions[f]))
  })

  # Repair uncovered samples using only IDs, groups, and planned validation counts.
  # Preserve per-class training/validation counts within each split; do not use expression values or fitted results.
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
            "Group %s cannot cover all samples: %d validation samples per split require at least %d outer repeats; each validation class also needs at least 2 samples.",
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
      message(sprintf("Adjusted %d validation assignments to cover all samples; the training fraction and outer repeat count are preserved.",
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

  # Test the full feature set without preliminary MAD, missing-rate, variance, or topN filtering.
  # Test observed training values before imputation.
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
        reason <- "At least one group has fewer than 2 nonmissing observations."
      } else if (length(unique(c(a, b))) < 2L) {
        reason <- "All observed values are identical; the test cannot be performed."
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
          reason <- "The test did not return a valid p-value."
        }
      }
      data.frame(Feature = colnames(raw)[j], N_Negative = length(a),
                 N_Positive = length(b), P = p, Warning = note,
                 Reason = reason, stringsAsFactors = FALSE)
    })
    tab <- bind_list(rows)
    tab$Testable <- is.finite(tab$P)
    # Filter by raw p-values; exclude untestable features.
    tab$Significant <- tab$Testable & tab$P < pvalue_cutoff
    tab
  }

  # After differential testing, remove only columns that cannot be imputed or standardized; do not rank features here.
  preprocess_train <- function(raw) {
    med <- apply(raw, 2, stats::median, na.rm = TRUE)
    z <- raw[, is.finite(med), drop = FALSE]
    med <- med[colnames(z)]
    if (ncol(z) < 2L) stop("Fewer than 2 candidate features can be imputed.", call. = FALSE)
    for (j in seq_len(ncol(z))) z[is.na(z[, j]), j] <- med[j]
    sigma <- apply(z, 2, stats::sd)
    mu <- colMeans(z)
    keep <- is.finite(sigma) & sigma > 0 & is.finite(mu)
    z <- z[, keep, drop = FALSE]
    if (ncol(z) < 2L) stop("Fewer than 2 candidate features can be standardized.", call. = FALSE)
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
      stop("Invalid lambda_max; the path cannot be constructed.", call. = FALSE)
    }
    ratio <- exp(seq(0, log(ratio_min), length.out = n_lambda))
    list(lambda_max = lambda_max, lambda = lambda_max * ratio)
  }
  validate_path <- function(fit) {
    if (length(fit$jerr) != 1L || !is.finite(fit$jerr) || fit$jerr > 0L) {
      stop("glmnet returned a fatal or invalid error code.", call. = FALSE)
    }
    beta <- as.matrix(fit$beta)
    if (ncol(beta) != length(fit$lambda) || length(fit$a0) != ncol(beta)) {
      stop("The coefficient matrix does not match the returned lambda values.", call. = FALSE)
    }
    valid <- is.finite(fit$lambda) & fit$lambda > 0 &
      is.finite(fit$a0) & colSums(!is.finite(beta)) == 0L
    if (!any(valid)) stop("The path contains no valid models.", call. = FALSE)
    list(beta = beta, valid = valid)
  }

  diagnostic_rows <- final_diagnostic_rows <- list()
  # Count each lower ratio separately; fit one path per alpha on the complete outer training set.
  screen_frequency <- function(prep, yt, ratio_min, key, fraction, repeat_id) {
    candidates <- colnames(prep$x)
    counts <- matrix(0, nrow = length(candidates), ncol = length(alphas),
                     dimnames = list(candidates, as.character(alphas)))
    totals <- partial <- failed <- numeric(length(alphas))
    xs <- prep$x
    ys <- as.integer(yt == positive_class)
    for (a in seq_along(alphas)) {
      diag <- data.frame(
        Holdout = key, train_fraction = fraction, Repeat = repeat_id,
        alpha = alphas[a], lambda_ratio_min = ratio_min,
        N_Candidates = length(candidates), N_Usable = ncol(xs),
        Requested = n_lambda, Returned = 0L, Used = 0L,
        Jerr = NA_integer_, Status = "failed", Warning = "",
        Error = "", stringsAsFactors = FALSE
      )
      z <- capture_conditions({
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
      diag$Warning <- z$warning
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
        counts[ans$feats, a] <- counts[ans$feats, a] +
          rowSums(abs(beta) > coef_tol)
        totals[a] <- totals[a] + ncol(beta)
      }
      diagnostic_rows[[length(diagnostic_rows) + 1L]] <<- diag
    }
    list(counts = counts, totals = totals, partial = partial, failed = failed)
  }

  auc_roc <- function(response, probability) {
    stopifnot(length(response) == length(probability), all(is.finite(probability)))
    roc <- pROC::roc(response, as.numeric(probability), levels = levels(y),
                     direction = "<", quiet = TRUE)
    list(roc = roc, auc = as.numeric(pROC::auc(roc)))
  }
  # Fit the TOP_n input features on the outer training set; CV excludes the outer validation set.
  fit_final_alpha <- function(prep, feats, yt, alpha, ratio_min, foldid,
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
    # Compare only lambda values with valid held-out predictions for every training CV observation.
    stopifnot(is.matrix(fit$fit.preval), ncol(fit$fit.preval) == length(fit$lambda))
    valid <- is.finite(fit$cvm) & master_valid &
      colSums(!is.finite(fit$fit.preval)) == 0L
    if (!any(valid)) stop("The final model has no complete, valid training CV points.", call. = FALSE)
    best_deviance <- min(fit$cvm[valid])
    chosen <- max(fit$lambda[valid & fit$cvm == best_deviance])
    chosen_index <- match(chosen, fit$glmnet.fit$lambda)
    final_nz <- sum(abs(checked$beta[, chosen_index]) > coef_tol)
    fit$fit.preval <- NULL
    list(fit = fit, alpha = alpha, final_lambda = chosen,
         final_ratio = chosen / grid$lambda_max, cv_deviance = best_deviance,
         n_nonzero = final_nz, jerr = fit$glmnet.fit$jerr)
  }

  # Use shared training CV folds for all alpha values; select alpha/lambda before predicting validation samples.
  fit_final <- function(prep, feats, yt, raw_val, yv, ratio_min, foldid,
                        expected_top_n, parameter_id, key, repeat_id) {
    attempts <- lapply(alphas, function(alpha) capture_conditions(fit_final_alpha(
      prep, feats, yt, alpha, ratio_min, foldid, expected_top_n
    )))
    ok <- vapply(attempts, function(z) !inherits(z$value, "error"), logical(1))
    chosen_a <- NA_integer_
    if (any(ok)) {
      valid_a <- which(ok)
      dev <- vapply(attempts[valid_a], function(z) z$value$cv_deviance, numeric(1))
      # Break exact ties by larger alpha; within each alpha, prefer larger lambda.
      chosen_a <- valid_a[order(dev, -alphas[valid_a])[1L]]
    }
    for (a in seq_along(alphas)) {
      z <- attempts[[a]]
      final_diagnostic_rows[[length(final_diagnostic_rows) + 1L]] <<- data.frame(
        Parameter = parameter_id, Holdout = key, Repeat = repeat_id,
        alpha = alphas[a], CV_Deviance = if (ok[a]) z$value$cv_deviance else NA_real_,
        Final_Lambda = if (ok[a]) z$value$final_lambda else NA_real_,
        Selected = !is.na(chosen_a) && a == chosen_a,
        Warning = z$warning, Jerr = if (ok[a]) z$value$jerr else NA_integer_,
        Error = if (ok[a]) "" else conditionMessage(z$value),
        stringsAsFactors = FALSE
      )
    }
    if (!any(ok)) stop("Training CV failed for all alpha values; see final_diagnostics.csv.", call. = FALSE)
    ans <- attempts[[chosen_a]]$value
    notes <- vapply(seq_along(alphas), function(a) {
      z <- attempts[[a]]
      if (!ok[a]) sprintf("Training CV failed for alpha=%.2f: %s", alphas[a], conditionMessage(z$value)) else z$warning
    }, character(1))
    ans$warning <- paste(unique(notes[nzchar(notes)]), collapse = " | ")
    fit <- ans$fit
    chosen <- ans$final_lambda
    xt <- prep$x[, feats, drop = FALSE]
    xv <- transform_new(raw_val, prep, feats)
    stopifnot(ncol(xv) == expected_top_n, identical(colnames(xv), colnames(xt)))
    pt <- as.numeric(stats::predict(fit, xt, s = chosen, type = "response"))
    pv <- as.numeric(stats::predict(fit, xv, s = chosen, type = "response"))
    trn <- auc_roc(yt, pt)
    val <- auc_roc(yv, pv)
    ans$pt <- pt
    ans$pv <- pv
    ans$trn <- trn
    ans$val <- val
    ans
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
      if (verbose) message(sprintf("Outer split [%d/%d] train_fraction=%.2f, repeat=%d",
                                   h_id, n_holdouts, fraction, i))

      outer <- capture_conditions({
        if (min(table(yt)) < 3L || min(table(y[va])) < 2L) {
          stop("This fraction cannot provide at least 3 training and 2 validation samples per class.", call. = FALSE)
        }
        differential <- differential_train(raw_train, yt)
        sig <- differential$Feature[differential$Significant]
        if (keep_differential) differential_tables[[length(differential_tables) + 1L]] <-
          cbind(Holdout = key, train_fraction = fraction, Repeat = i, differential)
        if (length(sig) < 2L) stop(
          if (is.null(pvalue_cutoff)) "Fewer than 2 training candidate features are available." else
            sprintf("Fewer than 2 training features have p-value < %.3g.", pvalue_cutoff),
          call. = FALSE
        )
        prep <- preprocess_train(raw_train[, sig, drop = FALSE])
        candidates <- colnames(prep$x)
        # Use the same training CV folds for all TOP_n, alpha, and lambda lower-ratio settings.
        set.seed(as.integer(seed + 2000 + h_id))
        folds <- as.integer(min(cv_nfolds, min(table(yt))))
        foldid <- integer(length(yt))
        for (g in levels(yt)) {
          ii <- which(yt == g)
          foldid[ii[sample.int(length(ii))]] <- rep(seq_len(folds), length.out = length(ii))
        }
        list(differential = differential, prep = prep, candidates = candidates,
             foldid = foldid, n_diff = length(sig))
      })
      outer_error <- if (inherits(outer$value, "error")) conditionMessage(outer$value) else ""

      selection <- list()
      if (!nzchar(outer_error)) {
        o <- outer$value
        for (m in seq_along(minima)) {
          screened <- screen_frequency(o$prep, yt, minima[m], key, fraction, i)
          scopes <- list(seq_along(alphas))
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
              n_requested = length(aa) * n_lambda,
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
          alpha_fixed = NA_real_,
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
          sel <- selection[[m]][[1L]]
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
            row$Error <- "No valid models are available for this frequency scope."
          } else if (length(eligible) < pars$TOP_n) {
            row$Error <- sprintf("Only %d candidate features have positive frequency, fewer than TOP_n=%d.",
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
              pars$lambda_ratio_min, o$foldid,
              expected_top_n = pars$TOP_n, parameter_id = pars$Parameter,
              key = key, repeat_id = i
            ))
            row$Warning <- paste(c(row$Warning, z$warning)[
              nzchar(c(row$Warning, z$warning))], collapse = " | ")
            if (inherits(z$value, "error")) {
              row$Error <- conditionMessage(z$value)
            } else {
              ans <- z$value
              row$Warning <- paste(c(row$Warning, ans$warning)[
                nzchar(c(row$Warning, ans$warning))], collapse = " | ")
              row$alpha_fixed <- ans$alpha
              row$Final_Lambda <- ans$final_lambda
              row$Final_Lambda_Ratio <- ans$final_ratio
              row$CV_Deviance <- ans$cv_deviance
              row$N_Final_Nonzero <- ans$n_nonzero
              row[["TRN-AUC"]] <- ans$trn$auc
              row[["VAL-AUC"]] <- ans$val$auc
              # Save all repeats when requested; otherwise retain only plot_iter to limit memory use.
              if (keep_all_repeats || i == plot_iter) {
                model_key <- sprintf("P%03d_%s", pars$Parameter, key)
                models[[model_key]] <- list(
                  Parameter = pars$Parameter, Holdout = key, Repeat = i,
                  fit = ans$fit, feats = feats, alpha_fixed = ans$alpha,
                  final_lambda = ans$final_lambda,
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
            paste0("  Final training/validation [%d/%d] train_fraction=%.2f, repeat=%d, ",
                   "TOP_n=%d, alpha=%.3g, lambda_ratio_min=%.3g"),
            evaluation_id, n_evaluations, fraction, i,
            pars$TOP_n, row$alpha_fixed, pars$lambda_ratio_min
          )
          if (nzchar(row$Error)) {
            message(prefix, "; failed: ", row$Error)
          } else {
            message(prefix, sprintf(
              "; input features=%d, final nonzero features=%d, TRN-AUC=%.4f, VAL-AUC=%.4f",
              row$N_Used, row$N_Final_Nonzero, row[["TRN-AUC"]], row[["VAL-AUC"]]
            ))
          }
        }
        result_rows[[evaluation_id]] <- row
      }
      if (verbose && nzchar(outer_error)) message("  All parameter combinations failed for this split: ", outer_error)
    }
  }
  holdout_results <- bind_list(result_rows)
  # Successful models must use their specified TOP_n input features; fewer coefficients may remain nonzero.
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
    alpha = numeric(), lambda_ratio_min = numeric(),
    N_Candidates = integer(), N_Usable = integer(), Requested = integer(),
    Returned = integer(), Used = integer(), Jerr = integer(), Status = character(),
    Warning = character(), Error = character(), stringsAsFactors = FALSE
  )
  final_diagnostics <- bind_list(final_diagnostic_rows)
  if (!nrow(final_diagnostics)) final_diagnostics <- data.frame(
    Parameter = integer(), Holdout = character(), Repeat = integer(),
    alpha = numeric(), CV_Deviance = numeric(), Final_Lambda = numeric(),
    Selected = logical(),
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
                  alpha_grid = alphas, lambda_ratio_min_grid = minima,
                  n_repeats = n_repeats,
                  n_lambda = n_lambda,
                  pvalue_cutoff = pvalue_cutoff,
                  diff_method = diff_method,
                  selection_alpha_grid = alphas, cv_nfolds = cv_nfolds,
                  plot_iter = plot_iter, seed = seed, coef_tol = coef_tol,
                  fit_maxit = fit_maxit, keep_differential = keep_differential,
                  keep_all_repeats = keep_all_repeats,
                  ensure_oof_coverage = ensure_oof_coverage)
  )
}

# Predict new samples using a saved model; probabilities always refer to positive_class.
predict_frequency_model <- function(model, newdata) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Please install the glmnet R package.")
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

# Explicit single-repeat display; use plot_frequency_roc() for repeat stability.
plot_frequency_single_roc <- function(result, parameter_id = 1L) {
  eligible <- vapply(result$models, function(m) {
    m$Parameter == parameter_id && m$Repeat == result$config$plot_iter
  }, logical(1))
  keys <- names(result$models)[eligible]
  if (!length(keys)) stop("No successful saved model exists for this parameter in the plot_iter split.", call. = FALSE)
  r <- result$roc_examples[[keys[1L]]]
  graphics::plot(r$TRN, col = "#2878B5", lwd = 2, legacy.axes = TRUE,
                 main = sprintf("Parameter %d, repeat %d", parameter_id, result$config$plot_iter))
  graphics::plot(r$VAL, col = "#D95319", lwd = 2, add = TRUE)
  graphics::legend("bottomright", bty = "n", lwd = 2, col = c("#2878B5", "#D95319"),
                   legend = c(sprintf("TRN AUC = %.3f", as.numeric(pROC::auc(r$TRN))),
                              sprintf("VAL AUC = %.3f", as.numeric(pROC::auc(r$VAL)))))
  invisible(r)
}

# The default plotting entry point always reads every successful repeat.
plot_frequency_roc <- function(
    result, parameter_id = 1L,
    fpr_grid = seq(0, 1, length.out = 101L), show_iterations = TRUE) {
  records <- frequency_roc_records(result, fpr_grid, parameter_id)
  roc_results <- calculate_resampled_roc(
    records$results_df, fpr_grid, records$positive_class,
    reference_auc = records$auc_df, reference_curves = records$curve_df,
    expected_n_repeats = records$n_repeats
  )
  .frequency_report_roc_records(roc_results)
  roc_plot <- plot_frequency_resampled_roc(roc_results, show_iterations)
  print(roc_plot)
  invisible(list(roc_results = roc_results, roc_plot = roc_plot, records = records))
}


frequency_result_lists <- function(result, fpr_grid = seq(0, 1, length.out = 101L)) {
  if (!isTRUE(result$config$keep_all_repeats)) {
    stop("Rerun the analysis with keep_all_repeats=TRUE before exporting details for all repeats.",
         call. = FALSE)
  }
  stopifnot(is.numeric(fpr_grid), length(fpr_grid) >= 2L,
            all(is.finite(fpr_grid)), all(fpr_grid >= 0 & fpr_grid <= 1),
            all(diff(fpr_grid) > 0))
  add_parameters <- function(z) {
    cols <- intersect(c("train_fraction", "TOP_n", "alpha_grid", "alpha_fixed", "lambda_ratio_min"),
                        names(result$holdout_results))
    key <- paste(result$holdout_results$Parameter, result$holdout_results$Holdout, sep = "::")
    stopifnot(!anyDuplicated(key))
    ii <- match(paste(z$Parameter, z$Holdout, sep = "::"), key)
    stopifnot(!anyNA(ii))
    out <- cbind(z, result$holdout_results[ii, cols, drop = FALSE])
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
  # Include successful and failed splits; failed AUCs are NA and reasons appear in Error.
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
    # Include the intercept and zero coefficients; feature coefficients use the training-standardized scale.
    model_list[[length(model_list) + 1L]] <- add_parameters(data.frame(
      Parameter = m$Parameter, Holdout = m$Holdout, Repeat = m$Repeat,
      Feature = rownames(b), Coefficient = as.numeric(b[, 1L]),
      Final_Lambda = m$final_lambda, stringsAsFactors = FALSE
    ))
    rr <- result$roc_examples[[key]]
    stopifnot(!is.null(rr$TRN), !is.null(rr$VAL))
    for (set_name in c("TRN", "VAL")) {
      r <- rr[[set_name]]
      # Use the common FPR grid for curve display only; pROC calculates AUC from the original predictions.
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


# Use the same saved-list retrieval sequence as the formal analysis code.
frequency_roc_records <- function(
    result, fpr_grid = seq(0, 1, length.out = 101L), parameter_id = NULL) {
  if (!is.null(result$analysis)) result <- result$analysis
  output_lists <- frequency_result_lists(result, fpr_grid = fpr_grid)
  records <- list(
    results_df = .frequency_bind_rows(output_lists$results_list),
    model_df = .frequency_bind_rows(output_lists$model_list),
    auc_df = .frequency_bind_rows(output_lists$auc_list),
    curve_df = .frequency_bind_rows(output_lists$curve_list)
  )
  if (!is.null(parameter_id)) {
    if (!all(parameter_id %in% unique(records$auc_df$Parameter))) {
      stop("Unknown Parameter requested for ROC record retrieval.", call. = FALSE)
    }
    for (nm in names(records)) {
      records[[nm]] <- records[[nm]][records[[nm]]$Parameter %in% parameter_id, , drop = FALSE]
    }
  }
  if (!nrow(records$results_df)) stop("No successful saved repeat predictions are available.", call. = FALSE)
  records$positive_class <- result$config$positive_class
  records$n_repeats <- result$config$n_repeats
  records
}


# Follow the per-repeat ROC -> common FPR grid -> mean/quantile workflow.
# The ribbon is a pointwise 2.5%-97.5% empirical repeat interval, not a CI for the mean.
calculate_resampled_roc <- function(
    pred_data,
    fpr_grid = seq(0, 1, length.out = 101L),
    positive_class = "Tumor",
    reference_auc = NULL,
    reference_curves = NULL,
    expected_n_repeats = NULL) {
  required <- c("Set", "Repeat", "Parameter", "Group", "Probability")
  stopifnot(
    is.data.frame(pred_data), all(required %in% names(pred_data)),
    nrow(pred_data) > 0L, !anyNA(pred_data[, required, drop = FALSE]),
    length(positive_class) == 1L, positive_class %in% c("Control", "Tumor"),
    all(pred_data$Set %in% c("TRN", "VAL")),
    all(pred_data$Group %in% c("Control", "Tumor")),
    is.numeric(pred_data$Probability), all(is.finite(pred_data$Probability)),
    all(pred_data$Probability >= 0 & pred_data$Probability <= 1),
    is.numeric(fpr_grid), length(fpr_grid) >= 2L,
    all(is.finite(fpr_grid)), all(fpr_grid >= 0 & fpr_grid <= 1),
    all(diff(fpr_grid) > 0)
  )
  if (!is.null(expected_n_repeats)) {
    stopifnot(is.numeric(expected_n_repeats), length(expected_n_repeats) == 1L,
              is.finite(expected_n_repeats), expected_n_repeats >= 2L,
              expected_n_repeats == floor(expected_n_repeats))
  }
  for (pkg in c("pROC", "dplyr", "tidyr", "purrr", "tibble")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Please install R package: ", pkg)
  }
  `%>%` <- dplyr::`%>%`
  class_levels <- c(setdiff(c("Control", "Tumor"), positive_class), positive_class)
  grouped_roc <- pred_data %>%
    dplyr::group_by(Parameter, Set, Repeat) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      roc_obj = purrr::map(data, function(z) {
        if (!setequal(unique(as.character(z$Group)), class_levels)) {
          stop("Each Parameter / Set / Repeat must contain both classes.", call. = FALSE)
        }
        if ("ID" %in% names(z) && anyDuplicated(z$ID)) {
          stop("A sample occurs more than once within a Parameter / Set / Repeat.", call. = FALSE)
        }
        pROC::roc(response = z$Group, predictor = z$Probability,
                   levels = class_levels, direction = "<", quiet = TRUE)
      }),
      auc = purrr::map_dbl(roc_obj, ~ as.numeric(pROC::auc(.x))),
      tpr_grid = purrr::map(roc_obj, function(r) {
        stats::approx(x = 1 - r$specificities, y = r$sensitivities,
                       xout = fpr_grid, ties = max, rule = 2)$y
      })
    ) %>%
    dplyr::ungroup()

  # AUC always comes from each repeat's original predictions.
  auc_per_iter <- grouped_roc %>% dplyr::select(Parameter, Set, Repeat, auc)
  repeat_counts <- auc_per_iter %>%
    dplyr::count(Parameter, Set, name = "n_repeats")
  if (any(repeat_counts$n_repeats < 2L)) {
    bad <- repeat_counts[repeat_counts$n_repeats < 2L, , drop = FALSE]
    stop("A stability interval needs at least two recorded repeats per Parameter/Set. Found: ",
           paste(paste0("Parameter=", bad$Parameter, ", Set=", bad$Set, ", n=", bad$n_repeats),
                 collapse = " | "),
           ". Retrieve every successful repeat with keep_all_repeats=TRUE; a plot_iter-only record is insufficient.",
           call. = FALSE)
  }
  expected_success <- NULL
  if (!is.null(reference_auc)) {
    required_auc <- c("Parameter", "Repeat", "TRN-AUC", "VAL-AUC", "Error")
    stopifnot(is.data.frame(reference_auc), all(required_auc %in% names(reference_auc)),
              !anyDuplicated(reference_auc[, c("Parameter", "Repeat"), drop = FALSE]))
    expected_success <- reference_auc[
      !is.na(reference_auc$Error) & reference_auc$Error == "" &
        is.finite(reference_auc[["TRN-AUC"]]) & is.finite(reference_auc[["VAL-AUC"]]), , drop = FALSE]
    wanted <- dplyr::bind_rows(
      tibble::tibble(Parameter = expected_success$Parameter, Set = "TRN",
                     Repeat = expected_success$Repeat, Original_AUC = expected_success[["TRN-AUC"]]),
      tibble::tibble(Parameter = expected_success$Parameter, Set = "VAL",
                     Repeat = expected_success$Repeat, Original_AUC = expected_success[["VAL-AUC"]])
    )
    missing <- dplyr::anti_join(wanted, auc_per_iter, by = c("Parameter", "Set", "Repeat"))
    extra <- dplyr::anti_join(auc_per_iter, wanted, by = c("Parameter", "Set", "Repeat"))
    if (nrow(missing) || nrow(extra)) {
      stop("ROC prediction records do not match all successful model attempts: ",
             nrow(missing), " Parameter/Set/Repeat records are missing and ",
             nrow(extra), " are unexpected. Retrieve the complete results_list before plotting.", call. = FALSE)
    }
    checked_auc <- dplyr::left_join(auc_per_iter, wanted, by = c("Parameter", "Set", "Repeat"))
    if (any(!is.finite(checked_auc$Original_AUC)) ||
        max(abs(checked_auc$auc - checked_auc$Original_AUC)) > 1e-10) {
      stop("Recomputed per-repeat AUC does not match the saved auc_list records.", call. = FALSE)
    }
  }
  avg_roc <- grouped_roc %>%
    dplyr::group_by(Parameter, Set) %>%
    dplyr::summarise(
      n_repeats = dplyr::n(), auc_mean = mean(auc), auc_sd = stats::sd(auc),
      roc_avg = list(tibble::tibble(
        FPR = fpr_grid, TPR = rowMeans(do.call(cbind, tpr_grid))
      )),
      .groups = "drop"
    )
  ci_data <- grouped_roc %>%
    dplyr::group_by(Parameter, Set) %>%
    dplyr::summarise(
      ci_bounds = list({
        tpr_mat <- do.call(cbind, tpr_grid)
        tibble::tibble(
          FPR = fpr_grid, TPR_mean = rowMeans(tpr_mat),
          TPR_lower = apply(tpr_mat, 1, stats::quantile, probs = 0.025, names = FALSE),
          TPR_upper = apply(tpr_mat, 1, stats::quantile, probs = 0.975, names = FALSE)
        )
      }),
      .groups = "drop"
    )
  auc_summary <- auc_per_iter %>%
    dplyr::group_by(Parameter, Set) %>%
    dplyr::summarise(
      n_repeats = dplyr::n(),
      auc_mean = mean(auc, na.rm = TRUE), auc_sd = stats::sd(auc, na.rm = TRUE),
      auc_lower = unname(stats::quantile(auc, 0.025, na.rm = TRUE)),
      auc_upper = unname(stats::quantile(auc, 0.975, na.rm = TRUE)),
      .groups = "drop"
    )
  iteration_curves <- grouped_roc %>%
    dplyr::select(Parameter, Set, Repeat, tpr_grid) %>%
    dplyr::mutate(roc_curve = purrr::map(tpr_grid, ~ tibble::tibble(FPR = fpr_grid, TPR = .x))) %>%
    dplyr::select(Parameter, Set, Repeat, roc_curve) %>%
    tidyr::unnest(roc_curve)
  if (!is.null(reference_curves)) {
    keys <- c("Parameter", "Set", "Repeat", "FPR")
    stopifnot(is.data.frame(reference_curves), all(c(keys, "TPR") %in% names(reference_curves)),
              !anyDuplicated(reference_curves[, keys, drop = FALSE]))
    checked <- dplyr::full_join(
      iteration_curves, reference_curves[, c(keys, "TPR"), drop = FALSE],
      by = keys, suffix = c(".Computed", ".Recorded")
    )
    if (any(!is.finite(checked$TPR.Computed)) || any(!is.finite(checked$TPR.Recorded)) ||
        max(abs(checked$TPR.Computed - checked$TPR.Recorded)) > 1e-10) {
      stop("Per-repeat ROC curves do not match the saved curve_list on the requested FPR grid.", call. = FALSE)
    }
  }
  plot_data <- avg_roc %>% tidyr::unnest(roc_avg) %>%
    dplyr::left_join(
      ci_data %>% tidyr::unnest(ci_bounds) %>%
        dplyr::select(Parameter, Set, FPR, TPR_lower, TPR_upper),
      by = c("Parameter", "Set", "FPR")
    ) %>%
    dplyr::select(Parameter, Set, FPR, TPR, TPR_lower, TPR_upper)
  record_checks <- lapply(seq_len(nrow(auc_summary)), function(i) {
    g <- auc_summary[i, , drop = FALSE]
    aa <- auc_per_iter[auc_per_iter$Parameter == g$Parameter & auc_per_iter$Set == g$Set, , drop = FALSE]
    cc <- grouped_roc[grouped_roc$Parameter == g$Parameter & grouped_roc$Set == g$Set, , drop = FALSE]
    band <- plot_data[plot_data$Parameter == g$Parameter & plot_data$Set == g$Set, , drop = FALSE]
    width <- band$TPR_upper - band$TPR_lower
    unique_curves <- length(unique(purrr::map_chr(cc$tpr_grid, ~ paste(format(.x, digits = 17), collapse = ";"))))
    planned <- if (is.null(expected_n_repeats)) NA_integer_ else as.integer(expected_n_repeats)
    expected <- failed <- NA_integer_
    if (!is.null(reference_auc)) {
      attempts <- reference_auc[reference_auc$Parameter == g$Parameter, , drop = FALSE]
      if (is.na(planned)) planned <- nrow(attempts)
      expected <- sum(expected_success$Parameter == g$Parameter)
      failed <- nrow(attempts) - expected
    }
    note <- if (unique_curves == 1L) {
      "All recorded ROC curves have identical grid coordinates; the interval has zero width."
    } else if (max(width) == 0) {
      "The central 95% repeat interval has zero width despite variation in some recorded curves."
    } else ""
    data.frame(
      Parameter = g$Parameter, Set = as.character(g$Set),
      N_Planned_Repeats = planned, N_Recorded_Repeats = nrow(aa),
      N_Expected_Successful_Repeats = expected, N_Failed_Repeats = failed,
      Unrecorded_Repeat_IDs = if (is.na(planned)) "" else paste(setdiff(seq_len(planned), aa$Repeat), collapse = ";"),
      AUC_Mean = g$auc_mean, AUC_SD = g$auc_sd,
      AUC_Lower = g$auc_lower, AUC_Upper = g$auc_upper,
      N_Unique_ROC_Curves = unique_curves,
      Mean_TPR_Interval_Width = mean(width), Max_TPR_Interval_Width = max(width),
      Interval_Note = note, stringsAsFactors = FALSE
    )
  })
  record_check <- .frequency_bind_rows(record_checks)
  list(avg_roc = avg_roc, ci_data = ci_data, auc_per_iter = auc_per_iter,
       auc_summary = auc_summary, iteration_curves = iteration_curves,
       plot_data = plot_data, record_check = record_check,
       interval_type = "Pointwise 2.5%-97.5% empirical quantiles across completed repeats")
}


plot_frequency_resampled_roc <- function(roc_results, show_iterations = TRUE) {
  for (pkg in c("ggplot2", "dplyr", "tidyr", "tibble")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Please install R package: ", pkg)
  }
  `%>%` <- dplyr::`%>%`
  stopifnot(is.logical(show_iterations), length(show_iterations) == 1L,
            !is.na(show_iterations))
  add_origin <- function(z, key_columns) {
    origin <- z[!duplicated(z[, key_columns, drop = FALSE]), , drop = FALSE]
    for (nm in intersect(c("FPR", "TPR", "TPR_lower", "TPR_upper"), names(origin))) {
      origin[[nm]] <- 0
    }
    origin$.origin <- 0L
    z$.origin <- 1L
    z <- rbind(origin, z)
    # Keep the real FPR=0 coordinate after the artificial origin to preserve the vertical segment.
    ordering <- c(lapply(key_columns, function(nm) z[[nm]]), list(z$FPR, z$.origin))
    z <- z[do.call(order, ordering), , drop = FALSE]
    z$Set <- factor(z$Set, levels = c("TRN", "VAL"))
    z
  }
  # Reproduce the supplied unnest -> join -> origin -> mean-curve plotting sequence.
  # Retain Parameter in the join/group keys to avoid mixing parameter combinations.
  plot_data <- roc_results$avg_roc %>% tidyr::unnest(roc_avg) %>%
    dplyr::left_join(
      roc_results$ci_data %>% tidyr::unnest(ci_bounds) %>%
        dplyr::select(Parameter, Set, FPR, TPR_lower, TPR_upper),
      by = c("Parameter", "Set", "FPR")
    ) %>%
    dplyr::group_by(Parameter, Set) %>%
    dplyr::arrange(FPR, .by_group = TRUE) %>%
    dplyr::group_modify(~ tibble::tibble(
      FPR = c(0, .x$FPR), TPR = c(0, .x$TPR),
      TPR_lower = c(0, .x$TPR_lower), TPR_upper = c(0, .x$TPR_upper)
    )) %>%
    dplyr::ungroup()
  plot_data$Set <- factor(plot_data$Set, levels = c("TRN", "VAL"))
  annotations <- roc_results$auc_per_iter %>%
    dplyr::group_by(Parameter, Set) %>%
    dplyr::summarise(
      auc_mean = mean(auc, na.rm = TRUE),
      auc_lower = stats::quantile(auc, 0.025, na.rm = TRUE),
      auc_upper = stats::quantile(auc, 0.975, na.rm = TRUE),
      n_repeats = dplyr::n(), .groups = "drop"
    )
  annotations$Set <- factor(annotations$Set, levels = c("TRN", "VAL"))
  annotations$Label <- sprintf("%s: %.3f (%.3f\u2013%.3f); n=%d", annotations$Set,
                                annotations$auc_mean, annotations$auc_lower,
                                annotations$auc_upper, annotations$n_repeats)
  annotations$Label_Y <- ifelse(annotations$Set == "TRN", 0.125, 0.05)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = FPR, y = TPR,
                                             color = Set, fill = Set, group = Set)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = TPR_lower, ymax = TPR_upper),
                         alpha = 0.25, color = NA, show.legend = FALSE)
  if (show_iterations) {
    curves <- add_origin(roc_results$iteration_curves, c("Parameter", "Set", "Repeat"))
    p <- p + ggplot2::geom_path(data = curves,
                                ggplot2::aes(group = interaction(Parameter, Set, Repeat)),
                                alpha = 0.045, linewidth = 0.18, show.legend = FALSE)
  }
  p <- p + ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                          color = "grey50", linewidth = 0.25) +
    ggplot2::scale_color_manual(values = c(TRN = "red", VAL = "blue"), drop = FALSE) +
    ggplot2::scale_fill_manual(values = c(TRN = "red", VAL = "blue"), drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(-0.025, 1.025), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-0.025, 1.025), expand = c(0, 0)) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::labs(x = "1 - Specificity", y = "Sensitivity",
                    caption = "Shading: pointwise central 95% interval across recorded repeats.") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "top", legend.title = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey80", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 10, face = "plain"),
      plot.caption = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.35)
    ) +
    ggplot2::geom_text(data = annotations,
                        ggplot2::aes(x = 0.175, y = Label_Y, label = Label, color = Set),
                        inherit.aes = FALSE, hjust = 0, size = 3, fontface = "bold",
                        show.legend = FALSE)
  if (length(unique(plot_data$Parameter)) > 1L) {
    p <- p + ggplot2::facet_wrap(~ Parameter, labeller = "label_both")
  }
  p
}


# "Best" means the highest observed outer-validation AUC; training CV deviance breaks ties.
# This model is a descriptive example, not a replacement for the all-repeat ROC summary.
select_frequency_best_model <- function(result) {
  h <- result$holdout_results
  stopifnot(all(c("Parameter", "Holdout", "Repeat", "VAL-AUC", "CV_Deviance", "Error") %in% names(h)))
  keys <- sprintf("P%03d_%s", h$Parameter, h$Holdout)
  ok <- !is.na(h$Error) & h$Error == "" & is.finite(h[["VAL-AUC"]]) &
    is.finite(h$CV_Deviance) & keys %in% names(result$models)
  if (!any(ok)) stop("No successful saved model is available for the nomogram.", call. = FALSE)
  h <- h[ok, , drop = FALSE]
  h <- h[order(-h[["VAL-AUC"]], h$CV_Deviance, h$Parameter, h$Repeat), , drop = FALSE]
  summary <- h[1L, , drop = FALSE]
  key <- sprintf("P%03d_%s", summary$Parameter, summary$Holdout)
  summary$Model_Key <- key
  summary$Selection_Rule <- "Highest validation AUC; ties: lowest training CV deviance, Parameter, Repeat"
  model <- result$models[[key]]
  model$selection_summary <- summary
  list(model = model, summary = summary)
}


# Construct an exact logistic nomogram from the saved glmnet coefficients.
# Raw-scale coefficients remove training standardization without refitting the model.
frequency_nomogram_data <- function(model, data) {
  stopifnot(is.data.frame(data), "ID" %in% names(data), !anyDuplicated(data$ID),
            all(model$feats %in% names(data)),
            all(model$train_ID %in% data$ID))
  b <- as.matrix(stats::coef(model$fit, s = model$final_lambda))
  stopifnot(nrow(b) == length(model$feats) + 1L,
            all(model$feats %in% rownames(b)[-1L]))
  beta <- b[-1L, 1L][match(model$feats, rownames(b)[-1L])]
  sigma <- model$preprocess$sigma[model$feats]
  mu <- model$preprocess$mu[model$feats]
  med <- model$preprocess$median[model$feats]
  stopifnot(all(is.finite(beta)), all(is.finite(sigma)), all(sigma > 0),
            all(is.finite(mu)), all(is.finite(med)))
  training <- data[match(model$train_ID, data$ID), model$feats, drop = FALSE]
  stopifnot(all(vapply(training, is.numeric, logical(1))))
  raw <- as.matrix(training)
  for (j in seq_along(model$feats)) raw[is.na(raw[, j]), j] <- med[j]
  stopifnot(all(is.finite(raw)))
  raw_beta <- beta / sigma
  raw_intercept <- as.numeric(b[1L, 1L] - sum(beta * mu / sigma))
  lo <- apply(raw, 2L, min)
  hi <- apply(raw, 2L, max)
  lower_term <- pmin(raw_beta * lo, raw_beta * hi)
  span <- abs(raw_beta) * (hi - lo)
  # Retain every exactly nonzero coefficient so the point system matches the fitted model.
  active <- which(beta != 0 & span > 0)
  points_per_logit <- if (length(active)) 100 / max(span[active]) else NA_real_
  eta_min <- raw_intercept + sum(lower_term)
  total_max <- if (length(active)) sum(span[active]) * points_per_logit else 0
  if (length(active)) stopifnot(is.finite(points_per_logit), is.finite(total_max), total_max > 0)
  coefficient_table <- data.frame(
    Feature = model$feats, Standardized_Coefficient = as.numeric(beta),
    OriginalScale_Coefficient = as.numeric(raw_beta),
    Training_Median = as.numeric(med), Training_Mean = as.numeric(mu),
    Training_SD = as.numeric(sigma), Training_Min = as.numeric(lo), Training_Max = as.numeric(hi),
    Max_Points = if (length(active)) as.numeric(span * points_per_logit) else 0,
    stringsAsFactors = FALSE
  )
  intercept <- coefficient_table[1L, , drop = FALSE]
  intercept[1L, ] <- NA
  intercept$Feature <- "(Intercept)"
  intercept$Standardized_Coefficient <- as.numeric(b[1L, 1L])
  intercept$OriginalScale_Coefficient <- raw_intercept
  coefficient_table <- rbind(intercept, coefficient_table)
  # Check both the original-scale formula and the total-points formula against glmnet.
  expected <- predict_frequency_model(model, training)
  eta <- as.numeric(raw_intercept + raw %*% raw_beta)
  stopifnot(max(abs(stats::plogis(eta) - expected)) < 1e-8)
  if (length(active)) {
    total <- as.numeric((raw %*% raw_beta - sum(lower_term)) * points_per_logit)
    stopifnot(max(abs(stats::plogis(eta_min + total / points_per_logit) - expected)) < 1e-8)
  }
  list(coefficients = coefficient_table, active = active,
       feats = model$feats, raw_beta = raw_beta, raw_intercept = raw_intercept,
       training_min = lo, training_max = hi, lower_term = lower_term,
       points_per_logit = points_per_logit, eta_min = eta_min, total_max = total_max,
       constant_probability = stats::plogis(raw_intercept),
       positive_class = model$positive_class)
}


plot_frequency_nomogram <- function(model, nomogram_data) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install the ggplot2 R package.")
  nd <- nomogram_data
  s <- model$selection_summary
  subtitle <- sprintf("Parameter %d | Repeat %d | alpha %.2f | lambda %.3g | validation AUC %.3f",
                       model$Parameter, model$Repeat, model$alpha_fixed,
                       model$final_lambda, s[["VAL-AUC"]])
  if (!length(nd$active)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                         label = sprintf("Intercept-only model\nP(%s) = %.4f",
                                           model$positive_class, nd$constant_probability), size = 5) +
      ggplot2::labs(title = "FOREST Best-Model Nomogram", subtitle = subtitle))
  }
  n_rows <- length(nd$active) + 3L
  axis_rows <- tick_rows <- list()
  add_axis <- function(label, y, x_end, x, labels) {
    axis_rows[[length(axis_rows) + 1L]] <<- data.frame(
      Label = label, Y = y, X = 0, X_end = x_end, stringsAsFactors = FALSE
    )
    tick_rows[[length(tick_rows) + 1L]] <<- data.frame(
      X = x, Y = y, Label = labels, stringsAsFactors = FALSE
    )
  }
  tick_labels <- function(x) format(signif(x, 4L), trim = TRUE, scientific = FALSE)
  add_axis("Points", n_rows, 100, seq(0, 100, by = 20), as.character(seq(0, 100, by = 20)))
  for (i in seq_along(nd$active)) {
    j <- nd$active[i]
    limits <- c(nd$training_min[j], nd$training_max[j])
    ticks <- pretty(limits, n = 4L)
    ticks <- ticks[ticks >= limits[1L] & ticks <= limits[2L]]
    if (length(ticks) < 2L) ticks <- limits
    points <- (nd$raw_beta[j] * ticks - nd$lower_term[j]) * nd$points_per_logit
    add_axis(nd$feats[j], n_rows - i,
              abs(nd$raw_beta[j]) * diff(limits) * nd$points_per_logit,
              as.numeric(points), tick_labels(ticks))
  }
  ticks <- pretty(c(0, nd$total_max), n = 5L)
  ticks <- ticks[ticks >= 0 & ticks <= nd$total_max]
  if (length(ticks) < 2L) ticks <- c(0, nd$total_max)
  add_axis("Total Points", 2, 100, ticks / nd$total_max * 100, tick_labels(ticks))
  probability_ticks <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3,
                          0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995, 0.999)
  x <- (stats::qlogis(probability_ticks) - nd$eta_min) * nd$points_per_logit / nd$total_max * 100
  keep <- which(is.finite(x) & x >= 0 & x <= 100)
  probability_ticks <- probability_ticks[keep]
  x <- x[keep]
  if (length(x) < 2L) {
    x <- seq(0, 100, length.out = 5L)
    probability_ticks <- stats::plogis(nd$eta_min + x / 100 * nd$total_max / nd$points_per_logit)
  } else {
    spaced <- 1L
    for (j in seq.int(2L, length(x))) {
      if (x[j] - x[tail(spaced, 1L)] >= 12) spaced <- c(spaced, j)
    }
    x <- x[spaced]
    probability_ticks <- probability_ticks[spaced]
  }
  add_axis(paste(model$positive_class, "Probability"), 1, 100, x,
             format(signif(probability_ticks, 3L), trim = TRUE))
  axes <- .frequency_bind_rows(axis_rows)
  ticks <- .frequency_bind_rows(tick_rows)
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = axes,
                           ggplot2::aes(x = X, xend = X_end, y = Y, yend = Y), linewidth = 0.4) +
    ggplot2::geom_segment(data = ticks,
                           ggplot2::aes(x = X, xend = X, y = Y, yend = Y - 0.12), linewidth = 0.35) +
    ggplot2::geom_text(data = ticks, ggplot2::aes(x = X, y = Y - 0.34, label = Label), size = 2.8) +
    ggplot2::geom_text(data = axes, ggplot2::aes(x = -7, y = Y, label = Label),
                        hjust = 1, size = 3.2) +
    ggplot2::coord_cartesian(xlim = c(-55, 105), ylim = c(0.25, n_rows + 0.5), clip = "off") +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.margin = ggplot2::margin(12, 18, 12, 18),
                    plot.title = ggplot2::element_text(face = "bold"),
                    plot.subtitle = ggplot2::element_text(size = 9)) +
    ggplot2::labs(title = "FOREST Best-Model Nomogram", subtitle = subtitle)
}


# Accept either outfile itself or its parent; always save inside a folder named outfile.
.frequency_outfile_dir <- function(output_dir) {
  stopifnot(is.character(output_dir), length(output_dir) == 1L,
            !is.na(output_dir), nzchar(output_dir))
  path <- sub("[/\\\\]+$", "", output_dir)
  if (!nzchar(path)) path <- output_dir
  if (identical(basename(path), "outfile")) path else file.path(path, "outfile")
}

.frequency_roc_filenames <- function() {
  c("FOREST_ROC.pdf", "FOREST_ROC.png", "roc_results.rds",
    "ROC_auc_per_iter.csv", "ROC_auc_summary.csv", "ROC_curve_summary.csv",
    "ROC_curve_per_iter.csv", "ROC_record_check.csv")
}

.frequency_figure_filenames <- function() {
  c(.frequency_roc_filenames(), "FOREST_best_model_NORM.pdf", "FOREST_best_model_NORM.png",
    "best_model.rds", "best_model_summary.csv", "best_model_coefficients.csv")
}

.frequency_report_roc_records <- function(roc_results) {
  check <- roc_results$record_check
  print(check[, c("Parameter", "Set", "N_Planned_Repeats", "N_Recorded_Repeats",
                   "N_Failed_Repeats", "AUC_SD", "Max_TPR_Interval_Width"), drop = FALSE],
          row.names = FALSE)
  for (i in seq_len(nrow(check))) {
    if (nzchar(check$Interval_Note[i])) {
      message("ROC interval, Parameter ", check$Parameter[i], ", ", check$Set[i], ": ",
                check$Interval_Note[i])
    }
    if (is.na(check$N_Expected_Successful_Repeats[i]) &&
        !is.na(check$N_Planned_Repeats[i]) &&
        check$N_Recorded_Repeats[i] < check$N_Planned_Repeats[i]) {
      message("Only ", check$N_Recorded_Repeats[i], " of ", check$N_Planned_Repeats[i],
                " planned repeats are present for Parameter ", check$Parameter[i], ", ", check$Set[i],
                ". Use result= to distinguish failed attempts from missing saved records.")
    }
  }
  invisible(check)
}

.frequency_write_roc_outputs <- function(roc_results, roc_plot, output_dir) {
  tables <- list(ROC_auc_per_iter = roc_results$auc_per_iter,
                 ROC_auc_summary = roc_results$auc_summary,
                 ROC_curve_summary = roc_results$plot_data,
                 ROC_curve_per_iter = roc_results$iteration_curves,
                 ROC_record_check = roc_results$record_check)
  for (nm in names(tables)) {
    utils::write.csv(tables[[nm]], file.path(output_dir, paste0(nm, ".csv")),
                       row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  }
  saveRDS(roc_results, file.path(output_dir, "roc_results.rds"))
  n_parameters <- length(unique(roc_results$avg_roc$Parameter))
  width <- if (n_parameters == 1L) 5 else 5 * min(3L, n_parameters)
  height <- if (n_parameters == 1L) 5 else 5 * ceiling(n_parameters / min(3L, n_parameters))
  for (extension in c("pdf", "png")) {
    ggplot2::ggsave(file.path(output_dir, paste0("FOREST_ROC.", extension)),
                     plot = roc_plot, width = width, height = height,
                     units = "in", dpi = 300, bg = "white", limitsize = FALSE)
  }
  file.path(output_dir, .frequency_roc_filenames())
}

# Save the supplied reference ROC workflow directly from results_df.
# No saved model, input CSV, or previous result object is needed for this function.
save_frequency_roc <- function(
    pred_data = NULL, output_dir = "outfile",
    fpr_grid = seq(0, 1, length.out = 101L), positive_class = "Tumor",
    show_iterations = TRUE, overwrite = FALSE, install_missing_packages = TRUE,
    result = NULL, expected_n_repeats = NULL, parameter_id = NULL) {
  stopifnot(is.logical(overwrite), length(overwrite) == 1L, !is.na(overwrite))
  .frequency_require_packages(install_missing_packages,
                              pkgs = c("pROC", "ggplot2", "dplyr", "tidyr", "purrr", "tibble"))
  output_dir <- .frequency_outfile_dir(output_dir)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create output directory: ", output_dir, call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  paths <- file.path(output_dir, .frequency_roc_filenames())
  existing <- paths[file.exists(paths)]
  if (length(existing) && !overwrite) {
    stop("ROC outputs already exist. Choose another directory or set overwrite=TRUE: ",
           paste(basename(existing), collapse = ", "), call. = FALSE)
  }
  records <- NULL
  if (!is.null(result)) {
    records <- frequency_roc_records(result, fpr_grid, parameter_id)
    if (is.null(pred_data)) pred_data <- records$results_df
    positive_class <- records$positive_class
    if (is.null(expected_n_repeats)) expected_n_repeats <- records$n_repeats
  }
  if (is.null(pred_data)) stop("Provide result= or pred_data=results_df for ROC record retrieval.", call. = FALSE)
  roc_results <- calculate_resampled_roc(
    pred_data, fpr_grid, positive_class,
    reference_auc = if (is.null(records)) NULL else records$auc_df,
    reference_curves = if (is.null(records)) NULL else records$curve_df,
    expected_n_repeats = expected_n_repeats
  )
  .frequency_report_roc_records(roc_results)
  roc_plot <- plot_frequency_resampled_roc(roc_results, show_iterations)
  .frequency_write_roc_outputs(roc_results, roc_plot, output_dir)
  invisible(list(roc_results = roc_results, roc_plot = roc_plot,
                   record_check = roc_results$record_check, records = records, output_files = paths))
}

# Export plots from an existing analysis without repeating any split, fitting, or validation.
save_frequency_figures <- function(
    result, data, results_df = NULL, output_dir = "outfile",
    fpr_grid = seq(0, 1, length.out = 101L), show_iterations = TRUE,
    overwrite = FALSE, install_missing_packages = TRUE) {
  stopifnot(is.logical(overwrite), length(overwrite) == 1L, !is.na(overwrite))
  .frequency_require_packages(install_missing_packages)
  output_dir <- .frequency_outfile_dir(output_dir)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create output directory: ", output_dir, call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  paths <- file.path(output_dir, .frequency_figure_filenames())
  existing <- paths[file.exists(paths)]
  if (length(existing) && !overwrite) {
    stop("Figure outputs already exist. Choose another output directory or set overwrite=TRUE: ",
           paste(basename(existing), collapse = ", "), call. = FALSE)
  }
  if (is.null(results_df)) results_df <- result$predictions
  if (!is.data.frame(results_df) || !nrow(results_df)) {
    stop("No completed repeat predictions are available for ROC export.", call. = FALSE)
  }
  expected <- result$holdout_results
  expected <- expected[!is.na(expected$Error) & expected$Error == "" &
                         is.finite(expected[["TRN-AUC"]]) & is.finite(expected[["VAL-AUC"]]), , drop = FALSE]
  saved <- unique(results_df[, c("Parameter", "Repeat", "Set"), drop = FALSE])
  wanted <- do.call(rbind, lapply(c("TRN", "VAL"), function(s) {
    data.frame(expected[, c("Parameter", "Repeat"), drop = FALSE], Set = s)
  }))
  group_key <- function(z) paste(z$Parameter, z$Repeat, z$Set, sep = "::")
  if (!setequal(group_key(saved), group_key(wanted))) {
    stop("Predictions must include both sets for every successful repeat. Use keep_all_repeats=TRUE.", call. = FALSE)
  }
  saved_curves <- .frequency_bind_rows(frequency_result_lists(result, fpr_grid)$curve_list)
  roc_results <- calculate_resampled_roc(
    results_df, fpr_grid, result$config$positive_class,
    reference_auc = result$holdout_results, reference_curves = saved_curves,
    expected_n_repeats = result$config$n_repeats
  )
  .frequency_report_roc_records(roc_results)
  # Verify AUCs from the saved predictions against the original fitting results.
  aa <- roc_results$auc_per_iter
  ii <- match(paste(aa$Parameter, aa$Repeat, sep = "::"),
                paste(expected$Parameter, expected$Repeat, sep = "::"))
  original_auc <- ifelse(aa$Set == "TRN", expected[["TRN-AUC"]][ii], expected[["VAL-AUC"]][ii])
  stopifnot(!anyNA(ii), max(abs(aa$auc - original_auc)) < 1e-10)
  roc_plot <- plot_frequency_resampled_roc(roc_results, show_iterations)
  best <- select_frequency_best_model(result)
  best_rows <- results_df[results_df$Parameter == best$model$Parameter &
                           results_df$Repeat == best$model$Repeat & results_df$Set == "TRN", , drop = FALSE]
  training <- data[match(best_rows$ID, data$ID), , drop = FALSE]
  stopifnot(!anyNA(match(best_rows$ID, data$ID)))
  if (max(abs(predict_frequency_model(best$model, training) - best_rows$Probability)) >= 1e-8) {
    stop("The supplied data do not reproduce the saved model's training predictions; use the original input.", call. = FALSE)
  }
  nomogram_data <- frequency_nomogram_data(best$model, data)
  best$model$nomogram <- nomogram_data
  norm_plot <- plot_frequency_nomogram(best$model, nomogram_data)
  tables <- list(best_model_summary = best$summary,
                 best_model_coefficients = nomogram_data$coefficients)
  for (nm in names(tables)) {
    utils::write.csv(tables[[nm]], file.path(output_dir, paste0(nm, ".csv")),
                       row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  }
  saveRDS(best$model, file.path(output_dir, "best_model.rds"))
  .frequency_write_roc_outputs(roc_results, roc_plot, output_dir)
  norm_height <- max(4.5, 2 + 0.55 * (length(nomogram_data$active) + 3L))
  for (extension in c("pdf", "png")) {
    ggplot2::ggsave(file.path(output_dir, paste0("FOREST_best_model_NORM.", extension)),
                     plot = norm_plot, width = 9.5, height = norm_height,
                     units = "in", dpi = 300, bg = "white", limitsize = FALSE)
  }
  invisible(list(roc_results = roc_results, roc_plot = roc_plot,
                   best_model = best$model, best_model_summary = best$summary,
                   nomogram_data = nomogram_data, norm_plot = norm_plot,
                   output_files = paths))
}



# ============================ Input and execution entry point ============================
.frequency_bind_rows <- function(z) {
  out <- do.call(rbind, z)
  if (is.null(out)) return(data.frame())
  rownames(out) <- NULL
  out
}

read_frequency_input <- function(input_file) {
  if (!is.character(input_file) || length(input_file) != 1L ||
      is.na(input_file) || !file.exists(input_file)) {
    stop("Input file not found. Place input.csv in getwd() or set input_file to its full path.",
         call. = FALSE)
  }
  d <- utils::read.csv(
    input_file, header = TRUE, row.names = NULL, check.names = FALSE,
    stringsAsFactors = FALSE, colClasses = "character", fill = FALSE,
    na.strings = c("", "NA", "NaN"), fileEncoding = "UTF-8",
    comment.char = ""
  )
  if (nrow(d) < 1L || ncol(d) < 4L) {
    stop("The CSV needs at least 1 sample row and 4 columns: ID, Group, and at least 2 features.",
         call. = FALSE)
  }
  names(d)[1L] <- sub(intToUtf8(65279L), "", names(d)[1L], fixed = TRUE)
  if (anyNA(names(d)) || any(!nzchar(trimws(names(d)))) ||
      anyDuplicated(names(d))) {
    stop("CSV column names must be nonempty and unique.", call. = FALSE)
  }
  names(d)[1:2] <- c("ID", "Group")
  if (any(c("ID", "Group") %in% names(d)[-(1:2)])) {
    stop("Feature columns cannot be named ID or Group.", call. = FALSE)
  }
  d$ID <- trimws(d$ID)
  d$Group <- trimws(d$Group)
  if (anyNA(d$ID) || any(!nzchar(d$ID)) || anyDuplicated(d$ID)) {
    stop("Sample IDs in column 1 must be nonempty and unique; each sample must occupy a single row.",
         call. = FALSE)
  }
  if (anyNA(d$Group) || !setequal(unique(d$Group), c("Tumor", "Control"))) {
    stop("Column 2 must contain both Tumor and Control, with exact case and no missing groups.",
         call. = FALSE)
  }
  for (j in seq.int(3L, ncol(d))) {
    raw <- trimws(d[[j]])
    missing <- is.na(raw) | raw %in% c("", "NA", "NaN")
    value <- suppressWarnings(as.numeric(raw))
    bad <- !missing & (is.na(value) | !is.finite(value))
    if (any(bad)) {
      k <- which(bad)[1L]
      stop(sprintf("Feature %s contains a nonnumeric or infinite value for sample %s: %s",
                   names(d)[j], d$ID[k], raw[k]), call. = FALSE)
    }
    value[missing] <- NA_real_
    d[[j]] <- value
  }
  d
}

.frequency_require_packages <- function(
    install_missing = TRUE,
    pkgs = c("glmnet", "pROC", "ggplot2", "dplyr", "tidyr", "purrr", "tibble")) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) && install_missing) {
    message("Installing missing R packages: ", paste(missing, collapse = ", "))
    attempt <- tryCatch({
      utils::install.packages(missing, repos = "https://cloud.r-project.org",
                              dependencies = NA)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(attempt)) message("Automatic installation did not complete: ", attempt)
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". Install these packages before running the script.",
         call. = FALSE)
  }
  if (utils::packageVersion("ggplot2") < "3.4.0") {
    stop("ROC and nomogram exports need ggplot2 >= 3.4.0. Update it with install.packages(\"ggplot2\").",
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
      stop("Setting ", k, " must be a single nonnegative integer.", call. = FALSE)
    }
  }
  for (k in c("ensure_oof_coverage", "require_complete_oof",
               "install_missing_packages", "save_models", "overwrite", "verbose",
               "roc_show_iterations")) {
    v <- cfg[[k]]
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("Setting ", k, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  for (k in c("train_fraction", "lambda_ratio_min")) {
    v <- cfg[[k]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0 ||
        v >= 1) {
      stop("Setting ", k, " is outside the allowed range. Use one fixed configuration per run.",
           call. = FALSE)
    }
  }
  if (!is.numeric(cfg$alpha_grid) || length(cfg$alpha_grid) < 1L ||
      any(!is.finite(cfg$alpha_grid)) || anyDuplicated(cfg$alpha_grid) ||
      any(cfg$alpha_grid < 0.5 | cfg$alpha_grid > 1)) {
    stop("alpha_grid must contain unique numeric values in [0.5,1]; the default is seq(0.5,1,by=0.05).",
         call. = FALSE)
  }
  if (!identical(cfg$positive_class, "Tumor")) {
    stop("This entry point requires Tumor as the positive class to preserve the meaning of Tumor_Probability.",
         call. = FALSE)
  }
  if (cfg$TOP_n < 2L || cfg$TOP_n > ncol(d) - 2L ||
      cfg$n_repeats < 2L ||
      cfg$n_lambda < 2L || cfg$cv_nfolds < 3L || cfg$fit_maxit < 1L ||
      cfg$plot_iter < 1L || cfg$plot_iter > cfg$n_repeats) {
    stop("Check integer settings: TOP_n cannot exceed the feature count; repeats and lambda points must be at least 2; CV folds must be at least 3.",
         call. = FALSE)
  }
  if (!is.null(cfg$pvalue_cutoff) &&
      (!is.numeric(cfg$pvalue_cutoff) || length(cfg$pvalue_cutoff) != 1L ||
       !is.finite(cfg$pvalue_cutoff) || cfg$pvalue_cutoff <= 0 || cfg$pvalue_cutoff >= 1)) {
    stop("pvalue_cutoff must be a number in (0,1) or NULL.", call. = FALSE)
  }
  if (!is.numeric(cfg$coef_tol) || length(cfg$coef_tol) != 1L ||
      !is.finite(cfg$coef_tol) || cfg$coef_tol < 0) {
    stop("coef_tol must be a single nonnegative finite number.", call. = FALSE)
  }
  if (!is.character(cfg$output_dir) || length(cfg$output_dir) != 1L ||
      is.na(cfg$output_dir) || !nzchar(cfg$output_dir)) {
    stop("output_dir must be a nonempty directory path.", call. = FALSE)
  }
  if (!is.character(cfg$diff_method) || length(cfg$diff_method) != 1L ||
      !cfg$diff_method %in% c("wilcox", "welch")) {
    stop("diff_method must be wilcox or welch.",
         call. = FALSE)
  }
  if (!is.numeric(cfg$fpr_grid) || length(cfg$fpr_grid) < 2L ||
      any(!is.finite(cfg$fpr_grid)) || any(cfg$fpr_grid < 0 | cfg$fpr_grid > 1) ||
      any(diff(cfg$fpr_grid) <= 0) || cfg$fpr_grid[1L] != 0 || tail(cfg$fpr_grid, 1L) != 1) {
    stop("fpr_grid must be a strictly increasing numeric grid starting at 0 and ending at 1.", call. = FALSE)
  }
  n <- as.numeric(table(factor(d$Group, levels = c("Control", "Tumor"))))
  tr <- floor(n * cfg$train_fraction)
  va <- n - tr
  if (min(n) < 5L || min(tr) < 3L || min(va) < 2L) {
    stop("Insufficient sample counts or an unsuitable fraction: each class needs at least 5 total, 3 training, and 2 validation samples.",
         call. = FALSE)
  }
  if (cfg$ensure_oof_coverage && any(cfg$n_repeats * va < n)) {
    stop(sprintf("Covering all samples requires at least %d outer repeats; increase n_repeats.",
                 max(ceiling(n / va))), call. = FALSE)
  }
}

# Supply only configuration overrides; NULL is preserved so p-value filtering can be disabled.
run_frequency_pipeline <- function(config = list()) {
  defaults <- frequency_default_config()
  if (!is.list(config) || is.null(names(config)) && length(config) ||
      anyDuplicated(names(config))) {
    stop("config must be a named list with unique parameter names.", call. = FALSE)
  }
  unknown <- setdiff(names(config), names(defaults))
  if (length(unknown)) {
    stop("Unknown configuration names: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  cfg <- utils::modifyList(defaults, config, keep.null = TRUE)
  d <- read_frequency_input(cfg$input_file)
  .frequency_validate_config(cfg, d)
  cfg$output_dir <- .frequency_outfile_dir(cfg$output_dir)

  files <- c(
    "oof_predictions.csv", "oof_predictions_summary.csv",
    "results_df.csv", "model_df.csv", "auc_df.csv", "curve_df.csv",
    "indexfilterresult.csv", "selected_features.csv", "feature_frequency.csv",
    "differential_results.csv", "screen_diagnostics.csv", "final_diagnostics.csv",
    "indexfilter_errors.csv", "split_membership.csv", "split_adjustments.csv",
    "run_config.R", "session_info.txt", "run_status.txt", "run_log.txt",
    .frequency_figure_filenames()
  )
  if (cfg$save_models) files <- c(files, "frequency_analysis.rds")
  if (!dir.exists(cfg$output_dir) &&
      !dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create output directory: ", cfg$output_dir, call. = FALSE)
  }
  cfg$input_file <- normalizePath(cfg$input_file, winslash = "/", mustWork = TRUE)
  cfg$output_dir <- normalizePath(cfg$output_dir, winslash = "/", mustWork = TRUE)
  paths <- file.path(cfg$output_dir, files)
  if (cfg$input_file %in% paths) stop("The input file cannot share a path with an output file.")
  existing <- paths[file.exists(paths)]
  if (length(existing) && !cfg$overwrite) {
    stop("Named result files already exist. Change output_dir or explicitly set overwrite=TRUE: ",
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
      "Status: ERROR", "Error: Analysis or export was interrupted; see the R console error and run_log.txt."
    )), silent = TRUE)
  }, add = TRUE)
  on.exit(utils::capture.output(utils::sessionInfo(),
                               file = file.path(cfg$output_dir, "session_info.txt")), add = TRUE)
  analysis <- tryCatch(withCallingHandlers(
    run_frequency_holdouts(
      data = d, positive_class = cfg$positive_class,
      train_fraction_grid = cfg$train_fraction, TOP_n = cfg$TOP_n,
      alpha_grid = cfg$alpha_grid, lambda_ratio_min_grid = cfg$lambda_ratio_min,
      n_repeats = cfg$n_repeats,
      n_lambda = cfg$n_lambda, pvalue_cutoff = cfg$pvalue_cutoff,
      diff_method = cfg$diff_method, cv_nfolds = cfg$cv_nfolds,
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

  # OOF probabilities come only from outer VAL predictions; exclude TRN and training-CV fit.preval.
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
    # If all splits fail, frequency/feature tables may lack columns; still save an identifiable empty table.
    if (!ncol(tab)) tab <- data.frame(Note = character())
    utils::write.csv(tab, file.path(cfg$output_dir, paste0(nm, ".csv")),
                     row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  }
  missing_ids <- oof_summary$ID[!oof_summary$Has_OOF]
  n_success <- sum(!nzchar(auc_df$Error))
  figures <- list(roc_results = NULL, roc_plot = NULL, best_model = NULL, norm_plot = NULL,
                  output_files = character())
  figure_error <- NULL
  if (n_success > 0L) {
    if (cfg$verbose) message("Exporting all-repeat ROC and best-model nomogram to ", cfg$output_dir)
    figures <- tryCatch(
      save_frequency_figures(
        result = analysis, data = d, results_df = results_df,
        output_dir = cfg$output_dir, fpr_grid = cfg$fpr_grid,
        show_iterations = cfg$roc_show_iterations, overwrite = cfg$overwrite,
        install_missing_packages = cfg$install_missing_packages
      ),
      error = function(e) {
        figure_error <<- conditionMessage(e)
        cat("Figure export error: ", figure_error, "\n", file = log_path, append = TRUE)
        list(roc_results = NULL, roc_plot = NULL, best_model = NULL, norm_plot = NULL,
             output_files = character())
      }
    )
  }
  bundle <- list(config = cfg, analysis = analysis,
                 results_df = results_df, model_df = model_df,
                 auc_df = auc_df, curve_df = curve_df,
                 oof_predictions = oof, oof_predictions_summary = oof_summary,
                 figures = figures, roc_results = figures$roc_results,
                 roc_plot = figures$roc_plot, best_model = figures$best_model,
                 norm_plot = figures$norm_plot, figure_error = figure_error,
                 output_files = paths)
  if (cfg$save_models) saveRDS(bundle, file.path(cfg$output_dir, "frequency_analysis.rds"))
  write_status(c(
    if (!is.null(figure_error)) "Status: ERROR_FIGURE_EXPORT" else if (length(missing_ids)) "Status: INCOMPLETE_OOF" else "Status: COMPLETE",
    paste("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Input_MD5:", unname(tools::md5sum(cfg$input_file))),
    paste("Samples:", nrow(d)), paste("Outer_Repeats:", cfg$n_repeats),
    paste("Alpha_Grid:", paste(cfg$alpha_grid, collapse = ";")),
    paste("Requested_Frequency_Models_Per_Holdout:", length(cfg$alpha_grid) * cfg$n_lambda),
    paste("Successful_Models:", n_success),
    paste("Failed_Models:", nrow(auc_df) - n_success),
    paste("Figure_Export:", if (!is.null(figure_error)) "ERROR" else if (n_success) "COMPLETE" else "SKIPPED_NO_SUCCESSFUL_MODELS"),
    paste("Figure_Error:", if (is.null(figure_error)) "" else figure_error),
    paste("Samples_With_OOF:", sum(oof_summary$Has_OOF)),
    paste("Samples_Without_OOF:", length(missing_ids)),
    paste("Missing_IDs:", paste(missing_ids, collapse = ";"))
  ))
  status_finalized <- TRUE
  if (cfg$verbose) {
    message(sprintf("Results saved to %s; successful models %d/%d; OOF predictions cover %d/%d samples.",
                    cfg$output_dir, n_success, nrow(auc_df),
                    sum(oof_summary$Has_OOF), nrow(d)))
  }
  if (!is.null(figure_error)) {
    stop("Modeling tables were saved, but figure export failed: ", figure_error,
         ". See run_status.txt and run_log.txt.", call. = FALSE)
  }
  if (length(missing_ids) && cfg$require_complete_oof) {
    stop(sprintf(
      "%d samples have no valid OOF probability. Results are saved; see oof_predictions_summary.csv and indexfilter_errors.csv. Training probabilities are not used as substitutes.",
      length(missing_ids)
    ), call. = FALSE)
  }
  invisible(bundle)
}

# Run directly with source() or Rscript; the file does not call other project scripts.
.frequency_is_cli <- sys.nframe() == 0L
.frequency_args <- if (.frequency_is_cli) commandArgs(trailingOnly = TRUE) else character()
if (.frequency_is_cli && any(.frequency_args %in% c("--help", "-h"))) {
  cat("Usage: Rscript FOREST.R [input.csv] [output_directory]\n")
  cat("Default input: input.csv. Default output: ./outfile. A supplied parent directory receives an outfile subfolder. Set parameters through frequency_config or the defaults at the top of the script.\n")
} else if (!isTRUE(getOption("frequency_oof.functions_only", FALSE))) {
  .frequency_cfg <- if (exists("frequency_config", inherits = FALSE)) frequency_config else list()
  # Define n_repeats in the calling session; an explicit frequency_config value takes precedence.
  if (is.list(.frequency_cfg) && !("n_repeats" %in% names(.frequency_cfg)) &&
      exists("n_repeats", inherits = FALSE)) {
    .frequency_cfg$n_repeats <- n_repeats
  }
  if (length(.frequency_args) > 2L) stop("At most two positional arguments are accepted: input CSV and output directory.")
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
  roc_results <- frequency_run$roc_results
  roc_plot <- frequency_run$roc_plot
  best_model <- frequency_run$best_model
  norm_plot <- frequency_run$norm_plot
}

