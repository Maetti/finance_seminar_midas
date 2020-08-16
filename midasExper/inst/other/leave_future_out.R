library("loo")
library("brms")
library("bayesplot")
library("ggplot2")
color_scheme_set("brightblue")
theme_set(theme_default())

CHAINS <- 4
SEED <- 5838296
set.seed(SEED)



N <- length(LakeHuron)
df <- data.frame(
      y = as.numeric(LakeHuron),
      year = as.numeric(time(LakeHuron)),
      time = 1:N
)

# save plot labels to reuse them
plot_labs <- labs(
      y = "Water Level (ft)",
      x = "Year",
      title = "Water Level in Lake Huron (1875-1972)"
)

ggplot(df, aes(x = year, y = y)) +
      geom_point(size = 1) +
      plot_labs


control <- list(adapt_delta = 0.99)
fit <- brm(
      y ~ ar(time, p = 4),
      data = df,
      prior = prior(normal(0, 0.5), class = "ar"),
      control = control,
      seed = SEED,
      chains = CHAINS, 
      cores = 4
)


preds <- posterior_predict(fit)

preds <- cbind(
      Estimate = colMeans(preds),
      Q5 = apply(preds, 2, quantile, probs = 0.05),
      Q95 = apply(preds, 2, quantile, probs = 0.95)
)

ggplot(cbind(df, preds), aes(x = year, y = Estimate)) +
      geom_smooth(aes(ymin = Q5, ymax = Q95), stat = "identity", size = 0.5) +
      geom_point(aes(y = y)) +
      labs(subtitle = "Mean (blue) and 90% predictive intervals (gray) vs. observed data (black)") +
      plot_labs

L <- 90

loo_cv <- loo(log_lik(fit)[, (L + 1):N])
print(loo_cv)



#   ____________________________________________________________________________
#   Exact 1 Step ahead LFO                                                  ####

loglik_exact <- matrix(nrow = nsamples(fit), ncol = N)
for (i in L:(N - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      fit_i <- update(fit, newdata = df_past, recompile = FALSE, cores = 4)
      loglik_exact[, i + 1] <- log_lik(fit_i, newdata = df_oos, oos = oos)[, oos]
}


# some helper functions we'll use throughout

# more stable than log(sum(exp(x))) 
log_sum_exp <- function(x) {
      max_x <- max(x)
      max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
      log_sum_exp(x) - log(length(x))
}

# compute log of raw importance ratios
# sums over observations *not* over posterior samples
sum_log_ratios <- function(loglik, ids = NULL) {
      if (!is.null(ids)) loglik <- loglik[, ids, drop = FALSE]
      rowSums(loglik)
}

# for printing comparisons later
rbind_print <- function(...) {
      round(rbind(...), digits = 2)
}

exact_elpds_1sap <- apply(loglik_exact, 2, log_mean_exp)
exact_elpd_1sap <- c(ELPD = sum(exact_elpds_1sap[-(1:L)]))

rbind_print(
      "LOO" = loo_cv$estimates["elpd_loo", "Estimate"],
      "LFO" = exact_elpd_1sap
)



#   ____________________________________________________________________________
#   Approximate 1-step ahead                                                ####

k_thres <- 0.7


approx_elpds_1sap <- rep(NA, N)

# initialize the process for i = L
past <- 1:L
oos <- L + 1
df_past <- df[past, , drop = FALSE]
df_oos <- df[c(past, oos), , drop = FALSE]
fit_past <- update(fit, newdata = df_past, recompile = FALSE, cores = 4)
loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
approx_elpds_1sap[L + 1] <- log_mean_exp(loglik[, oos])

# iterate over i > L
i_refit <- L
refits <- L
ks <- NULL
maPred <- matrix(0, nrow = 4000, ncol = 7)
j <- 1
for (i in (L + 1):(N - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
      
      logratio <- sum_log_ratios(loglik, (i_refit + 1):i)
      psis_obj <- suppressWarnings(psis(logratio))
      k <- pareto_k_values(psis_obj)
      ks <- c(ks, k)
      if (k > k_thres) {
            # refit the model based on the first i observations
            i_refit <- i
            refits <- c(refits, i)
            fit_past <- update(fit_past, newdata = df_past, recompile = FALSE)
            loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
            approx_elpds_1sap[i + 1] <- log_mean_exp(loglik[, oos])
      } else {
            lw <- weights(psis_obj, normalize = TRUE)[, 1]
            approx_elpds_1sap[i + 1] <- log_sum_exp(lw + loglik[, oos])
      }
      
      ## posterior predict
      maPred[, j] <- posterior_predict(fit_past, newdata = df_oos[oos, c(1, 3)])
      j <- j + 1
}

# preds <- cbind(
#       Estimate = colMeans(maPred),
#       Q5 = apply(maPred, 2, quantile, probs = 0.05),
#       Q95 = apply(maPred, 2, quantile, probs = 0.95)
# )
# preds <- as.data.frame(preds)
# preds$year <- 1966:1972
# 
# ggplot(df, aes(x = year, y = y)) +
#       geom_line() +
#       stat_smooth(data = preds, aes(x = year, y = Estimate, ymin = Q5, ymax = Q95), stat = "identity", size = 0.5) +
#       labs(subtitle = "Mean (blue) and 90% predictive intervals (gray) vs. observed data (black)") +
#       plot_labs
# 
# ggplot(df, aes(x = year, y = y)) + 
#       geom_line(colour='blue') +
#       geom_smooth(aes(x=year, y=Estimate, ymax=Q95, ymin=Q5), 
#             colour='red', data=preds, stat='identity')
# 
# ggplot() + 
#       geom_line(data=df1, aes(x = time, y = M, color = isin)) + 
#       stat_smooth(data=df2, aes(x = time, y = M, color = isin))


## analysis
plot_ks <- function(ks, ids, thres = 0.6) {
      dat_ks <- data.frame(ks = ks, ids = ids)
      ggplot(dat_ks, aes(x = ids, y = ks)) +
            geom_point(aes(color = ks > thres), shape = 3, show.legend = FALSE) +
            geom_hline(yintercept = thres, linetype = 2, color = "red2") +
            scale_color_manual(values = c("cornflowerblue", "darkblue")) +
            labs(x = "Data point", y = "Pareto k") +
            ylim(-0.5, 1.5)
}

cat("Using threshold ", k_thres,
      ", model was refit ", length(refits),
      " times, at observations", refits)

## compare approximation with exact elpd_1sap
approx_elpd_1sap <- sum(approx_elpds_1sap, na.rm = TRUE)
rbind_print(
      "approx LFO" = approx_elpd_1sap,
      "exact LFO" = exact_elpd_1sap
)


dat_elpd <- data.frame(
      approx_elpd = approx_elpds_1sap,
      exact_elpd = exact_elpds_1sap
)

ggplot(dat_elpd, aes(x = approx_elpd, y = exact_elpd)) +
      geom_abline(color = "gray30") +
      geom_point(size = 2) +
      labs(x = "Approximate ELPDs", y = "Exact ELPDs")
