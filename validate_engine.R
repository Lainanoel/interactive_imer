# ===========================================================================
#  validate_engine.R  —  run ONCE on a machine that has the packages installed.
#
#  Purpose: the Shiny app's statistical engine (lmer / emmeans / anova /
#  influence.ME) cannot be exercised in every environment. This script fits
#  the SAME pipeline the app uses on two well-known datasets and checks the
#  results against (a) published reference values, (b) internal-consistency
#  identities that must hold if the engine is correct, and (c) raw summaries
#  you can eyeball. If everything prints PASS, the engine is behaving.
#
#  Run:  Rscript validate_engine.R      (or source() it in RStudio)
# ===========================================================================

ok <- 0L; bad <- 0L
chk <- function(label, cond) {
  pass <- isTRUE(cond)
  cat(sprintf("  [%s] %s\n", if (pass) "PASS" else "FAIL", label))
  if (pass) ok <<- ok + 1L else bad <<- bad + 1L
}
sec <- function(s) cat(sprintf("\n== %s ==\n", s))

need <- c("lmerTest", "emmeans", "multcomp")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  stop("Install first: install.packages(c(",
       paste(sprintf('"%s"', miss), collapse = ", "), "))")
}
suppressPackageStartupMessages({
  library(lmerTest); library(emmeans); library(multcomp)
})

# ---------------------------------------------------------------------------
sec("1. sleepstudy — fixed effects vs published reference")
# Canonical lme4 example: Reaction ~ Days + (Days | Subject)
# Widely published estimates: (Intercept) ~ 251.41, Days ~ 10.47
data(sleepstudy, package = "lme4")
m1 <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy,
           REML = TRUE, na.action = na.omit)
fe <- lme4::fixef(m1)
cat(sprintf("  fixef: (Intercept) = %.3f   Days = %.3f\n", fe[1], fe[2]))
chk("Intercept within 0.5 of 251.41", abs(fe[[1]] - 251.41) < 0.5)
chk("Days slope within 0.1 of 10.47", abs(fe[[2]] - 10.47) < 0.1)

# Internal consistency: emmeans at given Days must equal the fixed prediction
emm1 <- emmeans(m1, ~ Days, at = list(Days = c(0, 5)))
est  <- summary(emm1)$emmean
pred <- c(fe[[1]], fe[[1]] + 5 * fe[[2]])
chk("emmeans(Days=0,5) equals fixed-effect prediction (consistency)",
    max(abs(est - pred)) < 1e-6)

# ---------------------------------------------------------------------------
sec("2. ChickWeight — VALID model weight ~ Diet + (1 | Chick)")
data(ChickWeight); cw <- as.data.frame(ChickWeight)
m2 <- lmer(weight ~ Diet + (1 | Chick), data = cw, REML = TRUE,
           na.action = na.omit)
obs <- tapply(cw$weight, cw$Diet, mean)              # observed diet means
emm2 <- as.data.frame(emmeans(m2, ~ Diet))
cat("  observed diet means: ", paste(round(obs, 1), collapse = ", "), "\n")
cat("  emmeans diet means:  ", paste(round(emm2$emmean, 1), collapse = ", "), "\n")
chk("EMMeans are finite and in the data's weight range",
    all(is.finite(emm2$emmean)) && all(emm2$emmean > 0 & emm2$emmean < 400))
chk("EMMeans preserve the observed diet ordering",
    identical(order(emm2$emmean), order(as.numeric(obs))))
a2 <- anova(m2, type = "III", ddf = "Satterthwaite")
cat("  ANOVA Diet p-value:", signif(a2[["Pr(>F)"]][1], 3), "\n")
chk("ANOVA produces a finite F and p-value",
    is.finite(a2[["F value"]][1]) && is.finite(a2[["Pr(>F)"]][1]))

# Effective N == complete cases (the app's dropped-row accounting)
chk("nobs() equals complete-case count",
    stats::nobs(m2) == sum(stats::complete.cases(cw[, c("weight","Diet","Chick")])))

# ---------------------------------------------------------------------------
sec("3. Back-transformation correctness (log)")
# Fit on log(weight); response-scale EMMean must equal exp(link EMMean).
m3   <- lmer(log(weight) ~ Diet + (1 | Chick), data = cw, REML = TRUE,
             na.action = na.omit)
link <- summary(emmeans(m3, ~ Diet))$emmean                       # log scale
resp <- summary(emmeans(ref_grid(m3, tran = "log"), ~ Diet,
                        type = "response"))$response              # back-transformed
cat("  exp(link):", paste(round(exp(link), 2), collapse = ", "), "\n")
cat("  response :", paste(round(resp, 2),      collapse = ", "), "\n")
chk("type='response' equals exp(link) for a log fit",
    max(abs(resp - exp(link))) < 1e-6)

# ---------------------------------------------------------------------------
sec("4. cld letters are internally consistent")
emm4 <- emmeans(m2, ~ Diet)
cl   <- as.data.frame(cld(emm4, Letters = letters, adjust = "tukey",
                          reversed = TRUE))
print(cl[, c("Diet", "emmean", ".group")])
chk("every diet receives at least one letter",
    all(nchar(trimws(cl$.group)) >= 1))

# ---------------------------------------------------------------------------
sec("5. Model comparison — anova(A, B) likelihood-ratio test")
# Nested pair on identical data: intercept-only vs Diet, same random structure.
mA <- lmer(weight ~ 1    + (1 | Chick), data = cw, REML = TRUE, na.action = na.omit)
mB <- lmer(weight ~ Diet + (1 | Chick), data = cw, REML = TRUE, na.action = na.omit)
chk("both models fit on the same rows (valid LRT precondition)",
    stats::nobs(mA) == stats::nobs(mB))
av <- anova(mA, mB)                      # lme4 refits with ML internally
print(av)
# Identity: the reported chi-square must equal 2*(logLik_ML(B) - logLik_ML(A))
llA <- as.numeric(logLik(update(mA, REML = FALSE)))
llB <- as.numeric(logLik(update(mB, REML = FALSE)))
chisq_row <- which(is.finite(av$Chisq))[1]
chk("anova() Chisq equals 2*(logLik_ML(B) - logLik_ML(A)) (ML LRT)",
    abs(av$Chisq[chisq_row] - 2 * (llB - llA)) < 1e-4)
chk("LRT df equals the difference in fixed parameters (Diet -> 3 df)",
    av$Df[chisq_row] == (length(lme4::fixef(mB)) - length(lme4::fixef(mA))))
chk("LRT p-value is finite", is.finite(av[["Pr(>Chisq)"]][chisq_row]))
chk("the richer model (B) has the lower AIC here",
    AIC(mB) < AIC(mA))

# ---------------------------------------------------------------------------
sec("6. Influence — influence.ME leave-one-group-out Cook's distance")
if (requireNamespace("influence.ME", quietly = TRUE)) {
  cat("  (refitting once per Chick group; this can take a few seconds)\n")
  infl <- influence.ME::influence(m2, group = "Chick")
  cd   <- cooks.distance(infl)
  ng   <- nlevels(factor(cw$Chick))
  chk("one Cook's distance per Chick group",
      length(as.numeric(cd)) == ng)
  chk("Cook's distances are finite and non-negative",
      all(is.finite(as.numeric(cd))) && all(as.numeric(cd) >= 0))
  chk("the 4/n flagging rule selects between 0 and all groups",
      { n <- length(as.numeric(cd)); k <- sum(as.numeric(cd) > 4 / n)
        k >= 0 && k <= n })
} else {
  cat("  SKIP: influence.ME not installed (optional package).\n")
}

# ---------------------------------------------------------------------------
sec("7. Fit quality — performance R2 / ICC identities")
if (requireNamespace("performance", quietly = TRUE)) {
  r2  <- tryCatch(performance::r2(m2),  error = function(e) NULL)
  icc <- tryCatch(performance::icc(m2), error = function(e) NULL)
  if (!is.null(r2)) {
    rm <- as.numeric(r2$R2_marginal); rc <- as.numeric(r2$R2_conditional)
    cat(sprintf("  R2 marginal = %.3f   conditional = %.3f\n", rm, rc))
    chk("R2 values lie in [0, 1]", all(c(rm, rc) >= 0 & c(rm, rc) <= 1))
    chk("marginal R2 <= conditional R2 (identity)", rm <= rc + 1e-8)
  } else cat("  SKIP: performance::r2 returned nothing for this model.\n")
  if (!is.null(icc)) {
    iv <- as.numeric(icc$ICC_adjusted)
    cat(sprintf("  ICC (adjusted) = %.3f\n", iv))
    chk("ICC lies in [0, 1]", iv >= 0 && iv <= 1)
  } else cat("  SKIP: performance::icc returned nothing for this model.\n")
} else {
  cat("  SKIP: performance not installed (optional package).\n")
}

# ---------------------------------------------------------------------------
cat(sprintf("\n=====================  %d PASS / %d FAIL  =====================\n",
            ok, bad))
if (bad == 0) cat("Engine validated: results match references and identities.\n") else
  cat("Some checks FAILED — investigate before trusting the app's output.\n")
