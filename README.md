# Quick Mixed-Model Review

An R Shiny app for fast, interactive linear mixed-model analysis built on
[`lmerTest`](https://cran.r-project.org/package=lmerTest) and
[`emmeans`](https://cran.r-project.org/package=emmeans). Pick a response, up to
three fixed effects, random grouping factor(s), and an optional transformation
from drop-down menus, and the app fits the model and returns an ANOVA, residual
and random-effects diagnostics, estimated marginal means with post-hoc
comparisons and compact-letter groupings, fit statistics, influence
diagnostics, and an optional model comparison — each with copy-and-run R code.

Built-in validation guards stop common mistakes before they reach `lmer`, and
plain-language warnings explain singular fits, convergence trouble, and other
issues instead of leaking raw console errors.

---

## Quick start

```r
# 1. Install dependencies (once)
install.packages(c(
  "shiny","shinyjs","lmerTest","emmeans","pbkrtest","multcomp",
  "multcompView","ggplot2","DT","lattice",
  "performance","influence.ME"   # optional: R2/ICC and Cook's distance
))

# 2. Run the app
shiny::runApp("app.R")
```

`performance` (marginal/conditional R² and ICC) and `influence.ME`
(leave-one-group-out Cook's distance) are optional — the app degrades
gracefully and tells you if they are missing. `pbkrtest` is required for
Kenward-Roger denominator degrees of freedom.

No data of your own? The app ships with the built-in `ChickWeight` example so
you can explore it immediately.

---

## Features

- **Drop-down model builder** — response, up to three fixed effects, random
  grouping factor(s), optional random slope, REML/ML, and a reset button.
- **Transformations** — log, log1p, sqrt, inverse (1/y), and arcsine-sqrt, each
  with validity checks; optional, correct back-transformation of EMMeans.
- **Type II / III ANOVA** with Kenward-Roger or Satterthwaite df (auto-switching
  to Satterthwaite under ML, where Kenward-Roger is not defined).
- **Diagnostics** — four-panel residual plots with plain-language guidance, a
  random-effects caterpillar plot, a BLUP normal Q-Q, and optional group-level
  Cook's distance.
- **EMMeans & post-hoc** — estimated marginal means, pairwise comparisons with
  selectable adjustment (Tukey/Šidák/Bonferroni/Holm/none), compact letter
  display, and a letter-annotated means plot (1–3 factors).
- **Fit quality** — AIC/BIC/logLik, REML criterion or deviance, effective N vs
  total, marginal/conditional R², ICC, and variance components.
- **Model comparison** — likelihood-ratio test of a saved model against the
  current one, with validity caveats.
- **Reproducible code** — every analysis/plot tab has a "copy & run in R" panel
  that emits the exact code for that output.
- **Safety rails** — guards against empty/oversized specs, self-prediction,
  fixed/random role overlap, single-level factors, and missing-data pitfalls,
  plus clear warnings for singular fits and unreliable few-level random effects.

---

## Tabs

| Tab | What it shows |
|-----|---------------|
| **Data** | Dimensions, structure, per-column missing counts, preview |
| **Explore (raw data)** | Design-balance crosstab (empty cells flagged), response-vs-effect plot, raw-means interaction plot |
| **Model & code** | Fitted formula, copy-and-run R code, `summary(mod)` |
| **Fit & variance** | AIC/BIC/logLik, REML criterion/deviance, effective N, R²/ICC, variance components, caterpillar plot, BLUP Q-Q |
| **ANOVA** | Type II/III table with the chosen denominator df |
| **Residuals** | Four-panel diagnostics + optional Cook's distance |
| **EMMeans & post-hoc** | Estimated marginal means, pairwise comparisons, compact letter display (all downloadable) |
| **EMMeans plot** | Letter-annotated means plot |
| **Compare models** | Likelihood-ratio test of saved vs current model |

---

## Repository contents

| File | Purpose |
|------|---------|
| `app.R` | The Shiny application |
| `validate_engine.R` | One-command check that the statistical engine reproduces known results (run on a machine with the packages) |
| `stress_tests.R` | Base-R unit tests for the app's deterministic logic (no modelling packages needed) |
| `README.md` | This file |

---

## Validation

The app's logic is covered two ways. The deterministic logic (formula building,
validation guards, transform checks, NA handling, display formatting) is tested
with base R:

```r
Rscript stress_tests.R        # expect: all assertions pass
```

The statistical engine itself is checked against published reference values and
identities that must hold if it is correct — the `lmer` fit and fixed effects,
the ANOVA, EMMeans and the `type="response"` back-transform (must equal
`exp(link)` for a log fit), the compact letter display, the `anova(A, B)`
likelihood-ratio test (its chi-square must equal
`2 * (logLik_ML(B) - logLik_ML(A))` with the correct df), `influence.ME` Cook's
distance, and the `performance` R²/ICC identities:

```r
Rscript validate_engine.R     # expect: all checks PASS
```

The Cook's-distance and R²/ICC checks are skipped (not failed) if those optional
packages are absent. Re-run `validate_engine.R` after any major upgrade of
`lmerTest` or `emmeans`.

---

## Statistical notes (please read)

- **Roles matter.** The response cannot also be a predictor, and a variable
  should be used in only one role. In particular, do not enter a factor as a
  fixed effect when a *coarser* grouping of it is your random effect — they are
  nested and confounded. A classic trap on `ChickWeight` is
  `weight ~ Chick + (1 | Diet)`, which is degenerate because each chick belongs
  to one diet; the sensible model is `weight ~ Diet * Time + (1 | Chick)`.
- **Random effects need enough levels.** A grouping factor with fewer than ~5
  levels gives an unreliable variance estimate; the app warns, and such a factor
  usually belongs in the model as a *fixed* effect instead.
- **Type III and contrasts.** Type III main-effect tests are contrast-dependent
  when interactions are present, so the app applies sum-to-zero contrasts in
  that case (which also re-codes the `summary()` coefficients as deviations from
  the grand mean). With no interactions, Type III equals Type II.
- **Kenward-Roger requires REML.** Under ML, the app auto-switches to
  Satterthwaite and says so.
- **Back-transformed means are not arithmetic means.** For a log fit the
  back-transformed EMMean is a geometric mean and differences appear as ratios;
  comparisons are always tested on the model (link) scale.
- **Singular / NaN-p-value warnings are diagnoses, not bugs.** They indicate an
  over-parameterized model; simplify the random structure, drop a sparse factor,
  or center/scale continuous predictors.

---

## Known limitations

- Random-effects structure is limited to random intercepts plus one optional
  random slope; nested/crossed specifications and uncorrelated slopes are not
  exposed.
- No alternative-optimizer retry on non-convergence, observation weights/offsets,
  or one-click report export — intentional omissions for a *quick review* tool.

---

## Acknowledgements

Powered by `lme4`/`lmerTest` for fitting and inference, `emmeans` and
`multcomp` for estimated marginal means and letter groupings, and `performance`
and `influence.ME` for fit quality and influence diagnostics.
