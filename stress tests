# ===========================================================================
#  Stress test: deterministic logic lifted verbatim from app.R
#  (helpers + data coercion + validation/transform guards + complete-case
#   handling). The lmer/emmeans engine is reasoned about separately.
# ===========================================================================

pass <- 0L; fail <- 0L
check <- function(label, expr) {
  res <- tryCatch(isTRUE(expr), error = function(e) structure(FALSE, msg = conditionMessage(e)))
  if (isTRUE(res)) { pass <<- pass + 1L; cat(sprintf("  ok   | %s\n", label)) }
  else { fail <<- fail + 1L
    cat(sprintf("  FAIL | %s%s\n", label,
                if (!is.null(attr(res,"msg"))) paste0("  [err: ", attr(res,"msg"), "]") else "")) }
}
section <- function(s) cat(sprintf("\n== %s ==\n", s))

# ---- helpers (copied verbatim from app.R) ---------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
bq <- function(x) sprintf("`%s`", x)

build_lhs <- function(response, transform) {
  r <- bq(response)
  switch(transform,
         "none"    = r,
         "log"     = sprintf("log(%s)", r),
         "log1p"   = sprintf("log1p(%s)", r),
         "sqrt"    = sprintf("sqrt(%s)", r),
         "inverse" = sprintf("I(1/%s)", r),
         "asin"    = sprintf("asin(sqrt(%s))", r),
         r)
}
build_random_part <- function(random, slope = NULL) {
  if (length(random) == 0) return(NULL)
  if (!is.null(slope) && nzchar(slope)) {
    paste(sprintf("(1 + %s | %s)", bq(slope), bq(random)), collapse = " + ")
  } else {
    paste(sprintf("(1 | %s)", bq(random)), collapse = " + ")
  }
}
build_formula_string <- function(response, transform, fixed, random,
                                 interactions, slope = NULL) {
  lhs <- build_lhs(response, transform)
  if (length(fixed) == 0) {
    fixed_part <- "1"
  } else if (isTRUE(interactions) && length(fixed) > 1) {
    fixed_part <- paste(bq(fixed), collapse = " * ")
  } else {
    fixed_part <- paste(bq(fixed), collapse = " + ")
  }
  random_part <- build_random_part(random, slope)
  paste(lhs, "~", fixed_part, "+", random_part)
}
tran_for_emmeans <- function(transform) {
  switch(transform, "none"=NULL, "log"="log", "log1p"="log1p", "sqrt"="sqrt",
         "inverse"="inverse", "asin"="asin.sqrt", NULL)
}
is_categorical <- function(df, vars) {
  vapply(vars, function(v) is.factor(df[[v]]) || is.character(df[[v]]), logical(1))
}
round_df <- function(d, digits = 4) {
  num  <- vapply(d, is.numeric, logical(1))
  is_p <- grepl("p\\.value|p_value|Pr\\(|^p$", names(d), ignore.case = TRUE)
  for (j in which(num))
    d[[j]] <- if (is_p[j]) signif(d[[j]], 3) else round(d[[j]], digits)
  d
}

# Reproduce the run-observer validation + transform guard logic -------------
validate_spec <- function(df, response, transform, fixed, random,
                          MAX_FIXED = 3L) {
  problems <- character(0)
  if (is.null(fixed)  || length(fixed)  == 0) problems <- c(problems, "no fixed")
  if (length(fixed) > MAX_FIXED)              problems <- c(problems, "too many fixed")
  if (is.null(random) || length(random) == 0) problems <- c(problems, "no random")
  if (length(problems)) return(list(ok = FALSE, problems = problems))
  y <- df[[response]]
  bad <- switch(transform,
    "log"     = if (any(y <= 0, na.rm = TRUE)) "log>0",
    "log1p"   = if (any(y <= -1, na.rm = TRUE)) "log1p>-1",
    "sqrt"    = if (any(y < 0,  na.rm = TRUE)) "sqrt>=0",
    "inverse" = if (any(y == 0, na.rm = TRUE)) "inv!=0",
    "asin"    = if (any(y < 0 | y > 1, na.rm = TRUE)) "asin[0,1]",
    NULL)
  if (!is.null(bad)) return(list(ok = FALSE, transform_error = bad))
  list(ok = TRUE)
}

# Reproduce dataset() coercion + type overrides -----------------------------
coerce_raw <- function(df) {
  df[] <- lapply(df, function(x) {
    if (is.character(x)) factor(x)
    else if (is.ordered(x)) factor(x, ordered = FALSE)
    else x
  })
  df
}
apply_overrides <- function(df, force_factor = NULL, force_numeric = NULL) {
  ff <- intersect(force_factor,  names(df))
  fn <- intersect(force_numeric, names(df))
  for (v in setdiff(ff, fn)) df[[v]] <- factor(df[[v]])
  for (v in fn) df[[v]] <- suppressWarnings(as.numeric(as.character(df[[v]])))
  df
}

# ===========================================================================
section("build_lhs / formula construction (now backticked)")
check("none passthrough",          build_lhs("y","none") == "`y`")
check("log wraps",                 build_lhs("y","log")  == "log(`y`)")
check("inverse uses I()",          build_lhs("y","inverse") == "I(1/`y`)")
check("asin nests",                build_lhs("y","asin") == "asin(sqrt(`y`))")
check("unknown transform falls back", build_lhs("y","bogus") == "`y`")

check("1 fixed + 1 random",
  build_formula_string("y","none","A","g",FALSE) == "`y` ~ `A` + (1 | `g`)")
check("3 fixed additive",
  build_formula_string("y","none",c("A","B","C"),"g",FALSE) == "`y` ~ `A` + `B` + `C` + (1 | `g`)")
check("interactions on >1",
  build_formula_string("y","none",c("A","B"),"g",TRUE) == "`y` ~ `A` * `B` + (1 | `g`)")
check("interactions ignored for single fixed",
  build_formula_string("y","none","A","g",TRUE) == "`y` ~ `A` + (1 | `g`)")
check("multiple random intercepts",
  build_formula_string("y","none","A",c("g1","g2"),FALSE) == "`y` ~ `A` + (1 | `g1`) + (1 | `g2`)")
check("random slope",
  build_formula_string("y","log","A","g",FALSE,slope="A") == "log(`y`) ~ `A` + (1 + `A` | `g`)")
check("empty slope string == intercept only",
  build_formula_string("y","none","A","g",FALSE,slope="") == "`y` ~ `A` + (1 | `g`)")
check("empty fixed -> intercept-only fixed part",
  build_formula_string("y","none",character(0),"g",FALSE) == "`y` ~ 1 + (1 | `g`)")
check("constructed formula parses",
  inherits(tryCatch(as.formula(build_formula_string("y","sqrt",c("A","B"),"g",TRUE)),
                    error=function(e) e), "formula"))

section("variable names with spaces / specials (now FIXED by backticking)")
df_sp <- data.frame(check.names = FALSE,
                    `resp val` = c(1,2,3,4),
                    `Trt Grp`  = factor(c("a","b","a","b")),
                    subj       = factor(c("s1","s1","s2","s2")))
fstr <- build_formula_string("resp val","none","Trt Grp","subj",FALSE)
cat("    produced:", fstr, "\n")
check("space-name formula NOW parses",
      inherits(tryCatch(as.formula(fstr), error=function(e) e), "formula"))
# and the model frame can actually be built from it
mf <- tryCatch(model.frame(as.formula("`resp val` ~ `Trt Grp`"), data = df_sp),
               error = function(e) e)
check("model.frame builds from backticked space-names", is.data.frame(mf))
check("leading-digit name parses too",
      inherits(tryCatch(as.formula(build_formula_string("1y","none","2x","g",FALSE)),
                        error=function(e) e), "formula"))

section("tran_for_emmeans mapping (back-transform correctness)")
check("none -> NULL",      is.null(tran_for_emmeans("none")))
check("log -> log",        identical(tran_for_emmeans("log"), "log"))
check("inverse -> inverse",identical(tran_for_emmeans("inverse"), "inverse"))
check("asin -> asin.sqrt", identical(tran_for_emmeans("asin"), "asin.sqrt"))
check("log1p -> log1p",    identical(tran_for_emmeans("log1p"), "log1p"))
# sanity: the named trans are the ones emmeans accepts as strings / links
check("inverse is a valid stats link",
      inherits(tryCatch(stats::make.link("inverse"), error=function(e) e), "link-glm"))

section("validation guards (run observer)")
d <- data.frame(y = c(1,2,3,4), A = factor(c("a","b","a","b")), g = factor(c(1,1,2,2)))
check("no fixed -> blocked",      isFALSE(validate_spec(d,"y","none",NULL,"g")$ok))
check("no random -> blocked",     isFALSE(validate_spec(d,"y","none","A",NULL)$ok))
check("4 fixed -> blocked",       isFALSE(validate_spec(d,"y","none",c("A","B","C","D"),"g")$ok))
check("valid spec -> ok",         isTRUE(validate_spec(d,"y","none","A","g")$ok))

section("transform guards on nasty responses")
dz  <- data.frame(y = c(0, 1, 2, 3),  A=factor(c("a","b","a","b")), g=factor(c(1,1,2,2)))
dng <- data.frame(y = c(-1, 1, 2, 3), A=factor(c("a","b","a","b")), g=factor(c(1,1,2,2)))
dp  <- data.frame(y = c(0.1,0.5,0.9,1.5), A=factor(c("a","b","a","b")), g=factor(c(1,1,2,2)))
dna <- data.frame(y = c(NA, 1, 2, 3), A=factor(c("a","b","a","b")), g=factor(c(1,1,2,2)))
check("log blocks zero",          !is.null(validate_spec(dz,"y","log","A","g")$transform_error))
check("log1p blocks < -1",        !is.null(validate_spec(dng,"y","log1p","A","g")$transform_error))
check("sqrt blocks negative",     !is.null(validate_spec(dng,"y","sqrt","A","g")$transform_error))
check("inverse blocks zero",      !is.null(validate_spec(dz,"y","inverse","A","g")$transform_error))
check("asin blocks > 1",          !is.null(validate_spec(dp,"y","asin","A","g")$transform_error))
check("asin ok on [0,1]",         isTRUE(validate_spec(data.frame(y=c(0,.3,.7,1),A=factor(c("a","b","a","b")),g=factor(c(1,1,2,2))),"y","asin","A","g")$ok))
check("NA in response: log guard ignores NA, still blocks 0? (no 0 here)",
      isTRUE(validate_spec(dna,"y","log","A","g")$ok))   # NA dropped by na.rm; values 1,2,3 > 0

section("type coercion + overrides")
raw <- data.frame(id = c(1,2,3), grp = c("x","y","x"), val = c(1.5,2.5,3.5),
                  stringsAsFactors = FALSE)
cr  <- coerce_raw(raw)
check("char -> factor",            is.factor(cr$grp))
check("numeric stays numeric",     is.numeric(cr$val))
ov  <- apply_overrides(cr, force_factor = "id")
check("force id -> factor",        is.factor(ov$id))
ov2 <- apply_overrides(cr, force_numeric = "grp")  # "x"/"y" -> NA with warning
check("force non-numeric -> all NA (documents data-loss risk)",
      all(is.na(ov2$grp)))
ov3 <- apply_overrides(cr, force_factor = "val", force_numeric = "val")
check("same var in both lists: numeric wins (setdiff guard)",
      is.numeric(ov3$val))
ord <- data.frame(o = ordered(c("lo","hi","lo"), levels=c("lo","hi")))
check("ordered factor demoted to plain factor",
      is.factor(coerce_raw(ord)$o) && !is.ordered(coerce_raw(ord)$o))

section("interaction-plot complete-case logic")
di <- data.frame(y = c(1,2,NA,4,5,6),
                 A = factor(c("a","a","b","b",NA,"a")),
                 B = factor(c("x","y","x","y","x",NA)))
cols <- c("y","A","B"); cc <- complete.cases(di[cols])
check("complete.cases drops rows with any NA in used cols", sum(cc) == 3)
ag <- aggregate(di[cc,"y"], by = list(A=di[cc,"A"], B=di[cc,"B"]),
                FUN = function(z) mean(z, na.rm=TRUE))
check("aggregate yields finite means", all(is.finite(ag$x)))

section("round_df doesn't choke on factor/character columns")
mixed <- data.frame(f = factor(c("a","b")), n = c(1.23456, 2.34567),
                    c = c("p","q"), stringsAsFactors = FALSE)
rd <- round_df(mixed)
check("numeric rounded",           rd$n[1] == 1.2346)
check("factor untouched",          identical(rd$f, mixed$f))
check("character untouched",       identical(rd$c, mixed$c))

section("%||% null-coalesce")
check("NULL -> default",  (`%||%`)(NULL, "x") == "x")
check("value kept",       (`%||%`)("y", "x") == "y")

# ---- single-level / degenerate factors (new guard) ------------------------
section("degenerate factor guard (new)")
n_levels_used <- function(v) nlevels(droplevels(factor(v[!is.na(v)])))
degenerate_check <- function(df, fixed, random) {
  is_cat <- function(v) is.factor(df[[v]]) || is.character(df[[v]])
  deg <- character(0)
  for (v in fixed[vapply(fixed, is_cat, logical(1))])
    if (n_levels_used(df[[v]]) < 2) deg <- c(deg, paste0("fixed:", v))
  for (g in random)
    if (n_levels_used(factor(df[[g]])) < 2) deg <- c(deg, paste0("random:", g))
  deg
}
dd <- data.frame(y = rnorm(8),
                 A1 = factor(rep("only", 8)),                 # 1 level
                 A2 = factor(rep(c("a","b"), 4)),             # 2 levels
                 g1 = factor(rep("grp", 8)),                  # 1 group
                 g2 = factor(rep(c("p","q"), each = 4)),      # 2 groups
                 cont = rnorm(8))
check("single-level fixed flagged",     "fixed:A1"  %in% degenerate_check(dd,"A1","g2"))
check("two-level fixed not flagged",    length(degenerate_check(dd,"A2","g2")) == 0)
check("single-group random flagged",    "random:g1" %in% degenerate_check(dd,"A2","g1"))
check("continuous fixed never flagged", length(degenerate_check(dd,"cont","g2")) == 0)
check("NA-only-second-level collapses to 1",
      n_levels_used(factor(c("x","x",NA,NA))) == 1)

section("Type II/III contrasts, ddf coupling, balance (new)")
make_contr <- function(atype, cat_fixed) {
  if (atype == "3" && length(cat_fixed))
    stats::setNames(as.list(rep("contr.sum", length(cat_fixed))), cat_fixed)
  else NULL
}
ca <- make_contr("3", c("A","B"))
check("Type III builds contr.sum list",
      identical(ca, list(A = "contr.sum", B = "contr.sum")))
check("Type II uses default contrasts (NULL)", is.null(make_contr("2", c("A","B"))))
check("Type III with no categorical fixed -> NULL", is.null(make_contr("3", character(0))))

eff_ddf <- function(use_reml, ddf) if (!use_reml && ddf == "Kenward-Roger") "Satterthwaite" else ddf
check("ML + KR auto-switches to Satterthwaite",
      eff_ddf(FALSE, "Kenward-Roger") == "Satterthwaite")
check("REML + KR keeps Kenward-Roger",
      eff_ddf(TRUE, "Kenward-Roger") == "Kenward-Roger")
check("ML + Satterthwaite unchanged",
      eff_ddf(FALSE, "Satterthwaite") == "Satterthwaite")

roman <- function(a) if (a == "3") "III" else "II"
check("type 3 -> III", roman("3") == "III")
check("type 2 -> II",  roman("2") == "II")

db <- data.frame(A = factor(c("a","a","b"), levels = c("a","b")),
                 B = factor(c("x","y","x"), levels = c("x","y")))
bt <- as.data.frame(table(db[c("A","B")]), responseName = "n")
check("balance table enumerates all level combinations", nrow(bt) == 4)
check("empty cell (b,y) shows n = 0",
      bt$n[bt$A == "b" & bt$B == "y"] == 0)
check("balance table flags >=1 empty cell", any(bt$n == 0))

section("REML-safe deviance + convergence-warning classifier (new)")
# mirror the fit_table deviance/REMLcrit name selection
dev_name <- function(is_reml) if (is_reml) "REML_criterion" else "deviance"
check("REML fit -> REML_criterion label", dev_name(TRUE) == "REML_criterion")
check("ML fit -> deviance label",          dev_name(FALSE) == "deviance")

# mirror the 'trouble' detector for convergence/identifiability warnings
trouble <- function(warn_all, singular) {
  singular || any(grepl("Hessian|gradient|converge|singular|isSingular|not uniquely",
                        warn_all, ignore.case = TRUE))
}
check("scaled-gradient warning flagged",
      trouble("unable to evaluate scaled gradient", FALSE))
check("Hessian-singular warning flagged",
      trouble("Hessian is numerically singular: parameters are not uniquely determined", FALSE))
check("isSingular alone flags trouble",  trouble(character(0), TRUE))
check("benign warning not flagged",
      !trouble("some unrelated message about levels", FALSE))
# the deviance-deprecation warning must NOT be misread as convergence trouble
check("deviance-deprecation is not convergence trouble",
      !trouble("deviance() is deprecated for REML fits; use REMLcrit", FALSE))

section("ANOVA NaN p-value / df detector (new)")
detect_bad <- function(a, warns) {
  pcol    <- intersect(c("Pr(>F)", "p.value"), names(a))
  bad_p   <- length(pcol) && any(!is.finite(a[[pcol[1]]]))
  bad_ddf <- "DenDF" %in% names(a) && any(!is.finite(a[["DenDF"]]))
  bad_p || bad_ddf || any(grepl("NaN", warns, ignore.case = TRUE))
}
good <- data.frame(check.names = FALSE, DenDF = c(10, 12),
                   "Pr(>F)" = c(0.01, 0.2))
nanp <- data.frame(check.names = FALSE, DenDF = c(10, 12),
                   "Pr(>F)" = c(0.01, NaN))
nanddf <- data.frame(check.names = FALSE, DenDF = c(NaN, 12),
                     "Pr(>F)" = c(0.01, 0.2))
check("clean table not flagged",          !detect_bad(good, character(0)))
check("NaN p-value flagged",               detect_bad(nanp, character(0)))
check("non-finite DenDF flagged",          detect_bad(nanddf, character(0)))
check("pf NaN warning flagged even if table clean",
      detect_bad(good, "NaNs produced"))

section("contrast gating + complete-case degeneracy + comparison logic (new)")
# contr.sum should apply only for Type III WITH interactions actually present
apply_sum <- function(atype, interactions, n_fixed, n_cat) {
  has_inter <- isTRUE(interactions) && n_fixed > 1
  atype == "3" && has_inter && n_cat > 0
}
check("III + interactions + categorical -> apply",  apply_sum("3", TRUE,  2, 2))
check("III + NO interactions -> do not apply",      !apply_sum("3", FALSE, 2, 2))
check("III + interactions but 1 fixed -> do not apply", !apply_sum("3", TRUE, 1, 1))
check("II + interactions -> do not apply",          !apply_sum("2", TRUE,  2, 2))
check("III + interactions but no categorical -> do not apply",
      !apply_sum("3", TRUE, 2, 0))

# complete-case degeneracy: factor collapsing to one level after NA drop
md <- data.frame(y = c(1,2,3,NA), A = factor(c("a","a","a","b")),
                 g = factor(c(1,2,1,2)))
mv <- c("y","A","g"); cc <- stats::complete.cases(md[, mv])
nlu <- function(v) nlevels(droplevels(factor(v)))
check("A has 2 levels in full data", nlevels(md$A) == 2)
check("A collapses to 1 level on complete cases (now caught)",
      nlu(md[cc, "A"]) == 1)
check("g still has 2 levels on complete cases", nlu(md[cc, "g"]) >= 2)

# comparison classifiers
a <- list(fixed=c("A"),     random="g", slope="", fml_str="y ~ A + (1|g)")
b1 <- list(fixed=c("A","B"),random="g", slope="", fml_str="y ~ A + B + (1|g)")
b2 <- list(fixed=c("A"),    random="g", slope="A",fml_str="y ~ A + (1+A|g)")
fixed_differ  <- function(a,b) !identical(sort(a$fixed), sort(b$fixed)) || !identical(a$fml_str, b$fml_str)
random_differ <- function(a,b) !identical(sort(a$random), sort(b$random)) || !identical(a$slope %||% "", b$slope %||% "")
check("differing fixed effects detected",  fixed_differ(a, b1))
check("same random when only fixed differ", !random_differ(a, b1))
check("differing random slope detected",    random_differ(a, b2))

section("display accuracy + few-level random-effect guard (new)")
# p-values must not collapse to 0 in the displayed tables
pf <- data.frame(contrast = c("a-b","a-c"),
                 estimate = c(1.23456, 2.34567),
                 p.value  = c(0.00000234, 0.5))
rp <- round_df(pf)
check("estimate rounded to 4 dp", rp$estimate[1] == 1.2346)
check("tiny p-value kept via signif (not 0)", rp$p.value[1] == 2.34e-06)
check("tiny p-value is non-zero", rp$p.value[1] > 0)
check("ordinary p-value preserved", rp$p.value[2] == 0.5)
# a column literally named with Pr( ) is also treated as a p-value
pf2 <- data.frame(check.names = FALSE, "Pr(>F)" = c(1e-9, 0.2))
check("Pr(>F) column kept via signif",
      round_df(pf2)[["Pr(>F)"]][1] == 1e-09)

# few-level random-effect detector mirrors the run observer
flag_few <- function(ngrps) names(ngrps[ngrps < 5])
check("4-level random effect flagged (e.g. ChickWeight Diet)",
      identical(flag_few(c(Diet = 4L, Chick = 50L)), "Diet"))
check(">=5-level random effect not flagged",
      length(flag_few(c(Subject = 18L))) == 0)
check("two sparse factors both flagged",
      setequal(flag_few(c(A = 2L, B = 3L, C = 9L)), c("A","B")))

section("self-prediction & role-overlap guards (new)")
spec_problems <- function(response, fixed, random, MAX_FIXED = 3L) {
  p <- character(0)
  if (is.null(fixed)  || length(fixed)  == 0) p <- c(p, "no fixed")
  if (length(fixed) > MAX_FIXED)              p <- c(p, "too many")
  if (is.null(random) || length(random) == 0) p <- c(p, "no random")
  if (!is.null(response) && response %in% fixed)  p <- c(p, "resp in fixed")
  if (!is.null(response) && response %in% random) p <- c(p, "resp in random")
  if (length(intersect(fixed, random)))           p <- c(p, "overlap")
  p
}
check("response as fixed effect blocked",
      "resp in fixed" %in% spec_problems("y", c("y","A"), "g"))
check("response as random effect blocked",
      "resp in random" %in% spec_problems("y", "A", c("y","g")))
check("fixed/random overlap blocked",
      "overlap" %in% spec_problems("y", c("A","g"), "g"))
check("clean spec has no problems",
      length(spec_problems("y", c("A","B"), "g")) == 0)

cat(sprintf("\n=================  %d passed, %d failed  =================\n", pass, fail))
