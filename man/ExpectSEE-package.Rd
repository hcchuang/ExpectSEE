\name{ExpectSEE-package}
\alias{ExpectSEE-package}
\alias{expectile_loss}
\alias{expectile_weight}
\alias{fit_expectile}
\alias{fit_expectile_lasso}
\alias{generated_se}
\alias{generated_se_i}
\alias{hdsals}
\alias{hdsfgals}
\alias{sals}
\alias{sandwich_se}
\alias{see_backtest}
\alias{see_build_design}
\alias{see_compare_forecasts}
\alias{see_diagnostics}
\alias{see_fit}
\alias{see_gof}
\alias{see_aicr_simple}
\alias{see_aicr_reference_scale}
\alias{see_lambda_grid}
\alias{see_linear_test}
\alias{see_model_selection}
\alias{see_predict}
\alias{see_specs}
\alias{see_spec_test}
\alias{see_traffic_light}
\alias{see_ts_splits}
\alias{see_wald_test}
\alias{sfgals}
\title{ExpectSEE: Systems of Equations for Expectiles}
\description{
ExpectSEE provides modular tools for systems of equations for expectiles.
Use \code{see_fit()} for the combined workflow, or call \code{sals()},
\code{sfgals()}, standard-error, diagnostic, selection, testing, and
backtesting helpers separately.
}
\usage{
see_fit(y, spec = "ABSQ", tau = 0.99855, q = 1,
        estimator = c("sfgals", "sals", "hdsfgals", "hdsals"),
        se = c("generated", "sandwich", "none"),
        own_only = FALSE, hac_lag = 0,
        coefficient_fit = c("lasso", "post_lasso"),
        nested_residual_cv = NULL,
        stage2_lag_candidates = c("full", "stage1_selected"),
        stage2_penalty = c("all", "stage1_lags_unpenalized"), ...)

sals(y, spec = "ABSQ", tau = 0.99855, q = 1,
     own_only = FALSE, hac_lag = 0)

sfgals(stage1, se = c("generated", "sandwich", "none"), hac_lag = 0)

hdsals(y, spec = "ABSQ", tau = 0.99855, q = 1,
       lambda_choice = c("lambda.min", "lambda.1se"),
       foldid = NULL, nfolds = 5, standardize = TRUE,
       own_only = FALSE,
       coefficient_fit = c("lasso", "post_lasso"), hac_lag = 0,
       cv = c("iid_kfold", "fixed_window", "expanding_window"),
       ts_splits = NULL, initial = NULL, assess = NULL,
       step = NULL, lambda = NULL, lambda2 = 0, zero_tol = 1e-10, ...)

hdsfgals(stage1, lambda_choice = c("lambda.min", "lambda.1se"),
         foldid = NULL, nfolds = 5, standardize = TRUE,
         coefficient_fit = NULL, se = NULL,
         hac_lag = 0,
         cv = c("iid_kfold", "fixed_window", "expanding_window"),
         ts_splits = NULL, initial = NULL, assess = NULL,
         step = NULL, lambda = NULL, lambda2 = 0, zero_tol = 1e-10,
         nested_residual_cv = NULL,
         stage2_lag_candidates = c("full", "stage1_selected"),
         stage2_penalty = c("all", "stage1_lags_unpenalized"), ...)

see_lambda_grid(nlambda = 100, lambda_max = 1, lambda_min = 1e-10)

see_predict(fit, y, stage = NULL)

see_aicr_simple(loss_sum, n, p_effective, reference_scale)

see_aicr_reference_scale(resid, tau)

see_gof(y, yhat, p, tau, reference_scale = NULL, p_effective = p)

see_model_selection(y, specs = see_specs(), taus = 0.99855,
                    estimator = c("sfgals", "sals", "hdsfgals", "hdsals"),
                    criterion = c("AICR_simple", "loss_mean", "R2", "adjR2"),
                    q = 1, own_only = FALSE, se = "generated", ...)

see_spec_test(fit, against = c("ABN", "ABP", "AB", "SQ"), equation = NULL,
              stage = NULL, se = c("generated", "sandwich"))
}
\arguments{
\item{y}{Numeric matrix or data frame with one column per equation.}
\item{spec}{One of \code{"AB"}, \code{"SQ"}, \code{"ABP"}, \code{"ABN"}, or \code{"ABSQ"}.}
\item{tau}{Expectile level in \code{(0, 1)}.}
\item{yhat}{Fitted values corresponding to a numeric response vector \code{y}.}
\item{p}{Number of fitted coefficients, including the intercept.}
\item{q}{Number of lags.}
\item{estimator}{One of \code{"sals"}, \code{"sfgals"}, \code{"hdsals"}, or \code{"hdsfgals"}.}
\item{se}{Standard-error type used for reported t statistics. Direct Lasso
uses \code{"none"}; post-Lasso may report conventional refit-based sandwich or
generated-residual standard errors.}
\item{own_only}{If \code{TRUE}, use only own lags in each equation.}
\item{hac_lag}{Optional HAC lag for sandwich covariance estimates.}
\item{stage1}{A fitted object from \code{sals()} or \code{hdsals()}.}
\item{fit}{A fitted ExpectSEE stage object or \code{see_fit()} wrapper.}
\item{lambda_choice}{\code{"lambda.min"} for minimum CV loss or \code{"lambda.1se"} for the largest lambda within one standard error.}
\item{foldid}{Optional fold labels for \pkg{SALES} random-fold CV.}
\item{nfolds}{Number of random CV folds when \code{foldid} is not supplied.}
\item{standardize}{Passed to the \pkg{SALES} lasso backend.}
\item{coefficient_fit}{For HDSALS/HDSFGALS, \code{"lasso"} reports penalized
coefficients and is the default; \code{"post_lasso"} refits the selected active
design by unpenalized expectile IWLS. In \code{hdsfgals()}, \code{NULL} inherits
the Stage-1 mode.}
\item{cv}{Lasso cross-validation scheme: \code{"iid_kfold"}, \code{"fixed_window"}, or \code{"expanding_window"}.}
\item{ts_splits}{Optional custom time-series split list with \code{train} and \code{test} row indices.}
\item{initial}{Initial training-window length for time-series CV.}
\item{assess}{Validation-window length for time-series CV.}
\item{step}{Step size between time-series validation windows.}
\item{lambda}{Optional manual lambda grid passed to \pkg{SALES}. Leave as \code{NULL} to use the native \pkg{SALES} lambda path.}
\item{lambda2}{Optional ridge penalty passed to \pkg{SALES}.}
\item{zero_tol}{Coefficient magnitude below which lasso terms are treated as zero.}
\item{nested_residual_cv}{For \code{hdsfgals()}, whether to recompute Stage 1 selections and generated residuals inside each time-series validation fold. The default is \code{TRUE} for \code{"fixed_window"} and \code{"expanding_window"}, and \code{FALSE} for \code{"iid_kfold"}.}
\item{stage2_lag_candidates}{For \code{hdsfgals()}, \code{"full"} uses the full lag design from Stage 1 as Stage 2 lag candidates; \code{"stage1_selected"} uses only the target equation's Stage 1 selected lags.}
\item{stage2_penalty}{For \code{hdsfgals()}, \code{"all"} penalizes both Stage 2 lag and generated-residual candidates. \code{"stage1_lags_unpenalized"} keeps Stage 1 selected lags unpenalized, reproducing the earlier audit-style behavior when combined with \code{stage2_lag_candidates = "stage1_selected"}.}
\item{nlambda}{Number of lambda values requested from the native \pkg{SALES} path, or returned by \code{see_lambda_grid()} for manual-grid checks.}
\item{lambda_max, lambda_min}{Endpoints for the optional manual log-spaced grid returned by \code{see_lambda_grid()}.}
\item{stage}{Optional stage name when \code{fit} is a \code{see_fit()} wrapper.}
\item{loss_sum}{Sum of asymmetric expectile losses for one fitted equation.}
\item{n}{Number of observations used to calculate the loss.}
\item{p_effective}{Design rank including the intercept; for Lasso, the rank of the active design.}
\item{reference_scale}{Positive weighted residual scale from the common ABSQ reference fit.}
\item{specs}{Candidate SEE specifications.}
\item{taus}{One or more expectile levels; model selection returns one best row per level.}
\item{criterion}{Model-selection measure. Lower is better except for \code{R2} and \code{adjR2}.}
\item{resid}{Residual vector or matrix.}
\item{against}{Nested specification compared with ABSQ: \code{"ABN"},
\code{"ABP"}, \code{"AB"}, or \code{"SQ"}.}
\item{equation}{Optional equation name or index.}
\item{...}{Extra arguments passed to lower-level fitting functions; for lasso, this includes \code{cv}, \code{initial}, \code{assess}, \code{step}, and \code{ts_splits}.}
}
\details{
The five built-in specifications are available from \code{see_specs()}.
Older aliases are accepted: \code{"ABS"} maps to \code{"AB"},
\code{"PSQ"} to \code{"ABP"}, \code{"NSQ"} to \code{"ABN"}, and
\code{"Both"} to \code{"ABSQ"}.

The core IWLS weight is \code{abs(tau - (resid < 0))}, matching the
paper-code convention. User-facing generated-regressor corrected covariance
functions are named \code{generated_se()} and \code{generated_se_i()}.

The high-dimensional functions \code{fit_expectile_lasso()}, \code{hdsals()},
and \code{hdsfgals()} require the \pkg{SALES} backend. Their default lambda
rule is \code{"lambda.min"} from cross-validation; \code{"lambda.1se"} is also
supported. With \code{lambda = NULL}, \pkg{SALES} builds the lambda path from
the current fitting window, and ExpectSEE reuses that realized path across the
corresponding validation folds. The default \code{cv = "iid_kfold"} uses
ordinary row-wise K-fold CV for cross-sectional or iid rows. Set
\code{cv = "fixed_window"} or \code{cv = "expanding_window"} to use forward
time-series validation splits from \code{see_ts_splits()}. For
\code{cv = "iid_kfold"}, choose the number of folds with \code{nfolds}, or pass
explicit labels through \code{foldid}. For the time-ordered choices,
\code{nfolds} is ignored; use \code{initial}, \code{assess}, \code{step}, or
\code{ts_splits}. In \code{hdsfgals()}, Stage 2 defaults to a fresh
equation-by-equation LASSO-SEE using the full lag candidate matrix plus all
other-equation Stage 1 residuals; both lag and residual candidates are selected
by the Stage 2 lasso. Set \code{stage2_lag_candidates = "stage1_selected"} and
\code{stage2_penalty = "stage1_lags_unpenalized"} to reproduce the older
behavior in which Stage 1 selected lags are the only lag inputs and are kept in
the Stage 2 refit. For HDSFGALS with time-series CV,
\code{nested_residual_cv = TRUE} recomputes Stage 1 selected lags and generated
residuals inside each Stage 2 validation fold, so lambda selection does not use
generated residuals fitted on future validation rows. \code{see_lambda_grid()}
is available only when an explicit manual grid is desired for robustness or
debugging. HDSALS and HDSFGALS support \code{coefficient_fit = "lasso"}
(default direct penalized coefficients) and \code{"post_lasso"} (unpenalized
IWLS on the Lasso-selected active design). Both modes retain the same Lasso
selection and tuning audit trail. Post-Lasso standard errors are conventional
conditional-on-selected-model quantities, not selection-adjusted inference;
high-dimensional Wald and specification tests remain unavailable. For matched
\code{"iid_kfold"} mode comparisons, pass the same \code{foldid} to both fits.

The package's simple model-comparison criterion is
\deqn{AICR_{simple}=2L_{tau}/(n sigma^2_{ref})+2p/n,}
where \eqn{L_{tau}} is summed expectile loss, \eqn{p} includes the intercept,
and one ABSQ reference scale is held fixed across the declared comparison set.
Within one \code{see_model_selection()} call, each expectile/equation comparison
uses the same estimator, lag order, coefficient mode, and usable outcome rows;
the candidate-level score is the finite mean of the equation-level criteria.
This transparent loss-plus-penalty approximation is not the exact general
matrix-trace AICR. Cross-call or cross-family comparisons are valid only when
the sample and ABSQ reference-scale construction are deliberately matched.

\code{see_gof()} defines the expectile pseudo-\eqn{R^2} by comparing fitted
expectile loss with the minimized loss from an intercept-only fit at the same
\eqn{\tau}. Its adjusted version is
\eqn{1-(1-R^2)(n-1)/(n-p)}, where \eqn{p} includes the intercept. Both are
undefined when the intercept-only loss is numerically zero.

\code{see_spec_test()} tests an ABSQ fixed-dimensional fit against ABN, ABP,
AB, or SQ using all source equations and lags. Impulse responses for expectiles
(IRFE) of SFGALS and HDSFGALS systems are computed by \code{see_irfe()}.

\code{see_predict()} applies a fitted stage object to a supplied multivariate
panel. For SFGALS and HDSFGALS it first predicts the Stage 1 equations and then
uses the resulting generated residuals for the Stage 2 fitted values.
}
\value{
Fitting functions return S3 lists with coefficient tables, fitted values,
residuals, goodness-of-fit measures, diagnostics, and covariance matrices.
High-dimensional fits also include lambda-path tables; nested HDSFGALS fits
include \code{nested_stage1_path} for auditing the Stage 1 candidate lambda
paths used inside each Stage 2 validation fold. HDSFGALS fits also include
\code{stage2_candidates} and \code{stage2_support} for auditing the Stage 2
candidate set and final selected support equation by equation.
High-dimensional objects and their coefficient/GOF tables record
\code{coefficient_fit} as \code{"lasso"} or \code{"post_lasso"}.
}
\examples{
set.seed(1)
y <- matrix(rnorm(180), 60, 3)
colnames(y) <- c("A", "B", "C")
fit <- see_fit(y, spec = "AB", tau = 0.8)
coef(fit)
see_diagnostics(fit)
}
\keyword{package}
