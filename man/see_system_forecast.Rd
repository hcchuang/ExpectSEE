\name{see_system_map}
\alias{see_system_map}
\alias{see_forecast}
\title{Reduced-System Forecasts}
\usage{
see_system_map(fit, rcond_tol = 1e-8)

see_forecast(fit, history = NULL, horizon = 1L,
  calibration = c("in_sample", "provided"),
  innovation_expectile = NULL, rcond_tol = 1e-8)
}
\arguments{
\item{fit}{A SALS, SFGALS, HDSALS or HDSFGALS estimate, or a \code{see_fit} wrapper.}
\item{rcond_tol}{Minimum reciprocal condition number of \eqn{I-D}. Systems at or below this threshold are rejected.}
\item{history}{Numeric matrix, oldest row first, containing at least \code{q} observations known at the forecast origin. Named columns are aligned by equation; unnamed columns use fitted equation order. Defaults to the final fitted outcome history.}
\item{horizon}{Must be one.}
\item{calibration}{\code{in_sample} estimates a constant transformed-innovation expectile from the fitted training residuals. \code{provided} uses \code{innovation_expectile}.}
\item{innovation_expectile}{Finite vector with one correction per equation, optionally named. Applies only with \code{calibration="provided"}.}
}
\description{
Builds a common lag basis for both stages, solves the contemporaneous system,
and provides one-step forecasts. Impulse responses for expectiles are computed
by \code{\link{see_irfe}}.
}
\details{
Let \eqn{F^S} be the first-stage lag map, \eqn{F^G} the second-stage lag map,
and \eqn{D} the matrix of generated-residual coefficients (zero diagonal).
The reduced map is \eqn{g=A(F^G-DF^S)}, where \eqn{A=(I-D)^{-1}}.
Missing selected lag terms are zero on a common basis. For first-stage-only
inputs, \eqn{D=0} and \eqn{g=F^S}. No pseudoinverse or automatic shrinkage
is used to repair an ill-conditioned system. The spectral radius of D and
the inverse norm describe contemporaneous conditioning, not dynamic stability.

The forecast is \eqn{g+c_\tau}, where \eqn{c_\tau} is the componentwise
expectile of \eqn{A\eta}. It is a conditional expectile only when the
conditional transformed-innovation expectile at the origin equals that
correction. Constant calibration is an assumption, not a diagnostic conclusion.
Default calibration uses only stored training residuals and must be recomputed
with a fit trained at each chosen origin. High-dimensional consistency needs
uniform inverse, prediction and calibration control; invertibility alone is
insufficient. If SALS is correctly specified for the same information set,
both population forecasting formulas equal its conditional expectile.

Only supported lag transformations are used. No future target outcomes are
required or accepted. The caller must ensure that the fitted model and supplied
history were available at the origin. Date fields record labels and do not
certify an externally fitted model's training provenance. Unknown future
exogenous regressors and direct multi-step conditional expectiles are outside
this interface. \code{see_predict} retains its existing fitted-value semantics;
its Stage-2 residual regressors use contemporaneous outcomes.

No forecast confidence intervals are provided. Joint coefficient covariance
from \code{see_vcov_joint} does not include estimated calibration uncertainty.
Penalized and post-selection inference require additional theory.
}
\value{
\code{see_system_map} returns a list with equation names, the common \code{basis},
\code{stage1_coefficients}, \code{stage2_coefficients}, \code{D}, \code{A},
\code{reduced_coefficients}, transformed training residuals, and conditioning
diagnostics. Coefficient matrices have target equations in rows and terms in columns.
The history and estimation-end label are retained for forecasting.

\code{see_forecast} returns a data frame with one row per equation, containing
\code{forecast}, \code{reduced_component}, \code{innovation_expectile}, origin and
estimation labels, estimator information, conditioning diagnostics and an
explicitly false \code{inference_available} flag.
}
\examples{
set.seed(24)
y <- matrix(rnorm(600), 200, 3, dimnames = list(NULL, c("A", "B", "C")))
fit <- see_fit(y, spec = "AB", tau = 0.8, estimator = "sfgals")
see_system_map(fit)$diagnostics
see_forecast(fit)
}
