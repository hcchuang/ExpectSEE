\name{see_irfe}
\alias{see_irfe}
\title{Impulse Responses for Expectiles (IRFE)}
\usage{
see_irfe(fit, source = NULL, shock = 0.01, horizon = 12L,
  sign = c("both", "positive", "negative"), ci = FALSE, level = 0.95,
  hac_lag = 4L, rcond_tol = 1e-8, divergence = c("error", "record"),
  fd_step = 1e-6, center = c("zero", "fixed_point"), fixed_point = NULL,
  fp_tol = 1e-12, ci_gradient = c("recursive", "finite_difference"))
}
\arguments{
\item{fit}{An SFGALS or HDSFGALS estimate, or a \code{see_fit} wrapper. Any number of equations and lags.}
\item{source}{Names or indices of the equations whose Stage-2 innovation is shocked. \code{NULL} uses every equation.}
\item{shock}{Positive size \eqn{d} of the innovation, in the units of the outcomes (0.01 is one percentage point for decimal returns).}
\item{horizon}{Number of horizons; \eqn{h=1} is the impact period.}
\item{sign}{Positive innovations \eqn{+d}, negative innovations \eqn{-d}, or both.}
\item{ci}{If \code{TRUE}, delta-method standard errors and intervals. Fixed-dimensional SFGALS fits only.}
\item{level}{Confidence level of the intervals.}
\item{hac_lag}{Bartlett lags of the joint covariance \code{see_vcov_joint()}; 0 gives the IID version.}
\item{rcond_tol}{Minimum reciprocal condition number of \eqn{I-D}.}
\item{divergence}{\code{error} stops at a nonfinite path; \code{record} keeps the finite prefix.}
\item{fd_step}{Relative step used when \code{ci_gradient = "finite_difference"}.}
\item{ci_gradient}{\code{"recursive"} propagates the analytic Jacobian through both paths;
\code{"finite_difference"} retains complete-path coefficient perturbations as a check.}
\item{center}{Initial history of both paths: \code{"zero"} starts every lag at zero
and evolves the baseline under the original transition; \code{"fixed_point"} starts
at a verified attracting fixed point \eqn{\bar y=G(\bar y)}.}
\item{fixed_point}{Optional fixed point \eqn{\bar y} (one value per equation) for
\code{center = "fixed_point"}; it is verified. If \code{NULL}, it is found by Newton's method.}
\item{fp_tol}{Tolerance for \eqn{\max_i|G(\bar y)_i-\bar y_i|}.}
}
\value{
A data frame of class \code{expectsee_irfe} with one row per source, sign, horizon and target:
\code{irfe}, \code{baseline}, \code{shocked}, and \code{initial_transition} (the target's
\eqn{G(0)}). Deprecated \code{removed_intercept} is always \code{NA}. With \code{ci = TRUE},
\code{se}, \code{lower}, \code{upper}, \code{p_value}, \code{ci_method}, \code{ci_status}
and \code{fd_step} describe inference. \code{fd_step} is \code{NA} for recursive gradients.
\code{inference_available} is true only for usable intervals.
Conditioning diagnostics for \eqn{I-D}, \code{path_status}, \code{first_nonfinite_horizon},
\code{divergent_path} and \code{path_boundary} are included. \code{fixed_point_residual} and
\code{fixed_point_spectral_radius} are \code{NA} for zero mode.
}
\description{
Computes model-implied response paths with zero future innovations under the reduced
transition \eqn{G=A(F^G-DF^S)}, where \eqn{A=(I-D)^{-1}}. A one-time Stage-2 innovation
\eqn{d_j^{\pm}=\pm d\,\iota_j} gives impact \eqn{\mathrm{IRFE}_{j,1}^{\pm}=A d_j^{\pm}}.

With \code{center = "zero"}, both q-lag histories start at zero. At impact,
\eqn{b_1=G(0)} and \eqn{x_1=G(0)+A d_j^{\pm}}. Thereafter both histories iterate the
original \eqn{G}, and the response is \eqn{x_h-b_h}. For one lag,
\eqn{\mathrm{IRFE}_{j,2}^{\pm}=G\{G(0)+A d_j^{\pm}\}-G\{G(0)\}}.
The intercept is retained. With \code{center = "fixed_point"}, both histories start
at \eqn{\bar y=G(\bar y)}, and the baseline stays there up to solver precision.
}
\details{
Version 0.2.9 replaced the former zero-normalized method with zero-initial-state
paired paths. Version 0.2.10 uses analytic paired-path Jacobians by default for
fixed-dimensional intervals. The methods agree at impact, but generally differ at
later horizons for nonlinear maps. These responses are not direct multi-step
conditional expectiles or identified causal effects.

The fixed point is found by damped Newton iteration on \eqn{G(y)-y=0} started from
\eqn{G(0)} and from the limit of iteration from zero when available, or verified if
supplied. It must be attracting (companion spectral radius below one) and the only
attracting root among those found. This limited search does not prove global uniqueness.
Zero mode requires no fixed point.

Intervals use the joint SALS-SFGALS coefficient covariance from \code{see_vcov_joint()}:
\eqn{\mathrm{se}=\sqrt{J\hat V J'}}. That covariance already includes sample-size scaling.
At impact the Jacobian is analytic, \eqn{dA=A\,dD\,A}. At later horizons the default
recursion propagates derivatives of both complete paths, including the evolving
baseline and the SALS-SFGALS generated-regressor covariance. For fixed-point starts,
the initial-state derivative is obtained from the implicit fixed-point equation.
These calculations apply to AB, SQ, ABP, ABN, and ABSQ. Optional central differences
recompute both complete paths and re-solve \eqn{\bar y} when relevant. Shock magnitude
is held fixed; uncertainty in an estimated shock calibration is excluded.

\code{path_boundary} conservatively flags active nonsquared positive/negative inputs
within \eqn{10^{-8}} of zero at parameter-dependent states and all subsequent horizons.
Fixed initial zero rows are excluded. Boundary intervals are withheld. \code{ci_status}
also identifies nonfinite paths or derivatives and degenerate variance. Finite-difference
comparison with finite-difference gradients and branch crossings under coefficient
perturbations remain useful diagnostics. Ordinary intervals require the joint CLT,
local path smoothness and nonzero variance. High-dimensional LASSO and post-LASSO
fits have no ordinary delta-method intervals.
}
\examples{
\dontrun{
irf <- see_irfe(fit, source = "JPM", shock = 0.01, horizon = 5)
irf_fp <- see_irfe(fit, source = "JPM", shock = 0.01, horizon = 5, center = "fixed_point")
irf5 <- see_irfe(fit5, ci = TRUE, hac_lag = 4)   # fixed-dimensional SFGALS
}
}
