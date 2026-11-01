\name{see_vcov_joint}
\alias{see_vcov_joint}
\title{Joint Covariance of Fixed-Dimensional SEE Coefficients}
\usage{
see_vcov_joint(fit, hac_lag = 0)
}
\arguments{
  \item{fit}{A fixed-dimensional SALS or SFGALS fit, or a \code{see_fit}
  wrapper whose active stage is one of these estimators.}
  \item{hac_lag}{Nonnegative integer Bartlett HAC bandwidth, smaller than
  the estimation sample size. Zero uses the uncentered IID outer product.}
}
\description{
Stacks coefficient influence functions across equations and, for SFGALS,
across both estimation stages. Retains all cross-equation and cross-stage
covariance blocks.
}
\details{
SALS coefficient influences use the expectile score and inverse weighted
design Hessian. SFGALS coefficient influences use the generated-residual
corrected score from \code{generated_se_i}, multiplied by its inverse
Hessian. The joint covariance is the IID or Bartlett HAC outer product of
the stacked coefficient influences, divided by the sample size. Influences
are not centered, consistently with the package's existing sandwich and
generated-residual covariance estimators. The requested bandwidth is
recomputed from the fitted data, independently of the bandwidth or standard
error option used when fitting.

The ordering is SALS before SFGALS, then the fitted outcome-column order,
then each equation's \code{term_names} order. SALS inputs return only the
SALS coefficients. SFGALS inputs include the complete SALS first stage.
Designs must be finite, have fewer columns than observations, have full
column rank and nonsingular weighted Hessians, and agree with the stored
coefficients and residuals. Both stages must use the same observations in
the same order. Invalid objects are rejected; a pseudoinverse is not used
to justify inference. High-dimensional Lasso and post-Lasso fits are rejected.

The asymptotic interpretation remains conditional on fixed-dimensional
regularity conditions, including an appropriate law of large numbers and
central limit theorem, moment conditions, and nonsingular population
Hessians. Numerical validation does not establish these assumptions. A
fixed HAC bandwidth alone does not establish consistency under unrestricted
serial dependence. This function supplies coefficient covariance only;
\code{see_irfe(ci = TRUE)} uses it for delta-method IRFE intervals, which also
require differentiability of the nonlinear response map.
It does not include the estimated innovation-expectile calibration vector;
calibrated forecast or response inference needs its additional influence terms.
}
\value{
A list containing:
\item{vcov}{Estimated covariance matrix of the stacked coefficients.}
\item{vcov_sqrt_n}{Estimated covariance at the square-root sample-size scale.}
\item{vcov_sqrt_t}{Alias of \code{vcov_sqrt_n}, matching existing covariance helpers.}
\item{influence}{One row per aligned observation and one column per coefficient.}
\item{index}{Data frame with \code{stage}, \code{equation}, and \code{term},
in covariance-column order.}
\item{coefficients}{Named vector of stacked coefficient estimates.}
\item{n}{Estimation sample size.}
\item{hac_lag}{Bandwidth used for the joint covariance.}
Matrix and coefficient names combine stage, equation, and term;
\code{index} is the authoritative mapping.
}
\examples{
set.seed(1)
y <- matrix(rnorm(400), 200, 2, dimnames = list(NULL, c("A", "B")))
fit <- see_fit(y, spec = "AB", tau = 0.8, estimator = "sfgals")
joint <- see_vcov_joint(fit, hac_lag = 2)
head(joint$index)
sqrt(diag(joint$vcov))
}
