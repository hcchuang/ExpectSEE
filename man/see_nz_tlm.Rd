\name{see_nz_tlm}
\alias{see_nz_score}
\alias{see_nz_tlm_test}
\alias{see_parzen_lrv}
\title{Nolde-Ziegel expectile traffic-light backtest}
\description{
Expectile comparative backtesting helpers following the scoring functions
and traffic-light construction used by Nolde and Ziegel (2017).
}
\usage{
see_nz_score(y, yhat, tau, h = 1, positive_guard = FALSE)

see_nz_tlm_test(internal_yhat, standard_yhat, y, tau, alpha = 0.05,
                h = 1, positive_guard = FALSE)

see_parzen_lrv(x)
}
\arguments{
\item{y}{Observed outcomes.}
\item{yhat}{Expectile forecasts.}
\item{tau}{Expectile level in \code{(0, 1)}.}
\item{h}{Scoring-function choice. Use \code{1} for the 2-homogeneous
expectile score in Nolde and Ziegel Eq. (2.21), or \code{0} for Eq. (2.22).}
\item{positive_guard}{If \code{TRUE}, reproduce the original
\code{sfEXP()} implementation, which scores only positive forecasts and
sets excluded observations to zero. If \code{FALSE}, apply Eq. (2.21)
directly to real-valued forecasts.}
\item{internal_yhat}{Forecasts from the internal model.}
\item{standard_yhat}{Forecasts from the standard benchmark model.}
\item{alpha}{Traffic-light test level.}
\item{x}{Score-difference series.}
}
\value{
\code{see_nz_score()} returns a score vector. \code{see_parzen_lrv()}
returns the Parzen HAC long-run variance used in the original code.
\code{see_nz_tlm_test()} returns a one-row data frame with mean score
difference, test statistic, one-sided p-values, and traffic-light decision.
}
\details{
The test compares the score difference
\code{S(internal, y) - S(standard, y)}. A negative mean difference favors
the internal model. The green, yellow, and red labels follow the original
traffic-light convention: green if the internal model significantly improves
on the standard model, red if it is significantly worse, and yellow otherwise.
}
\references{
Nolde, N. and Ziegel, J. F. (2017). Elicitability and backtesting:
Perspectives for banking regulation. \emph{Annals of Applied Statistics},
11(4), 1833--1874. arXiv:1608.05498.
}
