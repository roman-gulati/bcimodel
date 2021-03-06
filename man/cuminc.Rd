% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outcomes.R
\name{cuminc}
\alias{cuminc}
\title{Tally cumulative incidence of cancer}
\usage{
cuminc(futimes, times)
}
\arguments{
\item{futimes}{Vector of follow-up times to assess}

\item{times}{Matrix of times to cancer incidence}
}
\value{
Matrix where rows tally incident cases for each sim; columns represent the various follow-up times

Matrix of incidence counts for sims (rows) and follow-up times (columns)
}
\description{
Count incident cancer cases within a follow-up period
}
\examples{
# Assess 5- and 10-year follow-ups
fu <- c(5, 10)
times <- matrix(rexp(n=25, rate=1/5), nrow=5, ncol=5)
# Rows are simulations, columns are 5- and 10-year follow-up times, 
# and cells are incidence counts within the periods
cuminc(fu, times)

}

