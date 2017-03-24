% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initialize.R
\name{initialize_pop}
\alias{initialize_pop}
\title{Initialize a FEMALE population with an age structure, dates of cancer 
incidence, and dates of all-cause mortality}
\usage{
initialize_pop(pop_size, nsim, agesource, minage, maxage, incsource, mortsource,
  pop_chars = list(male = data.frame(male = c(0), prop = c(1))))
}
\arguments{
\item{agesource}{Country to use for age structure (see data(agestructure) )}

\item{minage}{Lower age limit for population at sim start}

\item{maxage}{Upper age limit for population at sim start}

\item{incsource}{Country to use for incidence rates (see data(incratesf) )}

\item{mortsource}{Country to use for life table (see data(allmortratesf) )}

\item{pop_chars}{A list of data frames that specify additional features
to simulate in the population. Defaults to giving the whole population
male=0, i.e. all female sex.}
}
\description{
Initialize a FEMALE population with an age structure, dates of cancer 
incidence, and dates of all-cause mortality
}
\examples{
pop <- initialize_pop(pop_size=100000,
                     nsim=2, 
                     agesource='Standard', 
                     minage=0, maxage=100, 
                     incsource='Uganda', 
                     mortsource='Uganda')

}
