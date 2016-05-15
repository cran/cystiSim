\name{do_pig_vac}
\alias{do_pig_vac}

\title{Do pig vaccination}

\description{
Intervention function that mimics the possible effects of pig vaccination.
}

\usage{
do_pig_vac(x, coverage, efficacy, immunity = Inf, interval = 4,
           min.age = 1, max.age = Inf)
}

\arguments{
  \item{x}{\code{cystiRun} object.}
  \item{coverage}{Presumed coverage of vaccination within eligible population, expressed as a decimal value.}
  \item{efficacy}{Presumed efficacy of vaccine, expressed as a decimal value.}
  \item{immunity}{Presumed duration of immunity following successful vaccination. Defaults to \code{Inf}.}
  \item{interval}{Maximum interval between two consecutive shots for immunity. Defaults to 4 months.}
  \item{min.age}{Minimum age for vaccination, in months. Defaults to 1 month.}
  \item{max.age}{Maximum age for vaccination, in months. Defaults to \code{Inf}, i.e., no upper age limit.}
}

\note{
For more details and examples, please visit the \pkg{cystiSim} Wiki pages on \url{https://github.com/brechtdv/cystiSim/wiki}.
}

\author{\email{brechtdv@gmail.com}}

\seealso{Other interventions: \code{\link{do_man_mda}}, \code{\link{do_pig_mda}}, \code{\link{do_pig_mda_vac}}}