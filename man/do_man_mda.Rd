\name{do_man_mda}
\alias{do_man_mda}

\title{Do human mass drug administration(MDA)}

\description{
Intervention function that mimics the possible effects of human mass drug administration.
}

\usage{
do_man_mda(x, coverage, efficacy, min.age = 0, max.age = Inf)
}

\arguments{
  \item{x}{\code{cystiRun} object.}
  \item{coverage}{Presumed coverage of drug administration within eligible population, expressed as a decimal value.}
  \item{efficacy}{Presumed efficacy of tapeworm treatment, expressed as a decimal value.}
  \item{min.age}{Minimum age for MDA, in months. Defaults to \code{0}, i.e., no lower age limit.}
  \item{max.age}{Maximum age for MDA, in months. Defaults to \code{Inf}, i.e., no upper age limit.}
}

\note{
For more details and examples, please visit the \pkg{cystiSim} Wiki pages on \url{https://github.com/brechtdv/cystiSim/wiki}.
}

\author{\email{brechtdv@gmail.com}}

\seealso{Other interventions: \code{\link{do_pig_mda}}, \code{\link{do_pig_vac}}, \code{\link{do_pig_mda_vac}}}