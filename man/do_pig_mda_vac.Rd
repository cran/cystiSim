\name{do_pig_mda_vac}
\alias{do_pig_mda_vac}

\title{Do pig mass drug administration (MDA) and vaccination}

\description{
Intervention function that mimics the possible effects of combined pig mass drug administration and vaccination.
}

\usage{
do_pig_mda_vac(x, coverage, efficacy.mda, efficacy.vac,
               immunity.mda = 3, immunity.vac = Inf, interval = 4,
               min.age = 1, max.age = Inf)
}

\arguments{
  \item{x}{\code{cystiRun} object.}
  \item{coverage}{Presumed coverage within eligible population, expressed as a decimal value.}
  \item{efficacy.mda}{Presumed efficacy of drug treatment, expressed as a decimal value.}
  \item{efficacy.vac}{Presumed efficacy of vaccine, expressed as a decimal value.}
  \item{immunity.mda}{Presumed duration of immunity following drug treatment of an infected pig. Defaults to 3 months.}
  \item{immunity.vac}{Presumed duration of immunity following successful vaccination. Defaults to \code{Inf}.}
  \item{interval}{Maximum interval between two consecutive shots for immunity. Defaults to 4 months.}
  \item{min.age}{Minimum age for intervention, in months. Defaults to 1 month.}
  \item{max.age}{Maximum age for intervention, in months. Defaults to \code{Inf}, i.e., no upper age limit.}
}

\note{
For more details and examples, please visit the \pkg{cystiSim} Wiki pages on \url{https://github.com/brechtdv/cystiSim/wiki}.
}

\author{\email{brechtdv@gmail.com}}

\seealso{Other interventions: \code{\link{do_man_mda}}, \code{\link{do_pig_mda}}, \code{\link{do_pig_vac}}}