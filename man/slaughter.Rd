\name{slaughter}
\alias{slaughter_binom}
\alias{slaughter_nbinom}

\title{Pig slaughter functions}

\description{
These functions simulate age-dependent slaughter of pigs. The default function is \code{slaughter_nbinom}.
}

\usage{
slaughter_binom(age, min, max, p)
slaughter_nbinom(age, min, max, size, mu)
}

\arguments{
  \item{age}{Age of the pigs.}
  \item{min}{Minimum age at slaughter, i.e., \emph{Pr(slaughter<min)=0}.}
  \item{max}{Age at which all pigs are definitely slaughtered, i.e., \emph{Pr(slaughter>=max)=1}.}
  \item{p}{Binomial probability of slaughter.}
  \item{size}{Size of Negative Binomial distribution of age-specific slaughter probability.}
  \item{mu}{Mean of Negative Binomial distribution of age-specific slaughter probability.}
}

\note{
For more details and examples, please visit the \pkg{cystiSim} Wiki pages on \url{https://github.com/brechtdv/cystiSim/wiki}.
}

\author{\email{brechtdv@gmail.com}}
