\name{fit}
\alias{fit}

\title{Fit cystiSim parameters}

\description{
This function simulates random \pkg{cystiSim} parameters, and return those scenarios that result in a baseline prevalence close to the specified target.
}

\usage{
fit(n.sim, n.update, target, limit,
    man, pig, ph2m, pl2m, m2p, e2p, age.coef = c(0, 0),
    slaughter = slaughter_nbinom,
    slaughter.args = list(min = 6,  max = 36, size = 0.70, mu = 80))
}

\arguments{
  \item{n.sim}{Desired number of retained scenarios.}
  \item{n.update}{Number of updates of the baseline model.}
  \item{target}{Named list of target prevalences for \code{ht}, \code{pc} and/or \code{pi}}
  \item{limit}{Maximum tolerated deviance.}
  
  \item{man}{Human population dataframe.}
  \item{pig}{Pig population dataframe.}
  \item{ph2m}{Pig (Heavy infection) to Man transmission probability.}
  \item{pl2m}{Pig (Light infection) to Man transmission probability.}
  \item{m2p}{Man to Pig transmission probability.}
  \item{e2p}{Environment to Pig transmission probability.}
  \item{age.coef}{Optional intercept and slope for the association between age and taeniosis.}
  \item{slaughter}{Function that defines the slaughter probability of pigs.}
  \item{slaughter.args}{Arguments to be passed to the slaughter function.}
}

\details{
The abbreviations used:
\itemize{
  \item{\code{ht}....}{Human Taeniosis prevalence}
  \item{\code{pc}....}{Porcine Cysticercosis prevalence}
  \item{\code{pi}....}{Pig Intensity proportion (i.e., proportion heavily infection pigs)}
}
}

\note{
For more details and examples, please visit the \pkg{cystiSim} Wiki pages on \url{https://github.com/brechtdv/cystiSim/wiki}.
}

\author{\email{brechtdv@gmail.com}}
