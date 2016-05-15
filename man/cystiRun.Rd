\name{cystiRun}
\alias{cystiRun}
\alias{initiate}
\alias{update.cystiRun}
\alias{print.cystiRun}
\alias{plot.cystiRun}
\alias{prevalence}

\title{'cystiRun' object}

\description{
Functions to initiate, update and explore 'cystiRun' objects. A 'cystiRun' object corresponds to a single run of a \pkg{cystiSim} model.
}

\usage{
initiate(man, pig,
         ph2m, pl2m, m2p, e2p, age.coef = c(0, 0),
         slaughter = slaughter_nbinom,
         slaughter.args = list(min = 6, max = 36, size = 0.70, mu = 80))

\S3method{update}{cystiRun}(object, n = 1200, verbose = TRUE, \dots)

\S3method{print}{cystiRun}(x, from = 200, to = NA, \dots)

\S3method{plot}{cystiRun}(x, y = NULL, show = c("PC", "PR", "HT", "EN"),
      start = 0, from = 1, to = NA, \dots)
      
prevalence(z)
}

\arguments{
  \item{man}{Human population dataframe.}
  \item{pig}{Pig population dataframe.}
  \item{ph2m}{Pig (Heavy infection) to Man transmission probability.}
  \item{pl2m}{Pig (Light infection) to Man transmission probability.}
  \item{m2p}{Man to Pig transmission probability.}
  \item{e2p}{Environment to Pig transmission probability.}
  \item{age.coef}{Optional intercept and slope for the association between age and taeniosis.}
  \item{slaughter}{Function that defines the slaughter probability of pigs.}
  \item{slaughter.args}{Arguments to be passed to the slaughter function.}

  \item{object}{Object of class 'cystiRun'.}
  \item{n}{Number of iterations (months).}
  \item{verbose}{If \code{TRUE}, a progress bar is shown.}

  \item{x}{Object of class 'cystiRun'.}
  \item{y}{Currently ignored.}
  \item{z}{Vector of infection indicators.}
  
  \item{start}{Origin of plot.}
  \item{from}{First iteration to be used in output.}
  \item{to}{Last iteration to be used in output. The default value \code{NA} corresponds to the last available iteration.}
  \item{show}{Which output should be plotted?}
  
  \item{\dots}{Arguments to be passed on to generic function.}
}

\value{
\code{initiate} and \code{update} return an object of S3 class 'cystiRun'.
}

\note{
For more details and examples, please visit the \pkg{cystiSim} Wiki pages on \url{https://github.com/brechtdv/cystiSim/wiki}.
}

\author{\email{brechtdv@gmail.com}}

\seealso{
Intervention functions: \code{\link{do_man_mda}}, \code{\link{do_pig_mda}}, \code{\link{do_pig_vac}}, \code{\link{do_pig_mda_vac}}
}

\examples{
## we will use the built-in Mbeya dataset
prevalence(pig_mbeya$cysti)
prevalence(man_mbeya$taenia)

## define transmission probabilities
ph2m <- 0.000174918
pl2m <- 0.000149501
m2p  <- 6.85E-05
e2p  <- 0.00022611

## first initiate the 'cystiRun' object
mod <- initiate(man_mbeya, pig_mbeya, ph2m, pl2m, m2p, e2p)

## update the model 240 cycles (=20 years)
## this is a burn-in period, needed to obtain steady state
mod <- update(mod, 240)

## apply human mass drug administration
mod <- do_man_mda(mod, efficacy = 0.70, coverage = 0.80)

## apply pig mass drug administration
mod <- do_pig_mda(mod, efficacy = 0.90, coverage = 0.90)

## update the model 120 more cycles (=10 years)
mod <- update(mod, 120)

## plot the cycles
plot(mod, from = 200, start = 40)
}