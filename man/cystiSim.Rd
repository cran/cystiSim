\name{cystiSim}
\alias{cystiSim}
\alias{report}
\alias{elim}
\alias{print.cystiSim}
\alias{summary.cystiSim}
\alias{plot.cystiSim}
\alias{report.cystiSim}
\alias{elim.cystiSim}

\title{'cystiSim' object}

\description{
Functions to initiate and explore 'cystiSim' objects. A 'cystiSim' object corresponds to multiple runs of a \pkg{cystiSim} model.}

\usage{
cystiSim(n = 100, mod, main = NULL)

report(x, \dots)

\S3method{print}{cystiSim}(x, \dots)

\S3method{summary}{cystiSim}(object, round = 3, \dots)

\S3method{plot}{cystiSim}(x, y, annotate = TRUE, \dots)

\S3method{report}{cystiSim}(x, name = "cystiSim", \dots)

\S3method{elim}{cystiSim}(x, show = c("m", "y"), \dots)
}

\arguments{
  \item{n}{Number of iterations (months).}
  \item{mod}{cystiSim model.}
  \item{main}{cystiSim model label.}
  
  \item{x}{Object of class 'cystiSim'.}
  \item{object}{Object of class 'cystiSim'.}
  \item{y}{Currently ignored.}
  \item{round}{Number of decimal digits to be printed.}
  \item{annotate}{Should plot be annotated with summary information?}
  \item{name}{Report name.}
  \item{show}{Show time till elimination in terms of months or years?}
  \item{\dots}{Other arguments to be passed to generics.}
}

%\details{}

\value{
\code{cystiSim} returns an object of S3 class 'cystiSim'.
}

\note{
For more details and examples, please visit the \pkg{cystiSim} Wiki pages on \url{https://github.com/brechtdv/cystiSim/wiki}.
}

\author{\email{brechtdv@gmail.com}}

\seealso{
\code{\link{cystiRun}}
}

\examples{
\dontrun{

## we will use the built-in Mbeya dataset
prevalence(pig_mbeya$cysti)
prevalence(man_mbeya$taenia)

## define transmission probabilities
ph2m <- 0.000174918
pl2m <- 0.000149501
m2p  <- 6.85E-05
e2p  <- 0.00022611

## set seed for reproducibility
set.seed(264)

## need to define coverage and efficacy of all interventions
cov_man_mda <- 0.80
cov_pig_mda <- 0.90
cov_pig_vac <- NULL
eff_man_mda <- 0.70
eff_pig_mda <- 0.90
eff_pig_vac <- NULL

## run the simulations
sim <-
cystiSim(
  n = 10,
  main = "example",
  mod = {
initiate(man_mbeya, pig_mbeya, ph2m, pl2m, m2p, e2p) \%>\%
  update(240) \%>\%
  do_man_mda(coverage = cov_man_mda, efficacy = eff_man_mda) \%>\%
  do_pig_mda(coverage = cov_pig_mda, efficacy = eff_pig_mda) \%>\%
  update(120)
  }
)

## summarize results
summary(sim)

## plot simulations (mean and uncertainty interval)
plot(sim)

## create PDF report and plot
report(sim)
}
}