\name{cystiSim-package}
\alias{cystiSim-package}

\docType{package}

\title{Agent-based model for \emph{Taenia solium} transmission and control.}

\description{The cystiSim model allows simulating \emph{Taenia solium} taeniosis/cysticercosis transmission in a virtual population of humans and pigs. It also allows evaluating the possible effects of human mass drug administration, pig mass drug administration, and pig vaccination. \pkg{cystiSim} was developed within the framework of CYSTINET, the European Network on taeniosis/cysticercosis, COST ACTION TD1302 (\url{http://www.cystinet.org/}).}

\details{
  \tabular{ll}{
    Package:    \tab  cystiSim\cr
    Version:    \tab  0.1.0\cr
    Date:       \tab  2016-05-15\cr
    Authors:    \tab  Brecht Devleesschauwer, Uffe Christian Braae\cr
    Maintainer: \tab  \email{brechtdv@gmail.com}\cr
    URL:        \tab  https://github.com/brechtdv/cystiSim\cr
    BugReports: \tab  https://github.com/brechtdv/cystiSim/issues\cr
    Depends:    \tab  R (>= 3.3.0), ggplot2\cr
    Imports:    \tab  magrittr, knitr\cr
    License:    \tab  GPL (>= 2)\cr
    LazyLoad:   \tab  yes\cr
}

Available functions in the \pkg{cystiSim} package:
  \tabular{ll}{
    \strong{baseline}\cr
    \code{\link{random_baseline_man}}  \tab Generate a random baseline human population.\cr
    \code{\link{random_baseline_pig}}  \tab Generate a random baseline pig population.\cr\cr
    
    \strong{model---S3 class '\code{\link{cystiRun}}'}\cr
    \code{\link{fit}}             \tab Fit parameters of a \code{cystiRun} model.\cr
    \code{\link{initiate}}        \tab Initiate a \code{cystiRun} model.\cr
    \code{\link{update}}          \tab Update a \code{cystiRun} model.\cr\cr

    \strong{interventions}\cr
    \code{\link{do_man_mda}}      \tab Do human mass drug administration.\cr
    \code{\link{do_pig_mda}}      \tab Do pig mass drug administration.\cr
    \code{\link{do_pig_vac}}      \tab Do pig vaccination.\cr
    \code{\link{do_pig_mda_vac}}  \tab Do pig mass drug administration AND vaccination.\cr\cr

    \strong{simulate---S3 class '\code{\link{cystiSim}}'}\cr
    \code{\link{cystiSim}}        \tab Simulate multiple \code{cystiRun} models.\cr
    \code{\link{report}}          \tab Generate a PDF report for a \code{cystiSim} object.\cr
    \code{\link{elim}}            \tab Derive when elimination occurred in a \code{cystiSim} object.\cr\cr
  }

For more details and examples, please visit the \pkg{cystiSim} Wiki pages on \url{https://github.com/brechtdv/cystiSim/wiki}.
}

\keyword{package}
