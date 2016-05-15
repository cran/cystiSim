## CYSTISIM INTERVENTIONS

## HUMANS - MASS DRUG ADMINISTRATION
do_man_mda <-
function(x, coverage, efficacy, min.age = 0, max.age = Inf) {
  ## identify HT that are eligible for treatment
  id <- (x$man$taenia == 1 | x$man$taenia_immature == 1) &
        x$man$age >= min.age & x$man$age <= max.age

  ## randomize treatment over eligible HT
  is_treated <- rbinom(sum(id), 1, coverage * efficacy) == 1

  ## both mature and immature tapeworms die
  x$man$taenia[id][is_treated] <- 0
  x$man$taenia_immature[id][is_treated] <- 0

  ## 'time_since_infection' gets set to zero
  x$man$time_since_infection[id][is_treated] <- 0

  ## 'time_since_contamination' gets set to one: initiate env contam decay
  x$man$time_since_contamination[id][is_treated] <- 1

  ## return updated 'cystiRun' object
  return(x)
}


## PIGS - MASS DRUG ADMINISTRATION
do_pig_mda <-
function(x, coverage, efficacy, immunity = 3, min.age = 1, max.age = Inf) {
  ## identify pigs that are eligible for treatment
  id <- x$pig$slaughtered == 0 &
        x$pig$age >= min.age & x$pig$age <= max.age

  ## randomize treatment over eligible pigs
  is_treated <- rbinom(sum(id), 1, coverage * efficacy) == 1

  ## generate immunity for treated positive pigs (mature & immature)
  ## .. if immunity is Inf (cf vaccination), no need to reset
  is_immune <-
    is.finite(x$pig$immunity[id][is_treated]) &
    (x$pig$cysti[id][is_treated] == 1 |
     x$pig$cysti_immature[id][is_treated] == 1)
  x$pig$immunity[id][is_treated][is_immune] <- immunity

  ## both mature and immature cysts die
  ## .. 'intensity' and 'time_since_infection' get set to zero
  x$pig$cysti[id][is_treated] <- 0
  x$pig$cysti_immature[id][is_treated] <- 0
  x$pig$intensity[id][is_treated] <- 0
  x$pig$time_since_infection[id][is_treated] <- 0

  ## return updated 'cystiRun' object
  return(x)
}


## PIGS - VACCINATION
do_pig_vac <-
function(x, coverage, efficacy, immunity = Inf, interval = 4,
         min.age = 1, max.age = Inf) {
  ## so far only lifelong immunity implemented !!!
  if (is.finite(immunity)) {
    stop("Currently only lifelong immunity is implemented.")
  }

  ## identify pigs that are eligible for vaccination
  id <- x$pig$slaughtered == 0 &
        x$pig$age >= min.age & x$pig$age <= max.age

  ## randomize vaccination over eligible pigs
  is_vaccinated <- rbinom(sum(id), 1, coverage * efficacy) == 1

  ## generate immunity for pigs that were vaccinated < XXX months before
  is_immune <- !is.na(x$pig$time_since_vaccination[id][is_vaccinated]) &
               x$pig$time_since_vaccination[id][is_vaccinated] > 0
  x$pig$immunity[id][is_vaccinated][is_immune] <- immunity

  ## reset counter for pigs that were vaccinated > XXX months before
  is_reset <- !is.na(x$pig$time_since_vaccination[id][is_vaccinated]) &
               x$pig$time_since_vaccination[id][is_vaccinated] <= 0
  x$pig$time_since_vaccination[id][is_vaccinated][is_reset] <- interval + 1

  ## initiate 'time_since_vaccination' if first vaccination
  is_first <- is.na(x$pig$time_since_vaccination[id][is_vaccinated])
  x$pig$time_since_vaccination[id][is_vaccinated][is_first] <- interval + 1

  ## return updated 'cystiSim' object
  return(x)
}


## PIGS - MDA + VACCINATION
do_pig_mda_vac <-
function(x, coverage, efficacy.mda, efficacy.vac,
         immunity.mda = 3, immunity.vac = Inf, interval = 4,
         min.age = 1, max.age = Inf) {
  ## so far only lifelong immunity implemented !!!
  if (is.finite(immunity.vac)) {
    stop("Currently only lifelong immunity is implemented.")
  }

  ## identify pigs that are eligible for treatment
  id <- x$pig$slaughtered == 0 &
        x$pig$age >= min.age & x$pig$age <= max.age

  ## randomize coverage over eligible pigs
  ## .. covered pigs receive both MDA and VAC
  ## .. hence perfect correlation between MDA and VAC coverage
  is_mda <- rbinom(sum(id), 1, coverage) == 1
  is_vac <- is_mda

  ## randomize effective MDA treatment over covered pigs
  is_mda_eff <-
    rbinom(sum(is_mda), 1, efficacy.mda) == 1

  ## randomize effective VAC treatment over covered pigs
  is_vac_eff <-
    rbinom(sum(is_vac), 1, efficacy.vac) == 1

  ## MDA

  ## generate immunity for treated positive pigs (mature & immature)
  ## .. if immunity is Inf (cf vaccination), no need to reset
  is_immune <-
    is.finite(x$pig$immunity[id][is_mda][is_mda_eff]) &
    (x$pig$cysti[id][is_mda][is_mda_eff] == 1 |
     x$pig$cysti_immature[id][is_mda][is_mda_eff] == 1)
  x$pig$immunity[id][is_mda][is_mda_eff][is_immune] <- immunity.mda

  ## both mature and immature cysts die
  ## .. 'intensity' and 'time_since_infection' get set to zero
  x$pig$cysti[id][is_mda][is_mda_eff] <- 0
  x$pig$cysti_immature[id][is_mda][is_mda_eff] <- 0
  x$pig$intensity[id][is_mda][is_mda_eff] <- 0
  x$pig$time_since_infection[id][is_mda][is_mda_eff] <- 0

  ## VAC

  ## generate immunity for pigs that were vaccinated < XXX months before
  ## .. note that VAC immunity has precendence over MDA immunity
  is_immune <-
    !is.na(x$pig$time_since_vaccination[id][is_vac][is_vac_eff]) &
    x$pig$time_since_vaccination[id][is_vac][is_vac_eff] > 0
  x$pig$immunity[id][is_vac][is_vac_eff][is_immune] <-
    immunity.vac

  ## reset counter for pigs that were vaccinated > XXX months before
  is_reset <-
    !is.na(x$pig$time_since_vaccination[id][is_vac][is_vac_eff]) &
    x$pig$time_since_vaccination[id][is_vac][is_vac_eff] <= 0
  x$pig$time_since_vaccination[id][is_vac][is_vac_eff][is_reset] <-
    interval + 1

  ## initiate 'time_since_vaccination' if first vaccination
  is_first <-
    is.na(x$pig$time_since_vaccination[id][is_vac][is_vac_eff])
  x$pig$time_since_vaccination[id][is_vac][is_vac_eff][is_first] <-
    interval + 1

  ## return updated 'cystiRun' object
  return(x)
}
