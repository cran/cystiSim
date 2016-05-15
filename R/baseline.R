### CYSTISIM: AGENT BASED TAENIA SOLIUM TRANSMISSION MODEL

## -------------------------------------------------------------------------#
## GENERATE RANDOM BASELINE ------------------------------------------------#

## simple random baseline for humans
random_baseline_man <-
function(n, p) {
  data.frame(age = round(runif(n, 0, 80) * 12),
             sex = sample(c("male", "female"), n, replace = TRUE),
             taenia = rbinom(n, 1, p),
             taenia_immature = 0,
             time_since_infection = 0,
             environment = 0,
             time_since_contamination = 0)
}

## slaughter function - default
slaughter_nbinom <-
function(age, min, max, size, mu) {
  ## if age below 'min', do not kill
  kill <- rep(0, length(age))

  ## if age larger or equal 'max', always kill
  kill[age >= max] <- 1

  ## if age in between, kill probabilistically
  kill_age <- age >= min & age < max
  kill[kill_age] <-
    rbinom(sum(kill_age), 1, pnbinom(age[kill_age] - min, size, mu = mu))

  ## return results
  return(kill)
}

## slaughter function - simplified
slaughter_binom <-
function(age, min, max, p) {
  ## if age below 'min', do not kill
  kill <- rep(0, length(age))

  ## if age larger or equal 'max', always kill
  kill[age >= max] <- 1

  ## if age in between, kill probabilistically
  kill_age <- age >= min & age < max
  kill[kill_age] <- rbinom(sum(kill_age), 1, p)

  ## return results
  return(kill)
}

## model pig age structure - default
pig_age_model <-
function(n, steps, size, mu) {
  ## create population with 'n' births per month
  pigs <- data.frame(age = rep(0, n))

  for (i in seq(5)) {
    ## ageing of pigs
    pigs$age <- pigs$age + 1

    ## new piglets get born
    pigs <- rbind(pigs, data.frame(age = rep(0, n)))
  }

  ## define progress bar
  pb <- txtProgressBar(max = steps, style = 3)

  ## run through cycles
  for (i in seq(steps)) {
    ## ageing of pigs
    pigs$age <- pigs$age + 1

    ## slaughter of pigs - default
    kill <- slaughter_nbinom(pigs$age, min = 6, max = 36, size, mu)

    ## remove slaughtered pigs from population
    pigs <- subset(pigs, !kill)

    ## new piglets get born
    ## number of births equal to number of killed pigs
    pigs <- rbind(pigs, data.frame(age = rep(0, sum(kill))))

    ## update progress bar
    setTxtProgressBar(pb, i)
  }

  ## close progress bar
  close(pb)

  ## return pigs dataframe
  return(pigs)
}

## random baseline pigs
random_baseline_pig <-
function(n, p, p.high) {
  ## model age structure of pig population
  pigs_age <- pig_age_model(n / 6, 500, 0.70, 80)$age

  ## identify pigs that are 3 months or older
  ## only these pigs are old enough to appear Ag positive
  pigs_inf <- pigs_age > 3

  ## randomly infect pigs
  cysti <- rep(0, length(pigs_inf))
  cysti[pigs_inf] <- rbinom(sum(pigs_inf), 1, p)

  ## randomly assign infection intensity to infected pigs
  intensity <- rep(0, length(pigs_inf))
  intensity[cysti == 1] <-
    sample(c("H", "L"), sum(cysti), TRUE, c(p.high, 1 - p.high))

  ## create 'pigs' dataframe
  pigs <-
    data.frame(age = pigs_age,
               cysti = cysti,
               cysti_immature = 0,
               time_since_infection = 0,
               intensity = intensity,
               immunity = 0,
               time_since_vaccination = NA,
               slaughtered = 0)

  ## return 'pigs' dataframe
  return(pigs)
}