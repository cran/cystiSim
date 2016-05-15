### CYSTISIM: AGENT BASED TAENIA SOLIUM TRANSMISSION MODEL

## -------------------------------------------------------------------------#
## HELPER FUNCTIONS --------------------------------------------------------#

## expit = inverse logit
expit <-
function(x) {
  exp(x)/(1 + exp(x))
}

## reporter function: prevalence
prevalence <-
function(z) {
  sum(z == 1, na.rm = TRUE) / length(z[!is.na(z)])
}

## reporter function: census
census <-
function(x) {
    ## porcine cysticercosis (mature) prevalence
  c(prevalence(x$pig$cysti),

    ## proportion high infection intensity
    prevalence(x$pig$intensity[x$pig$cysti == 1] == "H"),

    ## proportion low infection intensity
    prevalence(x$pig$intensity[x$pig$cysti == 1] == "L"),

    ## porcine cysticercosis (immature) prevalence
    prevalence(x$pig$cysti_immature),

    ## proportion immune pigs
    prevalence(x$pig$immunity > 0),

    ## taeniosis (mature) prevalence
    prevalence(x$man$taenia),

    ## taeniosis (mature) prevalence, infants
    prevalence(x$man$taenia[x$man$age <= 3*12]),

    ## taeniosis (mature) prevalence, children
    prevalence(x$man$taenia[x$man$age > 3*12 & x$man$age < 16*12]),

    ## taeniosis (mature) prevalence, adults
    prevalence(x$man$taenia[x$man$age >= 16*12]),

    ## environmental contamination prevalence
    prevalence(x$man$environment))
}

## exponential decay of environmental contamination
exp_decay <-
function(time) {
  p_decay <- rep(0.2682422, length(time))
  p_decay[time >= 9] <- 1     # 100% decay at month 9
  return(p_decay)
}


## -------------------------------------------------------------------------#
## TRANSMISSION FUNCTIONS --------------------------------------------------#

infect_man <-
function(man, pig, ph2m, pl2m, age.coef) {
  ## calculate number of Pigs Infected / High intensity
  PIH <- with(pig, sum(cysti == 1 & slaughtered == 1 & intensity == "H"))

  ## calculate number of Pigs Infected / Low intensity
  PIL <- with(pig, sum(cysti == 1 & slaughtered == 1 & intensity == "L"))

  ## identify Susceptible Humans
  HS <- with(man,
             taenia == 0 & taenia_immature == 0 &  # without tapeworm
             age >= 24)                            # no infants < 2yo

  ## randomly create new Human Tapeworm carriers
  ## probability of HT through any pig
  ## == 1 - probability of no infection through any pig
  p2m <- 1 - (1 - (1 - (1 - ph2m) ^ PIH)) *
             (1 - (1 - (1 - pl2m) ^ PIL))

  ## 'p2m' is the average probability of infection
  ## disaggregate 'p2m' according to age

  ## maximum age in population
  age_max <- max(man$age[HS] / 12)

  ## age-specific odds of infection defined from logistic regression
  odds <- expit(age.coef[1] + seq(0, age_max) * age.coef[2])

  ## tabulate age distribution of current dataframe
  age_dist <- table(cut(man$age[HS] / 12, seq(0, age_max + 1), right = FALSE))

  ## define ratio of odds and age-weighted average
  ratio <- odds / weighted.mean(odds, age_dist)

  ## define age-specific probability of HT
  p2m_age <- p2m * ratio[1 + (man$age[HS] / 12)]

  ## now randomly create new Human Tapeworm carriers
  ## based on age-specific transmission probability
  new_HT <- rbinom(sum(HS), 1, p2m_age)

  ## update 'taenia_immature' identifier
  man$taenia_immature[HS] <- new_HT

  ## update 'time_since_infection' identifier (redundant)
  man$time_since_infection[HS][new_HT == 1] <- 0

  ## return updated dataframe
  return(man)
}

infect_pig <-
function(pig, man, m2p, e2p) {

  ## DIRECT TRANSMISSION

  ## calculate number of Human Tapeworm carriers
  HT <- with(man, sum(taenia))

  ## identify Susceptible Pigs
  PS <- with(pig,
             cysti == 0 & cysti_immature == 0 &     # no cysts
             immunity == 0 &                        # no immunity
             slaughtered == 0)                      # not slaughtered

  ## randomly create new PC High intensity cases (from Susceptible cases)
  new_PCS2H <- rbinom(sum(PS), 1, 1 - (1 - m2p) ^ HT)

  ## identify Low intensity Pigs
  PSL <- with(pig,
              (cysti == 1 | cysti_immature == 1) &  # cysts
              immunity == 0 &                       # no immunity (redundant)
              slaughtered == 0 &                    # not slaughtered
              intensity == "L")                     # low intensity

  ## randomly create new PC High intensity cases (from Low intensity cases)
  new_PCL2H <- rbinom(sum(PSL), 1, 1 - (1 - m2p) ^ HT)


  ## INDIRECT TRANSMISSION

  ## calculate number of contaminated ENvironments
  EN <- with(man, sum(environment))

  ## randomly create new PC Low intensity cases (from Susceptible cases)
  new_PCS2L <- rbinom(sum(PS), 1, 1 - (1 - e2p) ^ EN)


  ## UPDATE IDENTIFIERS

  ## susceptible pigs becoming infected through HT or EN
  new_PC <- 1 - ((1 - new_PCS2L) * (1 - new_PCS2H))

  ## update 'cysti_immature' identifier
  pig$cysti_immature[PS] <- new_PC

  ## update 'intensity' identifier
  pig$intensity[PS][new_PCS2L == 1] <- "L"   # new Low intensity case
  pig$intensity[PS][new_PCS2H == 1] <- "H"   # new High intensity case
  pig$intensity[PSL][new_PCL2H == 1] <- "H"  # update to High from Low

  ## update 'time_since_infection' identifier (redundant)
  pig$time_since_infection[PS][new_PC == 1] <- 0

  ## return updated dataframe
  return(pig)
}


## -------------------------------------------------------------------------#
## INITIATE MODEL ----------------------------------------------------------#

initiate <-
function(man, pig,
         ph2m, pl2m, m2p, e2p, age.coef = c(0, 0),
         slaughter = slaughter_nbinom,
         slaughter.args = list(min = 6, max = 36, size = 0.70, mu = 80)) {
  ## create 'cystiRun' object
  x <-
  list(man = man,
       pig = pig,
       out = NULL,
       par = list(age_coef = age.coef,
                  ph2m = ph2m,
                  pl2m = pl2m,
                  m2p = m2p,
                  e2p = e2p),
       slaughter = list(f = slaughter,
                        args = slaughter.args))

  ## take baseline population census
  x$out <- matrix(census(x), ncol = 10)

  ## add S3 class 'cystiRun'
  class(x) <- "cystiRun"

  ## return 'cystiRun' object
  return(x)
}


## -------------------------------------------------------------------------#
## UPDATE MODEL ------------------------------------------------------------#

update.cystiRun <-
function(object, n = 1200, verbose = TRUE, ...) {

## MODEL
x <- object
out <- matrix(nrow = n, ncol = 10)

if (verbose) pb <- txtProgressBar(max = n, style = 3)

for (i in seq(n)) {

  ## STEP 01: maturation of tapeworms present in previous cycle
  man_infected <- x$man$taenia == 1 | x$man$taenia_immature == 1
  x$man$time_since_infection[man_infected] <-
    x$man$time_since_infection[man_infected] + 1

  ## STEP 02: immature tapeworms become infectious after 2 months
  ##          and starts contaminating environment
  man_infectious <- x$man$time_since_infection == 2
  x$man$taenia[man_infectious] <- 1             ## mature HT
  x$man$taenia_immature[man_infectious] <- 0    ## no immature HT
  x$man$environment[man_infectious] <- 1        ## envir contam

  ## STEP 03: tapeworms die after 13 months and human becomes susceptible
  tw_dies <- x$man$time_since_infection == 13
  x$man$taenia[tw_dies] <- 0                     # tw is gone
  x$man$time_since_infection[tw_dies] <- 0       # reset time since inf
  x$man$time_since_contamination[tw_dies] <- 1   # envir contamination!

  ## STEP 04: exponential decay of environmental contamination
  ## .. if no HT & time.since.contam > 0, exp decay
  ## .. if no HT & time.since.contam > 0 & no decay, time.since.contam++
  is_decay <- x$man$time_since_contamination > 0
  p_decay  <- exp_decay(x$man$time_since_contamination[is_decay])
  decay    <- rbinom(sum(is_decay), 1, p_decay) == 1
  x$man$time_since_contamination[is_decay][decay] <- 0
  x$man$environment[is_decay][decay] <- 0  ### SET ENVIRONMENT TO 0
  x$man$time_since_contamination[is_decay][!decay] <-
    x$man$time_since_contamination[is_decay][!decay] + 1

  ## STEP 05: new HT, based on slaughtered infectious pigs of previous cycle
  man_new <- infect_man(x$man, x$pig,
                        x$par$ph2m, x$par$pl2m, x$par$age_coef)

  ## STEP 06: remove slaughtered pigs of previous cycle
  x$pig <- subset(x$pig, x$pig$slaughtered == 0)

  ## STEP 07: maturation of cysts present in previous cycle
  is_infected <- x$pig$cysti_immature == 1
  x$pig$time_since_infection[is_infected] <-
    x$pig$time_since_infection[is_infected] + 1

  ## STEP 08: immature cysts become infectious after 3 months
  is_infectious <- x$pig$time_since_infection == 3
  x$pig$cysti[is_infectious] <- 1
  x$pig$cysti_immature[is_infectious] <- 0
  x$pig$time_since_infection[is_infectious] <- 0

  ## STEP 09: new PI, based on HT & ENVIR of previous cycle
  pig_new <- infect_pig(x$pig, x$man, x$par$m2p, x$par$e2p)

  ## STEP 10: slaughter pigs of previous cycle
  is_slaughtered <-
    do.call(x$slaughter$f, c(list(age = pig_new$age), x$slaughter$args))
  pig_new$slaughtered <- is_slaughtered

  ## STEP 11: ageing of pigs not slaughtered in current cycle
  pig_new$age <- pig_new$age + 1

  ## STEP 12: waning of immunity of pigs not slaughtered in current cycle
  is_immune <- pig_new$immunity > 0
  pig_new$immunity[is_immune] <-
    pig_new$immunity[is_immune] - 1

  ## STEP 13: update 'time_since_vaccination' counter for vaccinated pigs
  is_vaccinated <- !is.na(pig_new$time_since_vaccination)
  pig_new$time_since_vaccination[is_vaccinated] <-
    pig_new$time_since_vaccination[is_vaccinated] - 1

  ## STEP 14: birth of new pigs
  pig_new <-
    rbind(pig_new,
          data.frame(age = rep(0, sum(is_slaughtered)),
                     cysti = 0,
                     cysti_immature = 0,
                     time_since_infection = 0,
                     intensity = 0,
                     immunity = 0,
                     time_since_vaccination = NA,
                     slaughtered = 0))

  ## replace data frames with new ones
  x$pig <- pig_new
  x$man <- man_new

  ### take populaton census
  out[i, ] <- census(x)

  ## update progress bar
  if (verbose) setTxtProgressBar(pb, i)
}

## close progress bar
if (verbose) close(pb)

## combine old and new census
x$out <- rbind(x$out, out)

## return updated 'cystiRun' object
return(x)
}
