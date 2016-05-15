### CYSTISIM: AGENT BASED TAENIA SOLIUM TRANSMISSION MODEL

## -------------------------------------------------------------------------#
## FIT MODEL PARAMETERS ----------------------------------------------------#

fit <-
function(n.sim, n.update, target, limit,
         man, pig, ph2m, pl2m, m2p, e2p, age.coef = c(0, 0),
         slaughter = slaughter_nbinom,
         slaughter.args = list(min = 6,  max = 36, size = 0.70, mu = 80)) {
  ## evaluate 'target'
  if (!is.list(target)) stop("'target' must be a list")
  sapply(names(target), match.arg, c("ht", "pc", "pi"))

  ## predefine 
  out <- matrix(ncol = 16, nrow = 0)
  colnames(out) <- c("ph2m", "pl2m", "m2p", "e2p", 
                     "PC", "P(H)", "P(L)", "PCI", "PI",
                     "HT", "HT(inf)", "HT(ch)", "HT(ad)", "EN", "HT(tot)",
                     "DEV")

  ## initiate counter
  i <- 0

  ## initiate progress bar
  pb <- txtProgressBar(max = n.sim, style = 3)

  ## run simulations
  while (i < n.sim) {
    ph2m_sim <- runif(1, ph2m[1], ph2m[2])
    pl2m_sim <- runif(1, pl2m[1], pl2m[2])
    m2p_sim  <- runif(1, m2p[1], m2p[2])
    e2p_sim  <- runif(1, e2p[1], e2p[2])

    mod <-
      initiate(man, pig,
               ph2m_sim, pl2m_sim, m2p_sim, e2p_sim, age.coef,
               slaughter, slaughter.args) %>%
      update(n.update, verbose = FALSE)

    ht <- prevalence(mod$man$taenia_immature + mod$man$taenia)
    pc <- prevalence(mod$pig$cysti)
    pi <- sum(mod$pig$intensity == "H") / sum(mod$pig$intensity != "0")

    ht_diff <- ifelse("ht" %in% names(target), ht-target$ht, 0)
    pc_diff <- ifelse("pc" %in% names(target), pc-target$pc, 0)
    pi_diff <- ifelse("pi" %in% names(target), pi-target$pi, 0)
    if (is.nan(pi_diff)) pi_diff <- 0

    dev <- ht_diff^2 + pc_diff^2 + pi_diff^2

    if (dev < limit) {
      i <- i + 1
      out <- rbind(out,
                   c(ph2m_sim, pl2m_sim, m2p_sim, e2p_sim,
                     c(tail(mod$out, 1)), ht, dev))
      setTxtProgressBar(pb, i)
    }
  }

  cat("\n")

  return(out)
}
