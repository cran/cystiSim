### CYSTISIM: AGENT BASED TAENIA SOLIUM TRANSMISSION MODEL

## -------------------------------------------------------------------------#
## MAIN SIMULATION FUNCTION ------------------------------------------------#

cystiSim <-
function(n = 100, mod, main = NULL) {
  ## check if global variables are present
  ph2m <- get("ph2m", envir = .GlobalEnv)
  pl2m <- get("pl2m", envir = .GlobalEnv)
  m2p  <- get("m2p",  envir = .GlobalEnv)
  e2p  <- get("e2p", envir = .GlobalEnv)
  cov_man_mda <- get("cov_man_mda", envir = .GlobalEnv)
  cov_pig_mda <- get("cov_pig_mda", envir = .GlobalEnv)
  cov_pig_vac <- get("cov_pig_vac", envir = .GlobalEnv)
  eff_man_mda <- get("eff_man_mda", envir = .GlobalEnv)
  eff_pig_mda <- get("eff_pig_mda", envir = .GlobalEnv)
  eff_pig_vac <- get("eff_pig_vac", envir = .GlobalEnv)

  ## predefine list of simulations
  out <- vector("list", n)

  ## perform simulations
  for (i in seq(n)) {
    print(i)
    out[[i]] <- eval(substitute(mod))$out[, c(1, 5, 6)]
  }

  ## create 'cystiSim' object
  sim <-
    list(out = out,
         mod = substitute(mod),
         main = main,
         par = list(ph2m = ph2m,
                    pl2m = pl2m,
                    m2p = m2p,
                    e2p = e2p,
                    cov_man_mda = cov_man_mda,
                    cov_pig_mda = cov_pig_mda,
                    cov_pig_vac = cov_pig_vac,
                    eff_man_mda = eff_man_mda,
                    eff_pig_mda = eff_pig_mda,
                    eff_pig_vac = eff_pig_vac))

  ## add S3 class 'cystiSim'
  class(sim) <- "cystiSim"

  ## return 'cystiSim' object
  return(sim)
}


## -------------------------------------------------------------------------#
## cystiSim METHODS -------------------------------------------------------#

print.cystiSim <-
function(x, ...) {
  mod <- strsplit(as.character(x$mod)[2], "%>%")[[1]]
  mod <- gsub("\n", "", mod)
  mod <- sub("^ +", "", mod)
  mod <- sub(" +$", "", mod)
  for (i in seq_along(mod)) {
    if (i > 1) cat ("  ")
    cat(mod[i])
    if (i < length(mod)) cat(" %>%\n")
  }
  cat("\n")
}

summary.cystiSim <-
function(object, round = 3, ...) {
  sim_mx <- sapply(object$out, as.matrix)
  m <- nrow(sim_mx) / ncol(object$out[[1]])

  out <-
  rbind(
    c(mean(sim_mx[m, ] == 0), quantile(sim_mx[m, ], 0.95)),      # PC
    c(mean(sim_mx[3*m, ] == 0), quantile(sim_mx[3*m, ], 0.95)))  # HT

  rownames(out) <- c("PC", "HT")
  colnames(out) <- c("Pr(elim)", "95%")

  return(round(out, round))
}

plot.cystiSim <-
function(x, y = NULL, annotate = TRUE, ...) {
  lab <- factor(c("PC", "PR", "HT"), c("PC", "PR", "HT"))
  grp <- factor(c("pig", "pig", "human"), c("pig", "human"))
  col <- c(1, 3, 2)

  sim_mx <- sapply(x$out, as.matrix)  # ncol=it; nrow=par*months

  df <-
  data.frame(mean = rowMeans(sim_mx),
             lwr  = apply(sim_mx, 1, quantile, 0.025),
             upr  = apply(sim_mx, 1, quantile, 0.975),
             m = rep(seq(nrow(sim_mx) / ncol(x$out[[1]])),
                     ncol(x$out[[1]])),
             par = rep(lab, each = nrow(sim_mx) / ncol(x$out[[1]])),
             grp = rep(grp, each = nrow(sim_mx) / ncol(x$out[[1]])))

  stats <- summary(x)

  x_txt <- nrow(sim_mx) / ncol(x$out[[1]])
  y_txt <- tapply(apply(sim_mx, 1, max), df$grp, max)

  df_txt <-
    data.frame(lab = paste0(ncol(sim_mx), " simulations",
                            "\nPr(elim)=", stats[, 1]),
               grp = levels(df$grp),
               x = x_txt,
               y = y_txt)

  g <-
  ggplot(df, aes_string(x = "m", y = "mean")) +
    geom_ribbon(aes_string(ymin = "lwr", ymax = "upr", fill = "par"),
                alpha = .25) +
    geom_line(aes_string(y = "mean", col = "par"), size = 1) +
    geom_line(aes_string(y = "lwr", col = "par")) +
    geom_line(aes_string(y = "upr", col = "par")) +
    scale_fill_manual(values = col) +
    scale_colour_manual(values = col) +
    scale_x_continuous("month") +
    scale_y_continuous("prevalence") +
    facet_grid(grp~., scales = "free") +
    theme(legend.position = "none") +
    theme_bw() +
    ggtitle(x$main)

  if (annotate) {
    g <-
      g +
      geom_text(data = df_txt,
                aes_string(label = "lab", x = "x", y = "y"),
                size = 4, hjust = 1, vjust = 1)
  }

  return(g)
}

report <-
function(x, ...) {
  UseMethod("report")
}

report.cystiSim <-
function(x, name = "cystiSim", ...) {
  ## temporary 'mod' file
  ## write model
  tmp <- tempfile(fileext = ".Rnw")
  sink(tmp)
  cat("<<eval=FALSE>>=\n")  
  print(x)
  cat("@\n")
  sink()

  ## knit PDF
  knit2pdf(system.file("cystiSim.Rnw", package = "cystiSim"))

  ## rename PDF
  file.rename("cystiSim.pdf",
              paste0(name, "_", today(), ".pdf"))

  ## remove helper files
  unlink(tmp)
  file.remove(paste0("cystiSim.",
                     c("tex", "log", "aux")))
  unlink("figure", recursive = TRUE)

  ## generate PNG
  png(paste0(name, "_", today(), ".png"),
      width = 10, height = 5, units = "in", res = 300)
  print(plot(x))
  graphics.off()
}

elim <-
function(x, ...) {
  UseMethod("elim")
}

elim.cystiSim <-
function(x, show = c("m", "y"), ...) {
  show <- match.arg(show)
  cat(nrow(x$out[[1]]), "iterations\n")
  pc <- sapply(x$out, function(x) tail(which(x[, 1] != 0), 1)) - 201
  ht <- sapply(x$out, function(x) tail(which(x[, 3] != 0), 1)) - 201

  div <- ifelse(show == "m", 1, 12)

  Summary <-
  function(x) {
    c(mean = mean(x),
      min = min(x),
      max = max(x),
      quantile(x, c(0.025, 0.975)))
  }

  rbind(pc = Summary(pc), ht = Summary(ht)) / div
}

## -------------------------------------------------------------------------#
## HELPER FUNCTIONS --------------------------------------------------------#

today <-
function() {
  return(format(Sys.time(), "%Y%m%d"))
}