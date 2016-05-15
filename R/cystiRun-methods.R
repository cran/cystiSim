### CYSTIRUN METHODS

print.cystiRun <-
function(x, from = 200, to = NA, ...) {
  if (is.na(to)) to <- nrow(x$out)
  p <- colMeans(x$out[seq(from, to), ])
  cat("PC: ", round(p[1], 3),
      "| Pr(H):", round(p[2], 3), "\n")
  cat("PCI:", round(p[4], 3), "\n")
  cat("HT: ", round(p[6], 4),
      "| HT(inf):", round(p[7], 4),
      ", HT(chl):", round(p[8], 4),
      ", HT(adl):", round(p[9], 4), "\n")
  cat("E:  ", round(p[10], 4), "\n")
}

plot.cystiRun <-
function(x, y = NULL, show = c("PC", "PR", "HT", "EN"),
         start = 0, from = 1, to = NA, ...) {
  ## define census labels and groups
  label <- c("PC", "P(H)", "P(L)", "PCI", "PR",
             "HT", "HT(inf)", "HT(ch)", "HT(ad)", "EN")
  group <- c(rep("pig", 5), rep("man", 5))

  ## select labels and groups to show
  id_col <- match(show, label)
  lab <- factor(label[id_col], levels = label[id_col])
  grp <- factor(group[id_col], levels = unique(group[id_col]))

  ## define extraction points
  if (is.na(to)) to <- nrow(x$out)
  id_row <- seq(from, to)

  ## extract census data; force matrix
  out <- as.matrix(x$out[id_row, id_col])

  ## define months
  months <- seq(nrow(out)) - start

  ## create ggplot data.frame
  df <- data.frame(p = c(out),
                   m = rep(months, times = ncol(out)),
                   lab = rep(lab, each = nrow(out)),
                   grp = rep(grp, each = nrow(out)))

  ## build plot function
  ggplot(df, aes_string(x = "m", y = "p")) +
    geom_line(aes(col = lab)) +
    scale_colour_manual(values = seq(10)[id_col]) +
    facet_grid(grp ~ ., scales = "free", as.table = FALSE) +
    scale_x_continuous("month") +
    scale_y_continuous("proportion") +
    theme_bw()
}
