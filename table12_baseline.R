source('R/setup.R',local=TRUE)

raw <- read_excel(input_file, sheet = "Baseline characteristics")

required <- c("Subject ID", "Group", "Gender", "Age at enrolment",
              "Eosinophils", "IgE", "FeNO")

stopifnot(all(required %in% names(raw)))

allergen_col <- grep("^Specific IgE", names(raw), value = TRUE)

stopifnot(length(allergen_col) == 1L)

baseline <- data.frame(
  id = trimws(raw[["Subject ID"]]),
  group = trimws(raw[["Group"]]),
  sex = toupper(trimws(raw[["Gender"]])),
  age = raw[["Age at enrolment"]],
  eos = raw[["Eosinophils"]],
  ige = raw[["IgE"]],
  feno = raw[["FeNO"]],
  allergen = trimws(raw[[allergen_col]])
)

baseline <- baseline[!is.na(baseline$id) & nzchar(baseline$id), ]

stopifnot(!anyDuplicated(baseline$id), !anyDuplicated(infected_ids))

stopifnot(all(infected_ids %in% baseline$id))

dat <- baseline[baseline$id %in% infected_ids, ]

stopifnot(setequal(dat$id,infected_ids))

stopifnot(all(is.na(dat$sex) | dat$sex %in% c("F", "M")))

for (v in c("age", "eos", "ige", "feno")) {
  stopifnot(is.numeric(dat[[v]]), all(is.na(dat[[v]]) | is.finite(dat[[v]])))
}

dat$sensitised <- ifelse(is.na(dat$allergen) | dat$allergen == "", NA,
                         tolower(dat$allergen) != "all negative")

mean_sd <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("Not available; n=0")
  sprintf("%.1f ± %.1f; n=%d", mean(x), sd(x), length(x))
}

median_iqr <- function(x, digits = 1, unit = "") {
  x <- x[!is.na(x)]
  if (!length(x)) return("Not available; n=0")
  # R's default type 7 quantiles. IQR is displayed as Q1–Q3.
  q <- quantile(x, c(0.25, 0.5, 0.75), type = 7, names = FALSE)
  f <- function(z) formatC(z, format = "f", digits = digits)
  paste0(f(q[2]), unit, " (IQR ", f(q[1]), "–", f(q[3]), "); n=", length(x))
}

count_percent <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("Not available; n=0")
  sprintf("%d/%d (%.1f%%)", sum(x), length(x), 100 * mean(x))
}

summarise_group <- function(g) {
  d <- dat[dat$group == g, ]
  c(mean_sd(d$age), count_percent(d$sex == "F"), count_percent(d$sex == "M"),
    median_iqr(d$eos, 2), median_iqr(d$ige, 1),
    if (g == "Control" && all(is.na(d$feno))) "Not measured; n=0" else
      median_iqr(d$feno, 1, " ppb"),
    count_percent(d$sensitised))
}

table1 <- data.frame(
  Characteristic = c("Age at enrolment, years", "Female sex", "Male sex",
                     "Blood eosinophils, ×10⁹/L", "Total IgE, IU/mL",
                     "FeNO", "Allergen sensitisation"),
  Asthma = summarise_group("Asthma"),
  Control = summarise_group("Control"), check.names = FALSE
)

names(table1)[2:3] <- sprintf("%s (N = %d)", c("Asthma", "Control"),
                             c(sum(dat$group == "Asthma"), sum(dat$group == "Control")))

draw_table <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  # Wider first column, as in the original reference.
  xs <- c(0.012, 0.540, 0.850)
  top <- 0.982
  bottom <- 0.018
  edges <- seq(top, bottom, length.out = nrow(table1) + 2)
  centres <- (head(edges, -1) + tail(edges, -1)) / 2
  for (i in seq_along(edges)) {
    grid::grid.lines(x = c(0.006, 0.994), y = rep(edges[i], 2),
      gp = grid::gpar(col = "#D0D0D0", lwd = if (i %in% c(1, 2, length(edges))) 1.6 else 0.8))
  }
  for (j in 1:3) {
    grid::grid.text(names(table1)[j], x = xs[j], y = centres[1],
      just = if (j == 1) "left" else "centre",
      gp = grid::gpar(fontfamily = "Arial", fontsize = 21,
                      fontface = "bold", col = "#383838"))
    for (i in seq_len(nrow(table1))) {
      grid::grid.text(table1[i, j], x = xs[j], y = centres[i + 1],
        just = if (j == 1) "left" else "centre",
        gp = grid::gpar(fontfamily = "Arial", fontsize = 20,
                        fontface = if (j == 1) "bold" else "plain", col = "#383838"))
    }
  }
}

ragg::agg_png(file.path(output_dir, "Table12_baseline_characteristics.png"),
              width = 16.6, height = 6.05, units = "in", res = 300)

draw_table()

invisible(dev.off())

