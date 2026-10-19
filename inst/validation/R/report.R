# Figures and the markdown report of a validation run.

# reference palette (light): categorical slots, surface, ink, grid
val_pal <- list(
  series = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100"),
  surface = "#fcfcfb", ink = "#0b0b0b", ink2 = "#52514e", grid = "#e6e5e1",
  band = "#dcdad3"
)

val_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = val_pal$surface, colour = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = val_pal$grid, linewidth = 0.4),
      text = ggplot2::element_text(colour = val_pal$ink),
      axis.text = ggplot2::element_text(colour = val_pal$ink2),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      legend.position = "top", legend.justification = "left",
      plot.title.position = "plot"
    )
}

# 6.1: rejection rate per design and term, with and without the propensity
# model, against the Monte Carlo band around 0.05
val_fig_calibration <- function(rates) {
  d <- rates
  d$design <- factor(sprintf("K=%d\nT=%d", d$K, d$T),
    levels = unique(sprintf("K=%d\nT=%d", sort(d$K), d$T)[order(d$T, d$K)])
  )
  d$term <- factor(d$term, c("group", "time", "interaction"))
  d$model <- factor(ifelse(d$propensity, "with propensity_formula", "without"),
    c("with propensity_formula", "without")
  )
  ggplot2::ggplot(d, ggplot2::aes(.data$design, .data$rejection)) +
    ggplot2::annotate("rect",
      xmin = -Inf, xmax = Inf, ymin = d$lo[1], ymax = d$hi[1],
      fill = val_pal$band, alpha = 0.6
    ) +
    ggplot2::geom_hline(yintercept = 0.05, colour = val_pal$ink2, linewidth = 0.4, linetype = 2) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$model, shape = .data$model),
      size = 2.6, position = ggplot2::position_dodge(width = 0.5)
    ) +
    # y range per panel: 0-1 for the confounded group term, 0-0.15 (or the
    # data) for the within-subject terms, so the band reads as a band
    ggplot2::geom_blank(data = data.frame(
      term = factor(c("group", "time", "interaction"), levels(d$term)),
      design = vapply(c("group", "time", "interaction"), function(tm) {
        as.character(d$design[d$term == tm][1])
      }, ""),
      rejection = c(1, max(0.15, d$rejection[d$term == "time"]),
        max(0.15, d$rejection[d$term == "interaction"]))
    )) +
    ggplot2::scale_colour_manual(values = val_pal$series[1:2], name = NULL) +
    ggplot2::scale_shape_manual(values = c(16, 17), name = NULL) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::facet_wrap(~term, scales = "free") +
    ggplot2::labs(
      x = NULL, y = "Rejection rate at alpha = 0.05",
      title = "Calibration under a confounded null",
      subtitle = paste(
        "Shaded: Monte Carlo band around 0.05.",
        "Group labels depend on age; positivity depends on age only."
      )
    ) +
    val_theme()
}

# 6.2: pooled power of the four per-peptide statistics
val_fig_power <- function(power) {
  d <- power
  d$stat <- factor(d$stat, c("score", "srlr", "srlr_std", "wald_t"),
    c("score z", "srlr", "standardised srlr", "Wald t")
  )
  d$alt <- factor(d$alt, c("null", "up", "down"),
    c("null", "group 2 higher", "group 2 lower")
  )
  d$sizes <- factor(paste("n =", d$sizes), paste("n =", c("20/20", "30/10", "10/30")))
  ggplot2::ggplot(d, ggplot2::aes(.data$alt, .data$rejection)) +
    ggplot2::geom_hline(yintercept = 0.05, colour = val_pal$ink2, linewidth = 0.4, linetype = 2) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$stat, shape = .data$stat),
      size = 2.6, position = ggplot2::position_dodge(width = 0.6)
    ) +
    ggplot2::scale_colour_manual(values = val_pal$series, name = NULL) +
    ggplot2::scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::facet_wrap(~sizes) +
    ggplot2::labs(
      x = NULL, y = "Rejection rate at alpha = 0.05",
      title = "Pooled power by per-peptide statistic",
      subtitle = paste(
        "Same data and permutations; 500 sparse peptides, 20% shifted.",
        "Group 2 is the second size."
      )
    ) +
    val_theme()
}

# 6.5: estimated vs true rho
val_fig_rho <- function(recovery) {
  ggplot2::ggplot(recovery, ggplot2::aes(.data$rho, .data$mean)) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, colour = val_pal$ink2, linewidth = 0.4, linetype = 2
    ) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$mean - .data$sd, ymax = .data$mean + .data$sd),
      width = 0.02, colour = val_pal$series[1], linewidth = 0.5
    ) +
    ggplot2::geom_point(colour = val_pal$series[1], size = 2.6) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "True within-subject correlation", y = "Pooled estimate (mean, +/- 1 SD)",
      title = "Recovery of rho", subtitle = "60 subjects, 1-6 samples each, 200 peptides"
    ) +
    val_theme() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(colour = val_pal$grid, linewidth = 0.4)
    )
}

val_md_table <- function(df, digits = 3) {
  df <- as.data.frame(df)
  for (j in seq_along(df)) {
    if (is.numeric(df[[j]])) df[[j]] <- formatC(df[[j]], digits = digits, format = "fg")
    if (is.logical(df[[j]])) df[[j]] <- ifelse(df[[j]], "yes", "no")
  }
  c(
    paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  )
}

val_write_report <- function(results, out_dir, settings) {
  gates <- dplyr::bind_rows(lapply(results, `[[`, "gates"))
  lines <- c(
    "# compute_delta2() validation",
    "",
    sprintf("Run: %s. %s", format(Sys.time(), "%Y-%m-%d %H:%M"), settings),
    "",
    "## Gates",
    "",
    sprintf("**%d of %d gates pass.**", sum(gates$pass), nrow(gates)),
    "",
    val_md_table(dplyr::transmute(gates,
      section = .data$section, gate = .data$gate, value = .data$value,
      lo = .data$lo, hi = .data$hi, pass = .data$pass
    ))
  )
  if (!is.null(results$s1)) {
    lines <- c(lines, "", "## 6.1 Calibration per term", "",
      "![calibration](calibration.png)", "",
      val_md_table(results$s1$rates[c("K", "T", "propensity", "term", "n", "rejection")]))
  }
  if (!is.null(results$s2)) {
    lines <- c(lines, "", "## 6.2 Pooled power by statistic", "",
      "![power](pooled_power.png)", "",
      val_md_table(results$s2$power[
        c("sizes", "alt", "stat", "n", "rejection", "pooled_null_mean")
      ]),
      "", "Exact one-peptide permutation moments:", "",
      val_md_table(results$s2$moments))
  }
  if (!is.null(results$s3)) {
    lines <- c(lines, "", "## 6.3 Unbalanced designs", "", "(a) visits correlated with group:", "",
      val_md_table(results$s3$rates_a), "", "(b) informative dropout:", "",
      val_md_table(results$s3$rates_b))
  }
  if (!is.null(results$s4)) {
    lines <- c(lines, "", "## 6.4 Interaction robustness", "", val_md_table(results$s4$rates))
  }
  if (!is.null(results$s5)) {
    lines <- c(lines, "", "## 6.5 rho", "", "![rho](rho_recovery.png)", "",
      val_md_table(results$s5$recovery), "", "Power by rho option:", "",
      val_md_table(results$s5$power))
    if (!is.null(results$s5$strata)) {
      lines <- c(lines, "", "rho per prevalence stratum, phiper example data (phip_mixture):", "",
        val_md_table(results$s5$strata))
    }
  }
  if (!is.null(results$tests)) {
    lines <- c(lines, "", "## 6.6 / 6.7 Automated tests", "", val_md_table(results$tests$table))
  }
  writeLines(lines, file.path(out_dir, "report.md"))
}
