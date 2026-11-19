# poster_topic23_timecourse.R -- poster figure for topics 2 + 3
# "The 3 early phosphosites: do they change without their protein (topic 2),
#  and are they still changed after ~15 h (topic 3)?"
#
# Run from the repository root (RStudio: open the project, then "Source"),
# or in the Terminal:  Rscript scripts/poster_topic23_timecourse.R
#
# Makes one figure in reports/ (PDF = vector, best for PowerPoint; PNG = 300 dpi):
#   poster_topic23_timecourse : one panel per site; change vs TSB after 15 min
#                               and ~15 h PASN, phosphosite vs its protein
#
# Same rules as topics 1-3: changed = FDR < 0.05 and |log2FC| > 1, observed estimate.
# Colours as in poster_topic1_fractions.R: burgundy = protein, rose = phosphosite
# (the rose is a bit stronger than in the bar chart so the thin lines stay visible).

source("R/setup.R")

BASE_SIZE <- 22

# ---- 1. The 3 early sites (from the topic 2 / topic 3 notebooks) --------------
early_sites <- tribble(
  ~locus_tag,     ~residue, ~position, ~panel,                 ~pattern,
  "B7H15_01690",  "S",      33,        "5'-nucleotidase S33",  "stays up",
  "B7H15_05655",  "T",      21,        "RpoY T21",             "early, transient",
  "B7H15_14290",  "S",      194,       "IsaA S194",            "follows its protein"
) |>
  mutate(panel_label = paste0(panel, "\n", pattern),
         panel_label = factor(panel_label, levels = panel_label))

contr <- c("TSB_vs_SNpulse_at_anc", "TSB_vs_SN_at_anc")

# ---- 2. Fold changes of sites and proteins ------------------------------------
fc <- bind_rows(
  read_dea("phospho_STYH", contrasts = contr) |>
    inner_join(early_sites, by = c("locus_tag", "residue", "position")) |>
    mutate(layer = "Phosphosite"),
  read_dea("total", contrasts = contr) |>
    inner_join(select(early_sites, locus_tag, panel_label), by = "locus_tag") |>
    mutate(layer = "Protein")
) |>
  mutate(time = if_else(contrast == "TSB_vs_SNpulse_at_anc", "15 min\nPASN", "~15 h\nPASN"),
         # course contrast TSB_vs_X: diff < 0 = higher after PASN.
         # Flip the sign so that positive = higher after PASN (easier to read).
         log2fc = -diff, lo = -conf.high, hi = -conf.low,
         changed = coalesce(estimate_type == "observed" & FDR < 0.05 & abs(diff) > 1, FALSE))

# check: all estimates observed? (lod_imputed = presence/absence only, not a fold change)
fc |> count(layer, estimate_type) |> print()
fc |> select(panel_label, layer, time, log2fc, FDR, changed) |> print(n = Inf)

# TSB is the reference: change = 0 by definition
tsb <- fc |> distinct(panel_label, layer) |>
  mutate(time = "TSB", log2fc = 0, lo = NA_real_, hi = NA_real_, changed = NA)

plot_df <- bind_rows(tsb, fc) |>
  mutate(time  = factor(time, levels = c("TSB", "15 min\nPASN", "~15 h\nPASN")),
         layer = factor(layer, levels = c("Protein", "Phosphosite")),
         # TSB start point and significant points: filled in the line colour;
         # not significant: open (white) point
         point_fill = case_when(is.na(changed) | changed ~ as.character(layer),
                                .default = "open"),
         point_size = if_else(time == "TSB", 3.5, 5.5))

# ---- 3. Figure -----------------------------------------------------------------
LAYER_COLOURS <- c("Protein" = "#6A1424", "Phosphosite" = "#E08E9B")
FILL_COLOURS  <- c(LAYER_COLOURS, "open" = "white")
dodge <- position_dodge(width = 0.15)

fold_breaks <- -4:4
fold_labels <- c("16× ↓", "8× ↓", "4× ↓", "2× ↓", "0",
                 "2× ↑", "4× ↑", "8× ↑", "16× ↑")

fig <- ggplot(plot_df, aes(time, log2fc, colour = layer, group = layer)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1, ymax = 1,
           fill = "grey95") +                          # < 2-fold = not changed
  geom_hline(yintercept = 0, colour = "grey55") +
  geom_line(position = dodge, linewidth = 1.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = dodge, width = 0.12,
                linewidth = 1, na.rm = TRUE) +
  geom_point(aes(fill = point_fill, size = point_size), position = dodge,
             shape = 21, stroke = 1.8) +
  scale_size_identity() +
  facet_wrap(~ panel_label, nrow = 1) +
  scale_colour_manual(values = LAYER_COLOURS, name = NULL) +
  scale_fill_manual(values = FILL_COLOURS, guide = "none") +
  scale_y_continuous(breaks = fold_breaks, labels = fold_labels) +
  coord_cartesian(ylim = c(-3.6, 4.3)) +
  labs(x = NULL, y = "change vs TSB\n(↑ higher / ↓ lower after PASN)",
       caption = paste("wild type JE2, n = 3  |  bars = 95 % CI  |  grey band = < 2-fold",
                       "|  filled = significant (FDR < 0.05, > 2-fold), open = not significant")) +
  theme_minimal(base_size = BASE_SIZE) +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2, "lines"),
        strip.text = element_text(face = "bold", size = BASE_SIZE * 0.95),
        axis.text.x = element_text(colour = "black"),
        plot.caption = element_text(colour = "grey40", size = BASE_SIZE * 0.7))

ggsave(file.path(REPO_ROOT, "reports", "poster_topic23_timecourse.pdf"), fig,
       width = 36, height = 16, units = "cm")
ggsave(file.path(REPO_ROOT, "reports", "poster_topic23_timecourse.png"), fig,
       width = 36, height = 16, units = "cm", dpi = 300, bg = "white")

message("Done: reports/poster_topic23_timecourse.pdf / .png")
