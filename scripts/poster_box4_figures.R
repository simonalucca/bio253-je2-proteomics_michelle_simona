# poster_box4_figures.R -- figures for poster box 4 ("which proteins change after 15 min?")
#
# Run from the repository root (RStudio: open the project, then "Source"),
# or in the Terminal:  Rscript scripts/poster_box4_figures.R
#
# Makes two figures in reports/ (PDF = vector, best for PowerPoint; PNG = 300 dpi):
#   A) poster_box4_A_lollipop  : all 34 early-changed proteins, sorted by fold change
#   A-short) poster_box4_A_short: only 12 key proteins, grouped by function (smaller)
#   B) poster_box4_B_rawvalues : raw per-replicate values of 4 key proteins
#
# Input: reports/15min_changed_proteins_anc.csv (made in notebook 08:
# wild type, TSB_vs_SNpulse_at_anc, FDR < 0.05 and |log2FC| > 1, reliable 2/3).

source("R/setup.R")

# Big fonts: the figure is printed at roughly this size on the A0 poster.
BASE_SIZE <- 22

# ---- 1. Hit list ---------------------------------------------------------
hits <- read_csv(file.path(REPO_ROOT, "reports", "15min_changed_proteins_anc.csv"),
                 show_col_types = FALSE) |>
  # the 35th row is an FGCZ lab contaminant (no locus tag) -> drop it
  filter(!is.na(locus_tag)) |>
  # course contrast is TSB minus SNpulse, so diff < 0 = higher after PASN.
  # We flip the sign so that positive = higher after PASN (easier to read).
  mutate(log2fc_pasn = -diff)
stopifnot(nrow(hits) == 34)

# ---- 2. Short names and function groups (from notebook 08, section 9) -----
# Functions taken from JE2 annotation / STRING v12; confirm in UniProt before printing.
groups <- tribble(
  ~staa8_uniprot, ~short,                   ~group,
  # fermentation and anaerobic metabolism
  "Q2G1K9", "AdhE",                         "Fermentation / anaerobic",
  "Q2FV02", "NrdD",                         "Fermentation / anaerobic",
  "Q2G1D7", "PflA",                         "Fermentation / anaerobic",
  "Q2G218", "Ldh1",                         "Fermentation / anaerobic",
  "Q2FV03", "NrdG",                         "Fermentation / anaerobic",
  "Q2G0G1", "Adh",                          "Fermentation / anaerobic",
  # cell wall / division
  "Q2FWF8", "SceD",                         "Other",
  "O33599", "LytM",                         "Other",
  "Q2FV52", "IsaA (box 5)",                 "Other",
  "Q2FZ86", "SepF",                         "Other",
  # stress and regulation
  "Q2G0D2", "AraC-family regulator",        "Other",
  "Q2FUS7", "Sulfurtransferase",            "Other",
  "Q2FUQ9", "Cold-shock protein",           "Other",
  "Q2G009", "Cold-shock protein 2",         "Other",
  "Q2FWD0", "RpoE",                         "Other",
  "Q2G282", "PerR",                         "Other",
  # detoxification
  "Q2FV27", "VOC protein",                  "Other",
  "Q2FY61", "Glyoxalase-family protein",    "Other",
  "Q2G1U2", "VOC protein 2",                "Other",
  # other metabolism / transport / secretion
  "Q2FXQ9", "HemA",                         "Other",
  "Q2FV34", "Type II secretion protein",    "Other",
  "Q2FV23", "PanD",                         "Other",
  "Q2G171", "BrnQ-like transporter",        "Other",
  "Q2G0E8", "DhaM",                         "Other",
  "Q2G178", "EssG",                         "Other",
  # unknown (up)
  "Q2FZA8", "Hypothetical protein 1",       "Other",
  "Q2G2C3", "Hypothetical protein 2",       "Other",
  # down after PASN
  "Q2FWM8", "Hld (delta-hemolysin)",        "Toxin / virulence",
  "Q2FW51", "Map-like protein",             "Toxin / virulence",
  "Q2FW65", "Hypothetical protein 3",       "Other",
  "Q2G1V8", "Hypothetical protein 4",       "Other",
  "Q2G1H1", "Isochorismatase-family prot.", "Other",
  "Q2G2H8", "MnhE1",                        "Other",
  "Q2G2M0", "Tautomerase",                  "Other"
)

hits <- left_join(hits, groups, by = "staa8_uniprot")
stopifnot(!any(is.na(hits$short)))   # every protein got a name
print(count(hits, direction))         # should say 27 higher, 7 lower

GROUP_COLOURS <- c("Fermentation / anaerobic" = "#7A1F2B",  # poster dark red
                   "Toxin / virulence"        = "#1B7F79",  # teal
                   "Other"                    = "grey65")

# ---- 3. Figure A: lollipop of all 34 proteins -----------------------------
fold_breaks <- c(-3, -2, -1, 0, 1, 2, 3, 4)
fold_labels <- c("8x", "4x", "2x", "0", "2x", "4x", "8x", "16x")

fig_a <- hits |>
  mutate(short = reorder(short, log2fc_pasn),
         group = factor(group, levels = names(GROUP_COLOURS))) |>
  ggplot(aes(log2fc_pasn, short, colour = group)) +
  geom_vline(xintercept = 0, colour = "grey30") +
  geom_segment(aes(x = 0, xend = log2fc_pasn, yend = short), linewidth = 1.2) +
  geom_point(size = 4.5) +
  annotate("text", x = -3.4, y = 34.6, label = "lower after PASN", hjust = 0,
           size = BASE_SIZE / 3.2, colour = "grey30") +
  annotate("text", x = 4.6, y = 1, label = "higher after PASN", hjust = 1,
           size = BASE_SIZE / 3.2, colour = "grey30") +
  scale_colour_manual(values = GROUP_COLOURS, name = NULL) +
  scale_x_continuous(breaks = fold_breaks, labels = fold_labels,
                     limits = c(-3.5, 4.7)) +
  labs(x = "fold change after 15 min PASN vs TSB", y = NULL) +
  theme_minimal(base_size = BASE_SIZE) +
  theme(legend.position = "top",
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = "black"))

ggsave(file.path(REPO_ROOT, "reports", "poster_box4_A_lollipop.pdf"), fig_a,
       width = 26, height = 34, units = "cm")
ggsave(file.path(REPO_ROOT, "reports", "poster_box4_A_lollipop.png"), fig_a,
       width = 26, height = 34, units = "cm", dpi = 300, bg = "white")

# ---- 4. Figure B: raw values of 4 key proteins (course rule: raw values first) ----
key <- tribble(
  ~locus_tag,     ~panel,
  "B7H15_00850",  "AdhE (fermentation)",
  "B7H15_01295",  "Ldh1 (fermentation)",
  "B7H15_11280",  "Hld (toxin)",
  "B7H15_14290",  "IsaA (cell wall, box 5)"
)

raw <- read_abundances("total") |>
  filter(genotype == "anc") |>
  inner_join(key, by = "locus_tag") |>
  mutate(panel  = factor(panel, levels = key$panel),
         medium = factor(medium, levels = MEDIA,
                         labels = c("TSB", "15 min", "~15 h")))

# check: 3 of 3 replicates detected everywhere? (NA = not detected, never 0)
raw |> group_by(panel, medium) |>
  summarise(n_detected = sum(!is.na(normalized_abundance)), .groups = "drop") |>
  print(n = Inf)

MEDIUM_COLOURS <- c("TSB" = "grey55", "15 min" = "#3461C1", "~15 h" = "#D65A31")

fig_b <- ggplot(raw, aes(medium, normalized_abundance, colour = medium)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.6,
               na.rm = TRUE) +
  geom_point(size = 4, position = position_jitter(width = 0.08, height = 0, seed = 1),
             na.rm = TRUE) +
  facet_wrap(~panel, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = MEDIUM_COLOURS, guide = "none") +
  labs(x = NULL, y = "abundance (log2)") +
  theme_bw(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

ggsave(file.path(REPO_ROOT, "reports", "poster_box4_B_rawvalues.pdf"), fig_b,
       width = 40, height = 13, units = "cm")
ggsave(file.path(REPO_ROOT, "reports", "poster_box4_B_rawvalues.png"), fig_b,
       width = 40, height = 13, units = "cm", dpi = 300, bg = "white")


# ---- 5. Figure A-short: only the key proteins (for a smaller poster box) ----
# 12 of the 34 proteins, chosen by biology (the story of box 4), not only by size.
# 'role' = very short function, printed next to the bar (check in UniProt!).
key_short <- tribble(
  ~staa8_uniprot, ~block,             ~role,
  "Q2G1K9", "Fermentation",           "ethanol",               # AdhE
  "Q2G1D7", "Fermentation",           "formate",               # PflA
  "Q2G218", "Fermentation",           "lactate",               # Ldh1
  "Q2G0G1", "Fermentation",           "ethanol",               # Adh
  "Q2FV02", "Fermentation",           "DNA synthesis w/o O2",  # NrdD
  "Q2FV03", "Fermentation",           "activates NrdD",        # NrdG
  "Q2FWF8", "Cell wall",              "remodelling",           # SceD
  "O33599", "Cell wall",              "autolysin",             # LytM
  "Q2FV52", "Cell wall",              "autolysin",             # IsaA (box 5)
  "Q2G282", "Oxidative stress",       "H2O2 sensor",           # PerR
  "Q2FWM8", "Virulence",              "toxin (agr)",           # Hld
  "Q2FW51", "Virulence",              "immune evasion"         # Map-like
)
BLOCK_COLOURS <- c("Fermentation" = "#7A1F2B", "Cell wall" = "#5B4A8B",
                   "Oxidative stress" = "#6B6B6B", "Virulence" = "#1B7F79")

short_df <- hits |>
  inner_join(key_short, by = "staa8_uniprot") |>
  mutate(block = factor(block, levels = names(BLOCK_COLOURS)),
         short = str_remove(short, " \\(.*\\)"),        # "IsaA (box 5)" -> "IsaA"
         short = if_else(staa8_uniprot == "Q2FV52", "IsaA *", short),
         short = if_else(staa8_uniprot == "Q2FWM8", "Hld", short),
         short = reorder(short, log2fc_pasn),
         up    = log2fc_pasn > 0,
         fc_text = paste0(fold_change, "\u00d7"))            # e.g. "18.2x"

# group label in the top-left (or top-right for the 'down' block) of each panel
block_labels <- short_df |>
  group_by(block) |>
  summarise(up = all(up), .groups = "drop") |>
  mutate(x = if_else(up, -4.4, 5.3), hjust = if_else(up, 0, 1))

fig_a_short <- ggplot(short_df, aes(log2fc_pasn, short)) +
  # light tinted background per function group
  geom_rect(data = block_labels, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = block),
            alpha = 0.08) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.6) +
  geom_segment(aes(x = 0, xend = log2fc_pasn, yend = short, colour = block),
               linewidth = 2.2, lineend = "round") +
  geom_point(aes(fill = block), shape = 21, colour = "white", stroke = 1.2, size = 7) +
  # fold change printed at the end of each bar
  geom_text(aes(label = fc_text, hjust = if_else(up, -0.35, 1.35), colour = block),
            size = BASE_SIZE / 3.4, fontface = "bold") +
  # group name inside each panel
  geom_text(data = block_labels, inherit.aes = FALSE,
            aes(x = x, y = Inf, label = block, hjust = hjust, colour = block),
            vjust = 1.4, size = BASE_SIZE / 3, fontface = "bold") +
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = BLOCK_COLOURS, guide = "none") +
  scale_fill_manual(values = BLOCK_COLOURS, guide = "none") +
  scale_x_continuous(breaks = fold_breaks, labels = fold_labels,
                     limits = c(-4.5, 5.4)) +   # extra room for the fold-change labels
  scale_y_discrete(expand = expansion(add = c(0.7, 1.4))) +   # room for group name
  labs(x = "fold change after 15 min PASN vs TSB\n(left = lower, right = higher after PASN)",
       y = NULL,
       caption = "12 of 34 changed proteins shown  |  * phosphosite in box 5") +
  theme_minimal(base_size = BASE_SIZE) +
  theme(panel.grid = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(colour = "grey35"),
        plot.caption = element_text(colour = "grey45", size = BASE_SIZE * 0.6),
        panel.spacing.y = unit(0.5, "lines"))

ggsave(file.path(REPO_ROOT, "reports", "poster_box4_A_short.pdf"), fig_a_short,
       width = 28, height = 22, units = "cm")
ggsave(file.path(REPO_ROOT, "reports", "poster_box4_A_short.png"), fig_a_short,
       width = 28, height = 22, units = "cm", dpi = 300, bg = "white")

message("Done: figures written to reports/poster_box4_*")
