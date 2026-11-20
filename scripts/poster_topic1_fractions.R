# poster_topic1_fractions.R -- poster figure for topic 1
# "Do phosphosites change earlier than proteins after PASN exposure?"
#
# Run from the repository root (RStudio: open the project, then "Source"),
# or in the Terminal:  Rscript scripts/poster_topic1_fractions.R
#
# Makes one figure in reports/ (PDF = vector, best for PowerPoint; PNG = 300 dpi):
#   poster_topic1_fractions : share of reliably measured proteins / phosphosites
#                             that changed after 15 min and after ~15 h PASN
#
# Same definitions as the combined topic 1 notebook
# (analyses_15min_combined/01_combined_topic 1_final.qmd):
#   reliable = detected in >= 2/3 replicates in both groups AND estimate "observed"
#   changed  = FDR < 0.05 and |log2FC| > 1
# Colours match the poster (box 4 script and the BioRender method figure):
#   15 min = blue, ~15 h = orange.

source("R/setup.R")

BASE_SIZE <- 22   # big fonts: printed at roughly this size on the A0 poster

contr <- c("TSB_vs_SNpulse_at_anc", "TSB_vs_SN_at_anc")

# ---- 1. Flag every protein / site as reliable and changed ------------------
# key = column that identifies one feature: protein_key or site
flag_changes <- function(dataset, key) {
  # how many of the 3 replicates were measured per group (NA = not detected)
  n_obs <- read_abundances(dataset) |>
    filter(genotype == "anc") |>
    group_by(across(all_of(c(key, "medium")))) |>
    summarise(n = sum(!is.na(normalized_abundance)), .groups = "drop") |>
    pivot_wider(names_from = medium, values_from = n, values_fill = 0)

  read_dea(dataset, contrasts = contr) |>
    left_join(n_obs, by = key) |>
    mutate(
      reliable = case_when(
        contrast == "TSB_vs_SNpulse_at_anc" ~ TSB >= 2 & SNpulse >= 2,
        contrast == "TSB_vs_SN_at_anc"      ~ TSB >= 2 & SN >= 2,
        .default = FALSE),
      # lod_imputed estimates are presence/absence evidence, not a fold change
      reliable = coalesce(reliable, FALSE) &
                 coalesce(estimate_type == "observed", FALSE),
      changed  = coalesce(FDR < 0.05 & abs(diff) > 1, FALSE))
}

prot <- flag_changes("total", key = "protein_key") |>
  # FGCZ lab contaminants (zz_...) are not S. aureus proteins -> drop them
  filter(!str_starts(protein_key, "zz_"))
phos <- flag_changes("phospho_STYH", key = "site")

# ---- 2. Count: changed / reliably measured ----------------------------------
fractions <- bind_rows(
  prot |> mutate(layer = "Proteins"),
  phos |> mutate(layer = "Phosphosites")
) |>
  filter(reliable) |>
  group_by(layer, contrast) |>
  summarise(n_reliable = n(), n_changed = sum(changed), .groups = "drop") |>
  mutate(fraction = n_changed / n_reliable,
         exposure = factor(if_else(contrast == "TSB_vs_SNpulse_at_anc",
                                   "15 min PASN pulse", "~15 h PASN"),
                           levels = c("15 min PASN pulse", "~15 h PASN")),
         layer = factor(layer, levels = c("Proteins", "Phosphosites")),
         label = paste0(sprintf("%.1f", 100 * fraction), " %\n(",
                        n_changed, "/", n_reliable, ")"))
print(fractions)

# ---- 3. Test at 15 min: is the phosphosite share larger? ---------------------
early <- fractions |> filter(exposure == "15 min PASN pulse") |>
  arrange(desc(layer))                       # row 1 = phosphosites
early_counts <- cbind(changed     = early$n_changed,
                      not_changed = early$n_reliable - early$n_changed)
fisher_p <- fisher.test(early_counts, alternative = "greater")$p.value
message("one-sided Fisher test (15 min, phosphosites > proteins): p = ",
        signif(fisher_p, 2))

# ---- 4. Figure ---------------------------------------------------------------
# One colour per data layer (same in both exposures): burgundy tones,
# dark = proteins, light rose = phosphosites (stronger contrast for the printed poster).
LAYER_COLOURS <- c("Proteins" = "#6A1424", "Phosphosites" = "#EDB8BE")

ymax <- max(fractions$fraction) * 1.3
dodge <- position_dodge(width = 0.8)

fig <- ggplot(fractions, aes(exposure, fraction, fill = layer)) +
  geom_col(position = dodge, width = 0.75, colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label), position = dodge, vjust = -0.3,
            size = BASE_SIZE / 3.0, lineheight = 0.9, colour = "grey15",
            fontface = "bold") +
  scale_fill_manual(values = LAYER_COLOURS, name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, ymax),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "changed / reliably measured",
       caption = "wild type JE2, n = 3 per group  |  FDR < 0.05, > 2-fold, detected in \u2265 2/3 replicates") +
  theme_minimal(base_size = BASE_SIZE) +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        # bold throughout: the poster is read from ~1.5 m away
        axis.text.x = element_text(colour = "black", face = "bold"),
        axis.text.y = element_text(colour = "grey20", face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey30", face = "bold",
                                    size = BASE_SIZE * 0.8))

ggsave(file.path(REPO_ROOT, "reports", "poster_topic1_fractions.pdf"), fig,
       width = 24, height = 18, units = "cm")
ggsave(file.path(REPO_ROOT, "reports", "poster_topic1_fractions.png"), fig,
       width = 24, height = 18, units = "cm", dpi = 300, bg = "white")

message("Done: reports/poster_topic1_fractions.pdf / .png")
