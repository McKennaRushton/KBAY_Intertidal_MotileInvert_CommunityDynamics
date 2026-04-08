# BIOL 462 Manuscript: A Decade of Intertidal Community Dynamics in Alaska's Kachemak Bay
# Data analysis work flow and figure generation

# ============================================================
# Diversity Metrics per Site
# - Richness (S), Shannon (H), Simpson (1-D), evenness (J)
# - X axis = Year
# - Low Stratum Only
# - Separate panel for each index
# - Fill = Site
# Data file: KBAY2012-2022_Rocky_Intertidal_ Motile_Invert_Count.csv
# ============================================================
library(lubridate)
library(tidyverse)
library(vegan)

# Loading in Community Data ----------------------------------------------------
RI_Motile_Invert <- read.csv("KBAY2012-2022_Rocky_Intertidal_ Motile_Invert_Count.csv")

# DIVERSITY METRICS ------------------------------------------------------------
# Calculating the alpha diversity indices
# - Richness (S), Shannon (H), Simpson (1-D), evenness (J)
Diversity_Metrics <- RI_Motile_Invert |>
  group_by(Site, Year, Stratum, ScientificName_accepted) |>
  summarise(Abundance = sum(X.counts, na.rm = TRUE), .groups = "drop") |>
  group_by(Site, Year, Stratum) |>
  summarise(
    Species_richness = n_distinct(ScientificName_accepted),
    Total_individuals = sum(Abundance),
    shannon_H = diversity(Abundance, index = "shannon"),
    simpson_1minusD = diversity(Abundance, index = "simpson"),
    q0 = Species_richness, 
    q1 = exp(diversity(Abundance, index = "shannon")), 
    q2 = diversity(Abundance, index = "invsimpson"), 
    pielou_evenness = ifelse(
      Species_richness > 1,
      shannon_H / log(Species_richness),
      NA_real_
    ),
    .groups = "drop"
  )

# Make sure Year is ordered correctly
Diversity_Metrics$Year <- as.numeric(as.character(Diversity_Metrics$Year))
# (use this if Year is currently a factor)

# Filtering for only the low intertidal stratum 
Low_Diversity <- Diversity_Metrics |>
  filter(Stratum == "Low")

# Alpha diversity indices plots ------------------------------------------------
S <- ggplot(Low_Diversity,
            aes(x = Year,
                y = Species_richness,
                fill = Site)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7) +
  labs(x = "Year",
       y = "Species Richness (S)",
       fill = "Site") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(face = "bold")
  )

S

ggsave(filename = "Species_Richness_Grouped.svg")

H <- ggplot(Low_Diversity,
            aes(x = Year,
                y = shannon_H,
                fill = Site)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7) +
  labs(x = "Year",
       y = "Shannon Index (H)",
       fill = "Site") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(face = "bold")
  )

H

ggsave(filename = "Shannon_Index_Grouped.svg")

D <- ggplot(Low_Diversity,
            aes(x = Year,
                y = simpson_1minusD,
                fill = Site)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7) +
  labs(x = "Year",
       y = "Simpsons Index (1-D)",
       fill = "Site") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(face = "bold")
  )

D

ggsave(filename = "Simpsons_Index_Grouped.svg")

J <- ggplot(Low_Diversity,
            aes(x = Year,
                y = pielou_evenness,
                fill = Site)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7) +
  labs(x = "Year",
       y = "Pielou Eveness (J)",
       fill = "Site") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(face = "bold")
  )

J

ggsave(filename = "Pielou_Eveness_Grouped.svg")

# ============================================================
# Community Composition Analysis
# -	NMDS (Bray-Curtis) for Shanon Diversity
# Hypothesis test
# -	ANOSIM: tested the differences between communities
# -	SIMPER: which taxa contributed most to the observed differences between sites
# Data file: KBAY2012-2022_Rocky_Intertidal_ Motile_Invert_Count.csv
# ============================================================

library(dplyr)
library(tidyr)

# Visualizing Bray-Curtis dissimilarity with NMDS
# ---- 1) Read data ----
RI_Motile_Invert <- RI_Motile_Invert |>
  mutate(Year = as.factor(Year), Site = as.factor(Site))

# Filtering for only the low stratum (intertidal area)
dat <- RI_Motile_Invert |>
  filter(Stratum == "Low")

# Fixing rows with repeated species per Site-Year-Species
dat_summarised <- dat |>
  group_by(Site, Year, ScientificName_accepted) |>
  summarise(
    X.counts = sum(X.counts, na.rm = TRUE),
    .groups = "drop"
  )

dat_wide <- dat_summarised |>
  pivot_wider(
    id_cols = c(Site, Year),
    names_from = ScientificName_accepted,
    values_from = X.counts,
    values_fill = 0
  )

# Metadata columns 
meta_cols <- c("Site", "Year")
meta <- dat_wide[, meta_cols]
comm <- dat_wide[, setdiff(names(dat_wide), meta_cols)]

# ---- 2) Run Bray-Curtis NMDS ----
# Bray-Curtis is common for abundance data; try Jaccard for presence/absence
set.seed(123)
nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)

cat("NMDS stress:", nmds$stress, "\n")

# ---- 3) Extract scores ----
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$Site <- meta$Site
site_scores$Year <- meta$Year

sp_scores <- as.data.frame(scores(nmds, display = "species"))
sp_scores$ScientificName_accepted <- rownames(sp_scores)

write.csv(site_scores, "nmds_scores_sites.csv", row.names = FALSE)
write.csv(sp_scores, "nmds_scores_species.csv", row.names = FALSE)

# ---- 4 ) Plot NMDS (ggplot2)
# NMDS plotting  ----------------------------------------------------

# finding convex hulls per site
site_hulls <- site_scores |>
  group_by(Site) |>
  slice(chull(NMDS1, NMDS2))

s <- ggplot(site_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_polygon(data = site_hulls,
               aes(fill = Site, group = Site),
               alpha = 0.2,
               color = NA) +
  geom_point(aes(shape = Site, color = Year),
             size = 4,        # larger points for 18 cm width
             stroke = 1) +    # ensures outline is not too thin
  theme_bw(base_family = "Times New Roman") +
  labs(title = "NMDS (Bray-Curtis) — community composition",
       subtitle = paste0("Stress = ", round(nmds$stress, 3)),
       x = "NMDS1",
       y = "NMDS2") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_line(linewidth = 1)
  )

s

ggsave(
  filename = "NMDS_Site_Year_SitePolygons.svg",
  plot = s,
  width = 18,
  height = 14,      # keeps good aspect ratio
  units = "cm"
)

# ---- 6) Diagnostics: Shepard plot (stress plot) ----
png("nmds_shepard.png", width = 900, height = 700, res = 130)
stressplot(nmds)
dev.off()

# ---- 7) ANOSIM Hypothesis Testing
# ---- Distance matrix ----
d <- vegdist(comm, method = "bray")

# ---- Grouping variable ----
# Test differences among Sites
group_site <- meta$Site

# ---- Run ANOSIM - Site ----
set.seed(42)
ano_site <- anosim(d, grouping = group_site, permutations = 999)

print(ano_site)

# ---- 8) SIMPER Analysis
simp_site <- simper(comm, group_site, permutations = 999)
summary(simp_site)

# Summarize SIMPER by group-pair and export
pair_names <- names(simp_site)

# Build a compact CSV of top contributors per pair (top 10)
out_rows <- list()
txt_lines <- c("SIMPER top taxa per group pair (top 10)\n=======================================\n")

for (pn in pair_names) {
  s <- simp_site[[pn]]
  
  s_df <- as.data.frame(s)
  s_df <- as.data.frame(s)
  
  # If SIMPER stored species as a column
  if ("species" %in% colnames(s_df)) {
    s_df$taxon <- s_df$species
  }
  
  s_df <- s_df[order(-s_df$average), ]
  
  topN <- head(s_df, 10)
  topN$pair <- pn
  out_rows[[pn]] <- topN
  
  txt_lines <- c(txt_lines, paste0("\nPair: ", pn))
  txt_lines <- c(txt_lines, paste0("Top taxa:"))
  
  for (i in seq_len(nrow(topN))) {
    txt_lines <- c(txt_lines,
                   sprintf("  %2d) %s | avg=%.4f | cumsum=%.3f",
                           i,
                           as.character(topN$taxon[i]),
                           topN$average[i],
                           topN$cusum[i]))
  }
}

simper_top <- do.call(rbind, out_rows)
write.csv(simper_top, "simper_summary_by_pair.csv", row.names = FALSE)
writeLines(txt_lines, "simper_top_taxa_per_pair.txt")

# ============================================================
# Temperature effects on Shannon Diversity
# -	GLMM with site as a random effect
# Data file: KBAY2012-2022_Rocky_Intertidal_ Motile_Invert_Count.csv
# Data file: KBAY2012-2022_Intertidal_Temp.csv
# ============================================================
# TEMPERATURE EFFECTS ON SPECIES RICHNESS --------------------------------------
# GLM with temperature Shannon Diversity - with linear regression and trendline

# Loading in Community Data ----------------------------------------------------
RI_Motile_Invert <- read.csv("KBAY2012-2022_Rocky_Intertidal_ Motile_Invert_Count.csv")
# Convert species dates correctly
RI_Motile_Invert$Date <- mdy(RI_Motile_Invert$Date)
# filtering for only the low stratum
species <- RI_Motile_Invert |>
  filter(Stratum == "Low")

# Loading in Temperature Data --------------------------------------------------
combined <- read.csv("KBAY2012-2022_Intertidal_Temp.csv")
#filtering for Temperature data only when the HOBO logger was submerged in water and not exposed to air
combined_IT <- combined |>
  mutate(exposure = as.factor(exposure)) |>
  filter(exposure == "water")
temp <- combined_IT

# Making site and date columns match
temp <- temp |> 
  rename(Site = site, Date = date)
# Making site names match
temp <- temp |>
  mutate(
    Site = recode(Site,
                  "Bishop_Bch"  = "Bishop's Beach",
                  "Bluff_Pt"    = "Bluff Point",
                  "Cohen_Isl"   = "Cohen Island",
                  "Outside_Bch" = "Outside Beach",
                  "Port_Graham" = "Port Graham", 
                  "Elephant_Is" = "Elephant Island") )

# Ensure Date columns are proper Date class
species$Date <- as.Date(species$Date)
temp$Date    <- as.Date(temp$Date)

library(data.table)

# Convert to data.tables
species <- as.data.table(species)
temp <- as.data.table(temp)

# Set keys for efficient join
setkey(species, Site, Date)
setkey(temp, Site, Date)

# Finding temperature data for the closest day to invertebrate counts and sampling
closest_temp <- temp[species, roll = "nearest"]
# Finding an average temp value for the sampling date for each site-year
annual_temp <- closest_temp[,.(annual_mean_temp = mean(temperature, na.rm = TRUE)),
                            by = .(Site, Year)]

# DIVERSITY METRICS ------------------------------------------------------------
# Calculating the alpha diversity indices
# - Richness (S), Shannon (H), Simpson (1-D), evenness (J)
Diversity_Metrics <- RI_Motile_Invert |>
  group_by(Site, Year, Stratum, ScientificName_accepted) |>
  summarise(Abundance = sum(X.counts, na.rm = TRUE), .groups = "drop") |>
  group_by(Site, Year, Stratum) |>
  summarise(
    Species_richness = n_distinct(ScientificName_accepted),
    Total_individuals = sum(Abundance),
    shannon_H = diversity(Abundance, index = "shannon"),
    simpson_1minusD = diversity(Abundance, index = "simpson"),
    q0 = Species_richness, 
    q1 = exp(diversity(Abundance, index = "shannon")), 
    q2 = diversity(Abundance, index = "invsimpson"), 
    pielou_evenness = ifelse(
      Species_richness > 1,
      shannon_H / log(Species_richness),
      NA_real_
    ),
    .groups = "drop"
  )

# Make sure Year is ordered correctly
Diversity_Metrics$Year <- as.numeric(as.character(Diversity_Metrics$Year))
# (use this if Year is currently a factor)

# Filtering for only the low intertidal stratum 
Low_Diversity <- Diversity_Metrics |>
  filter(Stratum == "Low")

merged_data <- merge(
  Low_Diversity,
  annual_temp,
  by = c("Site", "Year"),
  all.x = TRUE
)

# Linear Regression ---------------------------------------------------------
library(lme4)

H_mixed_regression <- lmer(shannon_H ~ annual_mean_temp + (1 | Site), data = merged_data)
summary(H_mixed_regression)

# Plotting Linear Regression----------
t <- merged_data |>
  ggplot(aes(x = annual_mean_temp, y = shannon_H)) +
  geom_point(aes(color = Site), size = 3) +
  geom_smooth(color = "black", method = "lm", se = TRUE) +
  xlab("Temperature (°C)") +
  ylab("Shannon Index (H)") +
  theme_classic(base_size = 16) +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(hjust = 1),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
  )

ggsave(
  filename = "GLM_Shannon_Temp.svg",
  plot = t,
  width = 18,
  height = 14,      # keeps good aspect ratio
  units = "cm"
)

# Visualizing temp data
st <- ggplot(merged_data, aes(x = Year, y = annual_mean_temp)) +
  geom_line(aes(group = 1), linewidth = 1) +  # group=1 needed for lines with numeric/discrete x
  geom_point(size = 3) +
  facet_wrap(~Site, ncol = 2) +
  theme_classic(base_size = 14) +
  scale_x_continuous(breaks = seq(min(merged_data$Year), max(merged_data$Year), by = 1))+
  labs(
    x = "Year",
    y = "Sampling Day Water Temperature (°C)"
  ) +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 16)
  )
st

residuals(H_temp_regression)

# Next, let's visualize the residuals
plot(residuals(H_temp_regression) ~ annual_mean_temp, data = merged_data,
     col = "blue", las = 1, pch = 16, ylab = "Shannon Index Model Residuals", xlab=
       "Temperature (°C)") + abline(h=0)

# ============================================================
# Stacked relative abundance bar plot of Family composition
# - Overall relative abundance per Site
# - X axis = Year
# - Low Stratum Only
# - Separate panel for each Site
# - Fill = Species
# Data file: KBAY2012-2022_Rocky_Intertidal_ Motile_Invert_Count.csv
# ============================================================

# ----------------------------
# 0) Install + load packages
# ----------------------------
pkgs <- c("tidyverse", "scales")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)

library(tidyverse)
library(scales)

# ----------------------------
# 1) Read data
# ----------------------------
infile <- "KBAY2012-2022_Rocky_Intertidal_ Motile_Invert_Count.csv"
dat <- read_csv(infile, show_col_types = FALSE)
dat |>
  filter(Stratum == "Low")

# ----------------------------
# 2) Check required columns
# ----------------------------
req_cols <- c("Year", "Site", "ScientificName_accepted", "#counts")
missing_cols <- setdiff(req_cols, names(dat))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ----------------------------
# 3) Clean + summarize data
#    Overall counts per Site-Year-Species
# ----------------------------
plot_dat <- dat %>%
  mutate(
    Year   = as.factor(Year),
    Site   = as.factor(Site),
    Species = as.factor(ScientificName_accepted),
    counts = suppressWarnings(as.numeric(`#counts`))
  ) %>%
  filter(
    !is.na(Year),
    !is.na(Site),
    !is.na(Species),
    !is.na(counts)
  ) %>%
  group_by(Site, Year, Species) %>%
  summarise(total_count = sum(counts, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site, Year) %>%
  mutate(rel_abundance = total_count / sum(total_count, na.rm = TRUE)) %>%
  ungroup()

# ----------------------------
# 4) Optional: lump rare species into "Other"
#    Keep top N most abundant species overall
# ----------------------------
top_n_species <- 12

top_species <- plot_dat %>%
  group_by(Species) %>%
  summarise(grand_total = sum(total_count, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(grand_total)) %>%
  slice_head(n = top_n_species) %>%
  pull(Species)

plot_dat2 <- plot_dat %>%
  mutate(Species = if_else(Species %in% top_species, as.character(Species), "Other")) %>%
  group_by(Site, Year, Species) %>%
  summarise(total_count = sum(total_count, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site, Year) %>%
  mutate(rel_abundance = total_count / sum(total_count, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Species = as.factor(Species))

# ----------------------------
# 5) Plot stacked relative abundance
# ----------------------------
p <- ggplot(plot_dat2, aes(x = Year, y = rel_abundance, fill = Species)) +
  geom_col(width = 0.8, color = "grey20", linewidth = 0.15) +
  facet_wrap(~ Site) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Year",
    y = "Relative abundance",
    fill = "Species",
    title = "Decadal Relative Species Abundance",
    subtitle = "Low Intertidal Overall Composition per Site"
  ) +
  theme_classic(22) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
    legend.position = "right",
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )

print(p)

# ----------------------------
# 6) Save figure
# ----------------------------
ggsave(
  filename = "Species_relative_abundance_by_site_year.svg",
  plot = p,
  width = 14,
  height = 8,
  dpi = 300
)

# ----------------------------
# 7) Save summarized plotting data
# ----------------------------
write_csv(plot_dat2, "Species_relative_abundance_by_site_year_data.csv")

