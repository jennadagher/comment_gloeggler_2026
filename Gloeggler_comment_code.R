##### Comment on: Turkana Miocene to modern East Africa traits against environmental proxies
# Jenna Dagher 2026/19/07
####

# setwd("/Users/jennadagher/Desktop/Gloggler response/Data & Code")

pacman::p_load("tidyverse", "readxl", "palaeoverse", "yacca")

# destination_path <- "/Users/jennadagher/Desktop/Gloggler response/Data & Code"

# default_crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# theme_set(ggpubr::theme_pubclean()+theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                                         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')))

# color_scheme <- c("East" = "#1C9E77",
#                  "North" = "#7570B3",
#                  "West" = "#D95F02",
#                  "Tugen Hills"="#E72A8A",
#                  "Ileret" = "#6A3D9A",
#                  "Karari Ridge" = "#1F78B4",
#                  "Koobi Fora Ridge"="#33A02B",
#                  "Turkana - general" = "#00FFFF"
)

all_traits <- c("logBM", "HYP", "LOP", "AL", "ALX", "OL", "OLX", "GT", "BUN", "SF", "OT")

# genus occurrence list

genus_list <- read_csv("genus_number.csv") |>
  filter(
    is.na(Delete),
    Genus != "Homo" | is.na(Genus)
  ) |>
  mutate(Member = case_when(
    str_detect(Member, "Lomekwi") & !Member %in% c("Lower Lomekwi", "Upper Lomekwi") & min_ma > 3 ~ "Lower Lomekwi",
    str_detect(Member, "Lomekwi") & !Member %in% c("Lower Lomekwi", "Upper Lomekwi") & max_ma < 3 ~ "Upper Lomekwi",
    TRUE ~ Member
  ))


genus_traits <- read_xlsx("genus_traits.xlsx") |>
  mutate(across(all_of(c(all_traits, "Body_Mass_Kg")), as.numeric))

# Modifications to original code - Clade-only null model

genus_traits <- read_xlsx("genus_traits.xlsx") |>
  mutate(across(all_of(c(all_traits, "Body_Mass_Kg")), as.numeric))

genus_taxonomy <- read_csv("genus_number.csv") |>
  distinct(Genus, Family, Subfamily, Tribe)

genus_traits <- genus_traits |>
  left_join(genus_taxonomy, by = "Genus") |>
  mutate(clade = coalesce(Tribe, Subfamily, Family)) |>
  group_by(clade) |>
  mutate(
    across(
      all_of(c(all_traits, "Body_Mass_Kg")),
      ~ {
        if (all(is.na(.x))) {
          .x
        } else {
          mean(.x, na.rm = TRUE)
        }
      }
    )
  ) |>
  ungroup() |>
  select(-clade, -Family, -Subfamily, -Tribe)

# sanity check, printed immediately
genus_traits |> summarise(across(all_of(c(all_traits, "Body_Mass_Kg")), ~ sum(is.na(.))))


# Create time bins ---------------------------------------------------------

occdf <- genus_list |>
  select(Locality, Place, Member, Formation, max_ma, min_ma) |>
  filter(!is.na(min_ma) & !is.na(max_ma)) |>
  mutate(max_ma = if_else(max_ma - min_ma == 0, max_ma + 1e-7, max_ma)) # add 1/10 years in case of point dates, so that bin_time works with method = "majority" 

bins <- data.frame(bin = 1:10, min_ma = 0:9, max_ma = 1:10) |>
  rbind(c(0, 0, 0))

binned_Mb <- bin_time(occdf = occdf, bins = bins, method = "majority") |> #  "majority" gives percentage of overlap, can be used to filter out overlaps < 30-40%
  mutate(
    bin_assignment = case_when(Member == "modern" ~ 0, TRUE ~ as.numeric(bin_assignment)),
    bin_midpoint   = case_when(Member == "modern" ~ 0, TRUE ~ as.numeric(bin_midpoint))
  ) |>
  cbind(occdf$Locality) |>
  filter(overlap_percentage > 40)

# --- combine genus occurrences + traits, aggregate to locality level ---

min_body_mass_kg <- 1
min_genus_number_per_loc <- 5

loc_mb_mean_scores <- left_join(genus_list, genus_traits, by = join_by(Genus)) |>
  filter(Body_Mass_Kg >= min_body_mass_kg | is.na(Genus)) |>
  mutate(grouping_variable = coalesce(Tribe, Subfamily, Family)) |>
  group_by(grouping_variable) |>
  mutate(across(all_of(all_traits), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) |>
  group_by(Locality, Place, Member) |>
  mutate(
    across(c(all_of(all_traits), "Body_Mass_Kg"), function(.x) mean(.x, na.rm = TRUE), .names = "{col}_mean"),
    across(c(all_of(all_traits), "Body_Mass_Kg"), function(.x) sd(.x, na.rm = TRUE),   .names = "{col}_SD")
  ) |>
  group_by(Locality, Place, Member) |>
  mutate(n_genus = n()) |>
  ungroup() |>
  filter(n_genus >= min_genus_number_per_loc | Member == "modern") |>
  select(Region, Subregion, Locality, Place, Member, Formation, ends_with("_mean"), ends_with("_SD"),
         n_genus, min_ma, max_ma) |>
  distinct(Locality, Place, Member, .keep_all = TRUE) |>
  inner_join(binned_Mb[c("Locality", "Place", "Formation", "Member", "n_bins",
                         "bin_assignment", "bin_midpoint", "max_ma", "min_ma")],
             by = join_by(Locality, Place, Member, Formation),
             multiple = "first", suffix = c("", "_bin")) |>
  rename_with(.fn = \(x) sub("_mean$", "", x)) |>
  filter(!is.na(min_ma)) |>
  mutate(mean_ma = (max_ma + min_ma) / 2,
         Member = str_replace_all(Member, c("TuluBor" = "Tulu Bor",
                                            "UpperBurgi" = "Upper Burgi",
                                            "lower Kalochoro" = "Kalochoro")))

# Removed plotting code not needed for comment analysis

# read in environmental proxies -----------

pedogenic_isotopes <- read_xlsx("pedogenic_carbonate_isotopes.xlsx") |>
  mutate(
    Member = case_when(
      Formation == "Tugen Hills - modern" ~ "modern",
      Region == "Tugen Hills" ~ Formation,
      str_detect(Member, "Lomekwi") & !Member %in% c("Lower Lomekwi", "Upper Lomekwi") & mean_ma > 3 ~ "Lower Lomekwi",
      str_detect(Member, "Lomekwi") & !Member %in% c("Lower Lomekwi", "Upper Lomekwi") & mean_ma < 3 ~ "Upper Lomekwi",
      TRUE ~ Member
    ),
    Member = str_replace(Member, "middle", "Middle"),
    Formation = if_else(Formation == "Tugen Hills - modern", "modern", Formation),
    bin_assignment = ifelse(mean_ma == 0, 0, cut(mean_ma, breaks = seq(0, 10), labels = seq(1, 10)))
  ) |>
  group_by(Region, Member) |>
  filter(n() >= 5) |>
  ungroup()

region_mb_mean_CO_isotope <- pedogenic_isotopes |>
  mutate(
    mean_d13C = mean(d13C, na.rm = TRUE), SD_d13C = sd(d13C, na.rm = TRUE),
    mean_d18O = mean(d18O, na.rm = TRUE), SD_d18O = sd(d18O, na.rm = TRUE),
    .by = c(Region, Member)
  ) |>
  select(Region, Member, Formation, mean_d13C, SD_d13C, mean_d18O, SD_d18O) |>
  distinct(Region, Member, .keep_all = TRUE)

# Removed modern analysis and plotting code not necessary for comment

# Canonical correlation analysis (CCorA) with null model and compare with genus-level trait analysis ------

# Genus-level trait CCorA
ccora_dataframe_real <- region_mb_combined_CO |>                 # region_mb_combined_CO and region_mb_combined_enamel
  filter(Region != "Tugen Hills")
trait_metrics <- "mean"
isotope_metrics <- "mean"

traits_real <- ccora_dataframe_real[paste0(c("mean" = "mean_", "SD" = "SD_")[trait_metrics], all_traits)]
isotopes <- ccora_dataframe_real[paste0(c("mean" = "mean_", "SD" = "SD_")[isotope_metrics], c("d13C", "d18O"))]

cc <- yacca::cca(traits_real, isotopes, na.rm = TRUE)


CV1_real  <- as.matrix(traits_real)  %*% cc$xcoef[, 1]


# Clade-only null model CCorA
region_mb_mean_trait <- loc_mb_mean_scores |>
  group_by(Member) |>
  filter(max(max_ma) - min(min_ma) < 2) |>
  ungroup() |>
  mutate(
    across(c(all_of(all_traits), "Body_Mass_Kg"), function(.x) mean(.x, na.rm = TRUE), .names = "mean_{col}"),
    across(c(all_of(all_traits), "Body_Mass_Kg"), function(.x) sd(.x, na.rm = TRUE),   .names = "SD_{col}"),
    .by = c(Region, Member)
  ) |>
  select(Region, Subregion, Member, Formation, starts_with("mean_"), starts_with("SD_"),
         mean_ma, min_ma, max_ma, bin_assignment) |>
  distinct(Member, Region, .keep_all = TRUE)

region_mb_combined_CO <- inner_join(region_mb_mean_trait, region_mb_mean_CO_isotope,
                                    by = join_by(Region, Member, Formation),
                                    na_matches = "never")

ccora_dataframe_clade <- region_mb_combined_CO |>
  filter(Region != "Tugen Hills")

traits_clade   <- ccora_dataframe_clade[paste0("mean_", all_traits)]
isotopes <- ccora_dataframe_clade[paste0("mean_", c("d13C", "d18O"))]

cc_cl <- yacca::cca(traits_clade, isotopes, na.rm = TRUE)


CV1_clade  <- as.matrix(traits_clade)  %*% cc_cl$xcoef[, 1]

cc_clade <- cc_cl$corr[1]   

cc_clade

CCP::p.asym(cc$corr, dim(traits_real)[1], length(traits_real), length(isotopes), tstat = "Wilks") 
CCP::p.asym(cc_cl$corr, dim(traits_clade)[1], length(traits_clade), length(isotopes), tstat = "Wilks") 

cor(CV1_real, CV1_clade)



# PRODUCING FIG 1 -----------------------------------------------------------

library(tidyverse)
library(patchwork)

# Panel A: CV1_real vs CV1_clade scatter (same-axis evidence)

stopifnot(identical(ccora_dataframe_real$Region, ccora_dataframe_clade$Region),
          identical(ccora_dataframe_real$Member, ccora_dataframe_clade$Member))

cv1_df <- tibble(
  Region    = ccora_dataframe_real$Region,
  Member    = ccora_dataframe_real$Member,
  CV1_real  = as.numeric(CV1_real),
  CV1_clade = as.numeric(CV1_clade)
)

r_cv1 <- cor(cv1_df$CV1_real, cv1_df$CV1_clade)

panel_A <- ggplot(cv1_df, aes(CV1_real, CV1_clade, color = Region)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_color_manual(values = color_scheme) +   # reuse their published color scheme
  annotate("text", x = min(cv1_df$CV1_real), y = max(cv1_df$CV1_clade),
           label = paste0("r = ", round(r_cv1, 2)),
           hjust = 0, vjust = 1, size = 5, fontface = "bold") +
  labs(x = "CV1 (published genus-level traits)",
       y = "CV1 (clade-only traits)",
       title = "A") +
  ggpubr::theme_pubclean() +
  theme(legend.position = "right")


cancor_df <- tibble(
  model  = factor(c("Published\n(genus-level)", "Clade-only\n(no within-clade variation)"),
                  levels = c("Published\n(genus-level)", "Clade-only\n(no within-clade variation)")),
  cancor = c(cc$corr[1], cc_cl$corr[1])
)

pct_retained <- round(100 * cancor_df$cancor[2] / cancor_df$cancor[1], 0)

panel_B <- ggplot(cancor_df, aes(model, cancor, fill = model)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = round(cancor, 3)), vjust = -0.5, size = 5) +
  annotate("text", x = 1.5, y = max(cancor_df$cancor) * 1.12,
           label = paste0(pct_retained, "% retained"),
           size = 5, fontface = "italic") +
  scale_fill_manual(values = c("#1C9E77", "#D95F02")) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(x = NULL, y = "Canonical correlation (CV1)", title = "B") +
  ggpubr::theme_pubclean()


# Combine into one figure
final_figure <- panel_A + panel_B + plot_layout(widths = c(1.3, 1))
final_figure


# ggsave("/Users/jennadagher/Desktop/Gloggler response/figure1.png", final_figure, width = 11, height = 5)



