library(tidyverse)
library(cowplot)

data <- read_csv("00_data/all_gp_site_info.csv")

# Only plots that had more than 1 measurement over time.
remes_ids <- data %>% filter(meas_num %in% 1:2) %>% select(nfi_plot) %>% as.matrix
data_filt <- data %>% filter(nfi_plot %in% remes_ids)
data_filt$loc_id %>% table # loc_id 0 are sites with no location change. So, no sites with location change


# Removing volume columns
data2 <- data %>% 
  select(nfi_plot:land_pos, veg_type:ec_district, grossvol_mai,
         contains("plotbio"), contains("cc"), contains("sample_depth"), contains("avg_")) %>% 
  filter(rowSums(is.na(across(grossvol_mai:ncol(.)))) != length(across(grossvol_mai:ncol(.))))

# Check
setdiff(names(data), names(data2))
# aggregate variables following IPCC Good Practice Guidance and CBM3 standards
# CBM-CFS3: A model of carbon-dynamics in forestry and land-use change implementing 
# IPCC standards - doi:10.1016/j.ecolmodel.2008.10.018
var_class <- tibble(VARIABLE = names(data2),
                    VARIABLE_GPG = c(rep("ID", 4),
                                     rep("META", 12),
                                     rep("GROUPING", 8),
                                     "NOTES",
                                     rep("GROUPING", 2),
                                     "VOL",
                                     rep("DW", 2),# Stump and small stumps
                                     "AGB", "DW", # Tree
                                     "AGB", "DW", # Small Tree
                                     "AGB", "DW", # Large Shrub
                                     "AGB", "AGB", "AGB", # Small shrub to bryo
                                     "Litter", "Litter", # Fine and small wood debris
                                     "DW", "DW", "DW", # coarse, round and odd wood debris
                                     "Litter", "Litter",
                                     rep("Soil TBA", 13)
                    ),
                    VARIABLE_CBM = c(rep("ID", 4),
                                     rep("META", 12),
                                     rep("GROUPING", 8),
                                     "NOTES",
                                     rep("GROUPING", 2),
                                     "VOL",
                                     "Other Wood", "Snag branches DOM", 
                                     "Merch + Bark + Foliage", "Snag Stem + Snag Branches DOM",
                                     "Other Wood", "Snag branches DOM",
                                     NA, NA,
                                     NA, NA, NA,
                                     "Aboveground very fast DOM", "Aboveground fast DOM",
                                     "Aboveground fast DOM", "Aboveground fast DOM", "Aboveground fast DOM",
                                     "Aboveground very fast DOM", "Aboveground very fast DOM",
                                     rep("Belowground TBA", 13)
                    )
)

# and check completeness
# Individual sites

check_miss <- data2 %>% 
  filter(loc_id == 0) %>% # Excluding sites that changed locations.
  select(nfi_plot, meas_num, plotbio_stump:ncol(.)) %>% 
  pivot_longer(names_to = "VARIABLE", values_to = "VALUE", plotbio_stump:ncol(.)) %>% 
  pivot_wider(names_from = "meas_num", values_from = "VALUE") %>% 
  mutate(N_MEAS = rowSums(!is.na(across(`0`:`2`)))) %>% 
  left_join(var_class)

var_order <- check_miss %>% group_by(VARIABLE) %>% summarise(N_MEAS = sum(N_MEAS)) %>% arrange(-N_MEAS) %>% .$VARIABLE
nfi_order <- check_miss %>% group_by(nfi_plot) %>% summarise(N_MEAS = sum(N_MEAS)) %>% arrange(-N_MEAS) %>% .$nfi_plot

check_miss <- check_miss %>%
  mutate(VARIABLE = factor(VARIABLE, levels = var_order),
         nfi_plot = factor(nfi_plot, levels = nfi_order))

ggplot(check_miss) +
  geom_tile(aes(VARIABLE, factor(nfi_plot), fill = factor(N_MEAS))) +
  scale_fill_manual(values = c("red4", "salmon3", "palegreen3", "turquoise4")) +
  theme_half_open() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0)) + 
  labs(x = "Variable", y = "NFI plot", fill = "Number of \nmeasurements") +
  coord_flip()

smry_miss <- check_miss %>% group_by(VARIABLE, VARIABLE_CBM, VARIABLE_GPG) %>% count(N_MEAS) %>% 
  mutate(N_MEAS = factor(N_MEAS, levels = 0:3),
         VARIABLE = factor(VARIABLE, levels = rev(var_order)))

ggplot(smry_miss) +
  geom_col(aes(VARIABLE, n, fill = N_MEAS)) +
  scale_fill_manual(values = c("red4", "salmon3", "palegreen3", "turquoise4")) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 650, by = 50)) +
  facet_grid(VARIABLE_CBM~., scales = "free", space = "free", labeller = as_labeller(label_wrap_gen(14))) +
  coord_flip() +
  theme_half_open() +
  theme(axis.text.x = element_text(angle = 00, vjust = 0.5, hjust = 0.5),
        strip.text.y = element_text(angle = 0)) +
  labs(x = "Biomass variable calculated on NFI", y = "NFI plot", fill = "Number of \nmeasurements")

# Which variables to keep? and which plots to keep.