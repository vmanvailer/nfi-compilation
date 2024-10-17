library(tidyverse)

data <- read_csv("00_data/all_gp_site_info.csv")

# Starting exploration
# First concern is identifying variables of interest. All biomass estimates plus
# carbon content on forest floor.

# Sasha mentioned that carbon from soil samples may vary solely because for each
# campaign effort, one of four microplots is measure at each campaign and 
# therefore differences are largely varying because of spatial differences {CHECK THAT INFO BECAUSE I REMEBER THAT ALL OF THEM WERE MEASURED}
 
# Starting with basic summary statistics about the data
data$nfi_plot %>% unique %>% length #1060 plots

# There are four important variables relating to each plot: 
# loc_id            = if a given plot changed locations i.e. due to new danger to sampling or other reason
# meas_date         = when the sampling happened
# meas_num          = how many time a plot has been measured
# plot_completeness = if all parameters were measured or not 

# loc_id ---------------------------------------------------------------

# Number of sites that changed locations
data$loc_id %>% table
# Which sites were those 
loc_change_nfi_id <- data %>% filter(loc_id == 1) %>% select(nfi_plot) %>% as.matrix()
# Confirming they all have two measurements
data[data$nfi_plot %in% loc_change_nfi_id, "nfi_plot"] %>% table %>% View
# One of the sites only has measurements for the new site and not the original one 
data[data$nfi_plot == 1074436,]

data[data$nfi_plot %in% loc_change_nfi_id, ] %>% View


# meas_num --------------------------------------------------------------

# Check the number plots that did not have location changes
remes <- data[!data$nfi_plot %in% loc_change_nfi_id, "meas_num"] %>% table %>% as_tibble
# Adding the number of plots that had a change
remes <- rbind(remes, as_tibble(table(data[!data$loc_id != 1, "meas_num"]))) %>% mutate(LOC_CHANGED = c("No_change", "No_change", "No_change", "Changed"))

ggplot(remes) +
  geom_col(aes(meas_num, n, fill = interaction(meas_num, LOC_CHANGED))) +
  scale_fill_manual(values = c("#868686", "#5C5C5C", "#A29745", "#A24535")) +
  scale_y_continuous(breaks = seq(0,1100, by = 100)) +
  cowplot::theme_half_open() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey90"),
        legend.position = "none") +
  labs(x = "Measurement number", y = "Number of plots",
       fill = "Location",
       # title = "Measurements on sites that did not changed locations"
  )

# plot_completeness -----------------------------------------------------

data$plot_completion %>% table
# There a lot of plots for which the completeness was not recorded (N)
# Also Sasha mentioned that there shouldn't be more than 1 re-measurement. So maybe some are incomplete and are remeasured to be complete?

# Check N of plots per completion status per measurement
# Count sites that didn't change location separetly for plotting

# Sites that didn't change locations 
compl <- data[!data$nfi_plot %in% loc_change_nfi_id, c("meas_num", "plot_completion")] %>% table %>% as_tibble()

# Sites that didn't change location 
compl <- rbind(compl, as_tibble(table(data[!data$loc_id != 1, c("meas_num", "plot_completion")]))) %>% 
  mutate(LOC_COMP = paste(c(rep("No_change", 12), rep("Changed", 2)), plot_completion, sep = "_"))

# Adjust levels for plotting
compl$LOC_COMP <- factor(compl$LOC_COMP,
                         levels = c("Changed_F", "Changed_P",
                                    "No_change_F", "No_change_P", "No_change_N", "No_change_U"))

ggplot(compl) +
  geom_col(aes(meas_num, n, fill = LOC_COMP), linewidth = 1.5) +
  scale_fill_manual(values = c("palegreen3", "indianred3", "forestgreen",
                               "red4", "grey30", "black"),
                    labels = c("Full (location change)", "Partial (location change)",
                               "Full", "Partial", "Not reported", "Incomplete")) +
  scale_y_continuous(breaks = seq(0, 1050, by = 50)) +
  theme_minimal() +
  labs(x = "Measurement number",
       y = "Number of plots",
       title = "Completeness status",
       fill = "Measurement completeness")

# Make a table for quick access to numbers.
pcomp <- compl %>% 
  mutate(LOC_COMP = c(rep("No_change", 12), rep("Changed", 2))) %>%
  pivot_wider(names_from = "plot_completion", values_from = "n")

write_csv(pcomp, "NFI_S03_tally_of_plot_completion_codes.csv" )

# Just double check plots that did change location to include those numbers in the table
data[data$nfi_plot %in% loc_change_nfi_id, c("plot_completion", "meas_num")] %>% table

# Exploring NA structure to understand completeness status ------------------
library(naniar)

vis_miss(filter(data,
                plot_completion == "N",
                # nfi_plot %in% loc_change_nfi_id
                ), cluster = F) + 
  theme(axis.text.x = element_text(angle = 90))

# Some plots are missing a lot of variables. Which plots are missing all variables from small tree to avg bulk density
idx_many_miss <- apply(select(data, plotbio_smtr_live:avg_bulk_density_org), 1, function(x) all(is.na(x))) %>% which()
nfi_m_miss <- data[idx_many_miss, "nfi_plot"]

vis_miss(filter(data, nfi_plot %in% nfi_m_miss$nfi_plot)) + 
  theme(axis.text.x = element_text(angle = 90))
# They do have more than one measurement and seem to be mostly in quebec

# Measurement span by measurement number ---------------------

data2 <- separate(data, meas_date, into = c("YEAR", "MONTH", "DAY"), sep = "-", remove = FALSE)

# Selective Availability was lifted in May 1st 2000
barplot(table(select(filter(data2, YEAR <= 2000), MONTH))) # No sites measure before may

meas_dumb <- data2[,c("nfi_plot", "loc_id", "YEAR", "meas_num", "province")] %>% 
  pivot_wider(names_from = "meas_num", values_from = "YEAR", names_prefix = "X") %>% 
  mutate(REGION = ifelse(province %in% c("BC", "AB", "SK", "MB", "NT"), "WESTERN", "EASTERN"),
         nfi_plot_loc = paste(nfi_plot, loc_id, sep = "-"))

idx_order <- order(meas_dumb$province, desc(meas_dumb$X0), desc(meas_dumb$X1), desc(meas_dumb$X2))
fct_order <- meas_dumb$nfi_plot_loc[idx_order]
meas_dumb$nfi_plot_loc <- factor(meas_dumb$nfi_plot_loc,
                                levels = fct_order)

ggplot(filter(meas_dumb, REGION %in% c("EASTERN"))) +
  # First period
  geom_segment(aes(x = X0, xend = ifelse(is.na(X1), X0, X1),
                   y = nfi_plot_loc, yend = nfi_plot_loc),
               size = 1, color = "peachpuff3", alpha = 0.8) +
  # Firt measurment point
  geom_point(aes(X0, nfi_plot_loc), size = 3, color = "peachpuff2") +
  # Second period
  geom_segment(aes(x = X1, xend = ifelse(is.na(X2), X1, X2),
                   y = nfi_plot_loc, yend = nfi_plot_loc),
               size = 1, color = "peachpuff4", alpha = 0.8) +
  # Second measurement point
  geom_point(aes(ifelse(is.na(X1), X0, X1), nfi_plot_loc), size = 3, color = "peachpuff3") + 
  # Third point
  geom_point(aes(ifelse(is.na(X2), ifelse(is.na(X1), X0, X1), X2), nfi_plot_loc), size = 3, color = "peachpuff4") + 
  facet_grid(province~., scales = "free", space = "free") +
  labs(y = "NFI plot id", x = "Year", title = "Eastern Provinces") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey95", color = "transparent"))

# Measurement span by measurement type -------------------
# Include which had a change of location
# Include which were converted in Quebec
nfi_loc_change <- data2$nfi_plot[data2$loc_id == 1]
data2[data2$nfi_plot %in% nfi_loc_change,]
data3 <- data2 %>% mutate(REGION = ifelse(province %in% c("BC", "AB", "SK", "MB", "NT"), "WESTERN", "EASTERN"),
                          LOC_CHANGE = ifelse(nfi_plot %in% nfi_loc_change, "CHANGED", "UNCHANGED"),
                          CONVERTED = ifelse(is.na(user_info), "notconverted", 
                                             ifelse(str_detect(user_info, "Converted"), "converted", "notconverted")))
data3$PLOT_COMPL2 <- apply(data2, 1, function(x) sum(is.na(x["plotbio_lgtr_live"]),
                                                     is.na(x["plotbio_lgtr_dead"]),
                                                     is.na(x["plotbio_smtr_live"]),
                                                     is.na(x["plotbio_smtr_dead"]),
                                                     is.na(x["plotbio_lgshrub_live"]),
                                                     is.na(x["plotbio_lgshrub_dead"]),
                                                     is.na(x["plotbio_bryo"])
                                                     ))

meas_dumb2 <- data3[,c("nfi_plot", "loc_id", "YEAR", "meas_num", "province", "PLOT_COMPL2", "LOC_CHANGE", "CONVERTED", "REGION")] %>% 
  # pivot_wider(names_from = "meas_num", values_from = "YEAR", names_prefix = "X") %>%
  mutate(nfi_plot_loc = paste(nfi_plot, loc_id, sep = "-")) %>% 
  arrange(nfi_plot, meas_num, YEAR) %>%
  group_by(nfi_plot) %>% 
  mutate(nfi_plot = as.factor(nfi_plot),
         YEAR = as.integer(YEAR),
         YEAR_LEAD = lead(YEAR),
         MEAS_NUM2 = as.factor(row_number(nfi_plot)),
         MEAS_NUM3 = as.factor(as.numeric(MEAS_NUM2)+1))

md2_order <- meas_dumb2 %>% select(nfi_plot, province, MEAS_NUM2, YEAR) %>%  pivot_wider(names_from = MEAS_NUM2, values_from = YEAR, names_prefix = "X")
idx_order <- order(md2_order$province, desc(md2_order$X1), desc(md2_order$X2), desc(md2_order$X3))
fct_order <- md2_order$nfi_plot[idx_order]
meas_dumb2$nfi_plot <- factor(meas_dumb2$nfi_plot,
                                 levels = fct_order)

ggplot(data = filter(meas_dumb2, REGION == "EASTERN")) + 
  geom_segment(aes(y = nfi_plot, yend = nfi_plot, x = YEAR, xend = YEAR_LEAD, linetype = LOC_CHANGE, color = MEAS_NUM3, size = LOC_CHANGE)) +
  geom_segment(aes(y = nfi_plot, yend = nfi_plot, x = -Inf, xend = YEAR, color = CONVERTED), linetype = "dotted", size = 0.3) +
  geom_point(aes(YEAR, nfi_plot, color = MEAS_NUM2), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid")) +
  scale_size_manual(values = c("CHANGED" = 0.4, "UNCHANGED" = 1)) +
  scale_color_manual(values = c("peachpuff2", "peachpuff3", "peachpuff4", "red4", "grey80", "transparent")) +
  scale_x_continuous(breaks = c(seq(min(meas_dumb2$YEAR), max(meas_dumb2$YEAR), by = 1))) +
  facet_grid(province~., scales = "free", space = "free") +
  labs(y = "NFI plot id", x = "Year", title = "Eastern Provinces") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey95", color = "transparent"))

ggplot(data = filter(meas_dumb2, REGION == "WESTERN")) + 
  geom_segment(aes(y = nfi_plot, yend = nfi_plot, x = YEAR, xend = YEAR_LEAD, linetype = LOC_CHANGE, color = MEAS_NUM3, size = LOC_CHANGE)) +
  geom_point(aes(YEAR, nfi_plot, color = MEAS_NUM2), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid")) +
  scale_size_manual(values = c("CHANGED" = 0.4, "UNCHANGED" = 1)) +
  scale_color_manual(values = c("peachpuff2", "peachpuff3", "peachpuff4", "red4")) +
  scale_x_continuous(breaks = c(seq(min(meas_dumb2$YEAR), max(meas_dumb2$YEAR), by = 1))) +
  facet_grid(province~., scales = "free", space = "free") +
  labs(y = "NFI plot id", x = "Year", title = "Western Provinces") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey95", color = "transparent"))

# Measurement span by completeness   ------------------------------------------
meas_dumb3 <- meas_dumb2 %>% mutate(PLOT_COMPL3 = ifelse(PLOT_COMPL2 == 0, "COMPLETE", "INCOMPLETE"))

ggplot(data = filter(meas_dumb3, REGION == "EASTERN")) + 
  geom_segment(aes(y = nfi_plot, yend = nfi_plot, x = YEAR, xend = YEAR_LEAD, linetype = LOC_CHANGE, color = PLOT_COMPL3, size = LOC_CHANGE)) +
  geom_segment(aes(y = nfi_plot, yend = nfi_plot, x = -Inf, xend = YEAR, color = CONVERTED), linetype = "dashed", size = 0.6) +
  geom_point(aes(YEAR, nfi_plot, color = PLOT_COMPL3), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid")) +
  scale_size_manual(values = c("CHANGED" = 0.4, "UNCHANGED" = 1)) +
  scale_color_manual(values = c("forestgreen", "grey70", "red4", "transparent"), labels = c("Complete", "Converted plot", "Incomplete", "")) +
  scale_x_continuous(breaks = c(seq(min(meas_dumb2$YEAR), max(meas_dumb2$YEAR), by = 1))) +
  facet_grid(province~., scales = "free", space = "free") +
  labs(y = "NFI plot id", x = "Year", title = "Eastern Provinces", color = "Measurement\ncompleteness", linetype = "Location change", size = "Location change") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = "grey95", color = "transparent"),
        panel.background = element_blank())

ggplot(data = filter(meas_dumb3, REGION == "WESTERN")) + 
  geom_segment(aes(y = nfi_plot, yend = nfi_plot, x = YEAR, xend = YEAR_LEAD, linetype = LOC_CHANGE, color = PLOT_COMPL3, size = LOC_CHANGE)) +
  geom_point(aes(YEAR, nfi_plot, color = PLOT_COMPL3), size = 3) + 
  scale_linetype_manual(values = c("twodash", "solid"), labels = c("Changed", "Unchanged")) +
  scale_size_manual(values = c("CHANGED" = 0.4, "UNCHANGED" = 1), labels = c("Changed", "Unchanged")) +
  scale_color_manual(values = c("forestgreen", "red4"), labels = c("Complete", "Incomplete")) +
  scale_x_continuous(breaks = c(seq(min(meas_dumb2$YEAR), max(meas_dumb2$YEAR), by = 1))) +
  facet_grid(province~., scales = "free", space = "free") +
  labs(y = "NFI plot id", x = "Year", title = "Western Provinces", color = "Measurement\ncompleteness", linetype = "Location change", size = "Location change") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey95", color = "transparent"))

# Checking completeness per ecosystem compartment ------------------------------

compl_col <- data3 %>% 
  select(province,
         meas_num,
         LOC_CHANGE,
         REGION,
         plotbio_lgtr_live,
         plotbio_lgtr_dead,
         plotbio_smtr_live,
         plotbio_smtr_dead,
         plotbio_lgshrub_live,
         plotbio_lgshrub_dead,
         plotbio_herb,
         plotbio_bryo) %>%
  mutate(meas_num = as.factor(meas_num)) %>% 
  group_by(province, meas_num, LOC_CHANGE, REGION) %>% 
  summarise_all(function(x) sum(!is.na(x))) %>% 
  pivot_longer(names_to = "VAR", values_to = "N", plotbio_lgtr_live:ncol(.)) %>% 
  left_join(as.data.frame(table(data3[,c("province", "meas_num", "LOC_CHANGE", "REGION")]))) %>% 
  mutate(N_NA = Freq-N) %>% 
  select(-Freq) %>% 
  pivot_longer(names_to = "COMP_COUNT", values_to = "N", cols =  c("N", "N_NA"))

var_nms <- tibble(VAR = unique(compl_col$VAR),
                  VAR_NM = c("Lg Live Tree", "Lg Dead Tree", "Sm Live Tree", "Sm Dead Tree", "Lg Live Shrub", "Lg Dead Shrub", "Herbs", "Bryophytes"))

compl_col <- compl_col %>% 
  left_join(var_nms) %>% 
  mutate(VAR = factor(VAR, levels = unique(compl_col$VAR)),
         VAR_NM = factor(VAR_NM, levels = unique(var_nms$VAR_NM)))

ggplot(filter(compl_col, REGION == "WESTERN")) + 
  geom_col(aes(meas_num, N, fill = interaction(COMP_COUNT, LOC_CHANGE))) +
  scale_y_continuous(breaks = seq(0, 275, by = 25)) +
  scale_fill_manual(values = c("palegreen4", "tomato4", "forestgreen", "red4"), labels = c("Complete", "Missing", "Complete", "Missing")) +
  labs(x = "Variable", y = "Number of plots", fill = "Completeness") +
  facet_grid(province~VAR_NM, scales = "free", space = "free") +
  theme_minimal() +
  theme(panel.spacing.x = unit(0.25, "lines"))

ggplot(filter(compl_col, REGION == "EASTERN")) + 
  geom_col(aes(meas_num, N, fill = interaction(COMP_COUNT, LOC_CHANGE))) +
  scale_y_continuous(breaks = seq(0, 275, by = 25)) +
  scale_fill_manual(values = c("palegreen4", "tomato4", "forestgreen", "red4"), labels = c("Complete Changed", "Missing Changed", "Complete Unchanged", "Missing Unchanged")) +
  labs(x = "Variable", y = "Number of plots", fill = "Completeness") +
  facet_grid(province~VAR_NM, scales = "free", space = "free") +
  theme_minimal() +
  theme(panel.spacing.x = unit(0.25, "lines"))

# Final file contains extra metadata (region, completness, location change, converted plot)
write_csv(select(data3, -PLOT_COMPL2), "NFI_S03_metadata exploration/all_gp_site_info_MODIFIED.csv")
