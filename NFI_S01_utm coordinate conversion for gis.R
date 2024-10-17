library(tidyverse)
library(terra)

# Final file contains extra metadata (region, completness, location change, converted plot)
data <- read_csv("00_data/all_gp_site_info.csv")

# Conversion must be done for each utm zone. 
zcrs <- as.list(paste0("+proj=utm +zone=", sort(unique(data$utm_zone)), " +datum=WGS84 +units=m"))
zdata <- split(data, data$utm_zone)

zdata2 <- map2(zdata, zcrs, function(x, y) vect(x, geom = c("utm_e", "utm_n"), crs = y))
# Projecting on lambert canonical (https://proj.org/en/9.4/operations/projections/lcc.html)
# Don't remember where got parameters from. Likely a default version from QGIS.
zdata3 <- map(zdata2, project, "+proj=lcc +lat_0=63.390675 +lon_0=-91.8666666666667 +lat_1=49 +lat_2=77 +x_0=6200000 +y_0=3000000 +datum=NAD83 +units=m +no_defs")
zdata4 <- reduce(zdata3, rbind)
zdata5 <- cbind(terra::geom(zdata4, df = TRUE)[,c("x", "y")],
                terra::as.data.frame(zdata4)) %>% 
  rename("x_lambert" = "x", "y_lambert" = "y")

write_csv(zdata5, "2_UTM-Lambert conversion & retrieved ecoregions/all_gp_site_info_coord_lamb.csv")
write_csv(select(zdata5, x_lambert, y_lambert, nfi_plot, LOC_CHANGE, CONVERTED), "02_utm-lambert conversion & retrieved ecoregions/NFI_S01nfi_coord_lamb.csv")

# Project to latlong
zdata3 <- map(zdata2, project, "+proj=longlat +datum=WGS84 +no_defs")
zdata4 <- reduce(zdata3, rbind)
zdata5 <- cbind(terra::geom(zdata4, df = TRUE)[,c("x", "y")], terra::as.data.frame(zdata4)) %>% rename("x_long" = "x", "y_lat" = "y")

write_csv(select(zdata5, x_long, y_lat, nfi_plot, LOC_CHANGE, CONVERTED), "02_utm-lambert conversion & retrieved ecoregions/NFI_S01nfi_coord_latlon.csv")

# geodesic distance of sites that had location change

zdata6 <- zdata5 %>% filter(meas_num == 0)  %>% select(nfi_plot, x_long, y_lat, loc_id) %>%  pivot_wider(names_from = "loc_id", values_from = c("x_long", "y_lat"), names_prefix = "Z") %>% filter(!is.na(x_long_Z1))
dist <- geodist::geodist(zdata6[,c("x_long_Z0", "y_lat_Z0")],
                         zdata6[,c("x_long_Z1", "y_lat_Z1")],
                         paired = TRUE, measure = "geodesic") %>% as_tibble()

dist2 <- cbind(nfi_plot = zdata6$nfi_plot, dist) %>% left_join(data[,c("nfi_plot", "province")])

ggplot(filter(dist2, province != "ON")) + 
  geom_boxplot(aes(province, value, fill = province)) + 
  scale_fill_manual(values = c("grey30", "grey60")) +
  theme_minimal()

write_csv(dist2, "02_utm-lambert conversion & retrieved ecoregions/NFI_S01loc_change_geodist.csv")


# geodesic distance of sites that do not flag location change but have different coordinates


zdata6 <- zdata5 %>% 
  filter(loc_id == 0)  %>% select(nfi_plot, x_long, y_lat, meas_num) %>%  
  pivot_wider(names_from = "meas_num", values_from = c("x_long", "y_lat"), names_prefix = "Z") %>% 
  filter(!is.na(x_long_Z1))
dist01 <- geodist::geodist(zdata6[,c("x_long_Z0", "y_lat_Z0")], zdata6[,c("x_long_Z1", "y_lat_Z1")], paired = TRUE, measure = "geodesic") %>% as_tibble() %>% rename("DIST_0_1" = "value")
dist12 <- geodist::geodist(zdata6[,c("x_long_Z1", "y_lat_Z1")], zdata6[,c("x_long_Z2", "y_lat_Z2")], paired = TRUE, measure = "geodesic") %>% as_tibble() %>% rename("DIST_1_2" = "value")
dist2 <- cbind(nfi_plot = zdata6$nfi_plot, dist01) %>% cbind(dist12) %>% left_join(unique(data[,c("nfi_plot", "province")]))

dist2[,2:3] <- lapply(dist2[,2:3], function(x) {
  a <- format(x, scientific = FALSE, digits = 5) 
  b <- as.numeric(a)
}
  )

write_csv(dist2, "02_utm-lambert conversion & retrieved ecoregions/NFI_S01_loc_unchanged_geodist.csv")


ggplot(filter(dist2, DIST_0_1 > 1.5 & DIST_0_1 < 300),
       aes(province, DIST_0_1)) + 
  geom_boxplot(outlier.colour = "grey10", outlier.size = 2, outlier.shape = 15, fill = "grey95", color = "grey75") + 
  geom_jitter(alpha = 0.3) +
  labs(x = "Province", y = "Distance (m) between M0 and M1") +
  scale_y_continuous(breaks = c(seq(0,800, 50))) +
  theme_minimal() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggplot(filter(dist2, DIST_0_1 > 300),
       aes(province, DIST_0_1/1000)) + 
  geom_jitter(alpha = 0.7, size = 3, width = 0.2) +
  geom_line(alpha = 0.3, size = 2) +
  scale_y_continuous(breaks = c(seq(0,400, 50))) +
  labs(x = "Province", y = "Distance (km) between M0 and M1") +
  theme_minimal() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank())

filter(dist2, DIST_0_1 > 300) %>% mutate(DIST_0_1_KM = DIST_0_1/1000) %>% write_csv("02_utm-lambert conversion & retrieved ecoregions/NFI_S01_above300m_loc_diff.csv")

ggplot(filter(dist2, DIST_1_2>0),
       aes(province, DIST_1_2)) + 
  geom_boxplot(outlier.colour = "grey10", outlier.size = 2, outlier.shape = 15, fill = "grey95", color = "grey75") + 
  geom_jitter(alpha = 0.3) +
  labs(x = "Province", y = "Distance (m) between M1 and M2") +
  theme_minimal() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
