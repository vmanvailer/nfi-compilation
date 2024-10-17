library(data.table)
library(tidyverse)

ltp_tree <- fread("00_data/all_gp_trees/all_gp_trees/all_gp_ltp_tree.csv")
sp_list <- fread("00_data/treesplistV45-added names.csv") %>% 
  mutate(CODE_GN_SP = paste0(CODE_GENU, ".", CODE_SPEC)) %>% select(CODE_GN_SP, GROUP) %>% unique()
ecoz_prov <- fread("02_utm-lambert conversion & retrieved ecoregions/all_gp_site_info_ecoclasses.csv",
                   select = c("nfi_plot", "meas_num", "loc_id", "ecozone", "province")) %>% unique



# To build threshold table we want to estimate biomass for a given DBH and HEIGHT.
# For Large Trees DBH threshold is 9cm but we need to estimate height.
# We can model height using Chapman-Richard's equation but we need to prepare the data first.

# CBM3 Biomass models are based on volume which, in turn, are based on tree shape (assumption).
# Trees commonly have broken tops and thus measured height of broken tree cannot be used 
# to model shape. So we must use projected height instead and later correct biomass
# for the broken top (The NFI procedure) 

# 1. Wrangling data ----------------------------------------------------------
## Creating species code for id generation and to match Lambert functions
ht_est <- ltp_tree[, SPECIES := paste0(lgtree_genus, ".", lgtree_species)]

# Creating a single column for height that uses projected height when available.
ht_est <- ecoz_prov[ht_est, on = names(ecoz_prov)[1:3]
                    ][, HEIGHT_TOTA := fifelse(!is.na(height_prj), height_prj, fifelse(
                      height < 0, NA_real_, height))
                      ][,DBH := fifelse(dbh < 0, NA_real_, dbh)
                        ][!is.na(DBH),
                          ][!is.na(HEIGHT_TOTA),]

ht_est <- ht_est[,.(nfi_plot, loc_id, meas_num, province, ecozone, SPECIES, meas_est_height, HEIGHT_TOTA, DBH)]
ht_est <- setnames(sp_list[ht_est, on = .(CODE_GN_SP = SPECIES)], "CODE_GN_SP", "SPECIES")
ht_est2 <- ht_est[meas_est_height %in% c("M", "S"),]
# 1.a. Creating groups for modelling data ------------------------------------

# Great we can now model height. Two things to consider: (1) tree N and (2) dbh range.
# (1) Depending on the grouping I decide to make for modelling, I can get more 
#     accurate but have smaller N (e.g. nfi id and species), or have higher N 
#     but more error as well e.g. (e.g. species and ecozones). 
# (2) Sites that don't have trees <15cm may be significantly off when extrapolating values. 
#     It will affect the "c" parameter of Chapman-Richard's and likely overestimate height (and thus biomass)

# NFI sets the minimum N threshold to three and does not consider model accuracy.
# I am unsure if it is better to use a 3 tree specific model or a 10+ tree general model.

# Let's see if we can aim for at least 10 values.
# Site / species
# Site / Veg type
# Province / Ecozone / Species
# Province / Ecozone / Veg type

# Minimum number of trees to include in the model 
mintrees = 10 

# Create all groups
glist <- list(
  G1 = c("nfi_plot", "SPECIES", "meas_num"),            # Model data 1 - Site / species
  G2 = c("nfi_plot", "GROUP", "meas_num"),              # Model data 2 - Site / Veg type
  G3 = c("province", "ecozone", "SPECIES", "meas_num"), # Model data 3 - Province / Ecozone / Species
  G4 = c("province", "ecozone", "GROUP", "meas_num"),   # Model data 4 - Province / Ecozone / Veg type
  G5 = c("ecozone", "SPECIES", "meas_num"),             # Model data 5 - Ecozone / species
  G6 = c("ecozone", "GROUP", "meas_num"),               # Model data 6 - Ecozone / Group
  G7 = c("ecozone", "province", "meas_num"),            # Model data 7 - Ecozone / Province
  G8 = c("province", "GROUP")
)


# # Flagging only sites that require modelling. No need to model sites without ingrowth.
# ingrowth_sites <- ht_est[GROWTH_COMP %in% c("INGROWTH", "INGROWTH_DEATH"),.(nfi_plot)] %>% unique() %>% as.matrix

# For a stronger model let's keep only height values that were measured "M" or that are unspecified "S"
# Create grouping columns used for filtering

# Use lapply to apply paste() to each set of columns (ChatGPT based-code with tweaks to work)
# This code essentially builds ids to group the data
model_data <- ht_est[meas_est_height %in% c("M", "S"),
                     ][,paste0("G", 1:length(glist)) := lapply(glist, function(x) do.call(paste, c(.SD[, x, with = FALSE], sep = ".")))]

# 1.b. Splitting data set in it's different groups ---------------------------
# create one column for each grouping pasting the names of each grouping
# modify mod_data function to spit out which groups of a given grouping should be used for creating a dataset
# Create a function that filters and split the data according to groups, storing results into a list.
# that list will be used to create the model.

# Function for filtering data used on the model based on minimum number of trees, and groupings.
data_aggr_list <- function(x){
  filt <- unique(model_data[,.N, by = c(x)][N>mintrees, 1])                  # For a given model level variable (e.g. G1), which groups (e.g. 634156.PICE.MAR.0) have more than N trees?
  dat_filt <- model_data[get(x) %in% as.matrix(filt),]                       # Subset only those groups
  cols_to_drop <- setdiff(grep("G\\d", names(dat_filt), value = TRUE), x)    # Get grouping not being used for filtering
  dat_filt <- dat_filt[,!cols_to_drop, with = FALSE]                         # Drop those columns
  fact <- dat_filt[,get(x)]                                                  # Create a separate variable to use for splitting the dataset by group. 
  mod_list <- split(dat_filt, fact)                                          # Apply the split
  return(mod_list)
}

list_mod_dfs <- lapply(names(glist), data_aggr_list)   # For every model grouping variable, subset and apply the split. 
names(list_mod_dfs) <- names(glist)                    # Rename grouping variables for easy tracking 

# 1.c. Data checking ---------------------------------------------------------
# A little verification code to make sure that the all model created will cover all nfi plots.
miss_mod1 <- model_data[,.N, by = G1, ][N<mintrees,] # Which one don't have enough tree
miss_mod2 <- model_data[G1 %in% miss_mod1$G1, .N, by = G2, ][N<mintrees,]
miss_mod3 <- model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2, .N, by = G3, ][N<mintrees,]
miss_mod4 <- model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2 & G3 %in% miss_mod3$G3, .N, by = G4, ][N<mintrees,]
miss_mod5 <- model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2 & G3 %in% miss_mod3$G3 & G4 %in% miss_mod4$G4, .N, by = G5, ][N<mintrees,]
miss_mod6 <- model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2 & G3 %in% miss_mod3$G3 & G4 %in% miss_mod4$G4 & G5 %in% miss_mod5$G5, .N, by = G6, ][N<mintrees, ]
miss_mod7 <- model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2 & G3 %in% miss_mod3$G3 & G4 %in% miss_mod4$G4 & G5 %in% miss_mod5$G5 & G6 %in% miss_mod6$G6, .N, by = G7, ][N<mintrees,]
miss_mod8 <- model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2 & G3 %in% miss_mod3$G3 & G4 %in% miss_mod4$G4 & G5 %in% miss_mod5$G5 & G6 %in% miss_mod6$G6 & G7 %in% miss_mod7$G7, .N, by = G8][N<mintrees,]

# Just a little look at how many sites were left out in which modelling group.
model_data[,.N, by = G1, ][N<mintrees, ]
model_data[G1 %in% miss_mod1$G1]
model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2,.N]
model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2 & G3 %in% miss_mod3$G3,.N]
model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2 & G3 %in% miss_mod3$G3 & G4 %in% miss_mod4$G4,.N]
model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2 & G3 %in% miss_mod3$G3 & G4 %in% miss_mod4$G4 & G5 %in% miss_mod5$G5,.N]
model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2 & G3 %in% miss_mod3$G3 & G4 %in% miss_mod4$G4 & G5 %in% miss_mod5$G5 & G6 %in% miss_mod6$G6,]
model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2 & G3 %in% miss_mod3$G3 & G4 %in% miss_mod4$G4 & G5 %in% miss_mod5$G5 & G6 %in% miss_mod6$G6 & G7 %in% miss_mod7$G7,]
model_data[G1 %in% miss_mod1$G1 & G2 %in% miss_mod2$G2 & G3 %in% miss_mod3$G3 & G4 %in% miss_mod4$G4 & G5 %in% miss_mod5$G5 & G6 %in% miss_mod6$G6 & G7 %in% miss_mod7$G7 & G8 %in% miss_mod8,]

# 2.a. Defining initial parameters for modelling height ---------------------- 
# Non-linear equations require good starting values for the parameters.
start.fun <- function(dat, height, dbh){
  q25 <- quantile(dat[[dbh]], 0.25, na.rm = TRUE)
  q50  <- quantile(dat[[dbh]], 0.5, na.rm = TRUE)
  q75 <- quantile(dat[[dbh]], 0.75, na.rm = TRUE)
  q100 <- quantile(dat[[dbh]], 1, na.rm = TRUE)
  dat_middle <- dat[dat[[dbh]] >= q25 & dat[[dbh]] <= q75, ]
  dat50 <- dat[dat[[dbh]] >= q50 & dat[[dbh]] <= q100, ]
  dat75 <- dat[dat[[dbh]] >= q75 & dat[[dbh]] <= q100, ]
  
  # Height at which dbh growth flattens
  a <- c(min(dat75[[height]], na.rm = TRUE),
         max(dat[[height]], na.rm = TRUE))
  
  # Growth rate i.e. steepness of the curve
  # b <- c(1/q75, q50)
  b <- c(lm(paste(height, " ~ ", dbh), data = dat75)$coefficients[2]*.3,
         lm(paste(height, " ~ ", dbh), data = dat_middle)$coefficients[2])
  
  # Growth rate at early stages. High value means slow growth when young
  c <- c(0.001, 5)
  start <- list(a = a, b = b, c = c)
  return(start)
}

# b = log((height[5] - height[1]) / (height[1] * dbh[5]^m)) / dbh[5]

start <- lapply(list_mod_dfs, function(x) lapply(x, start.fun, "HEIGHT_TOTA", "DBH"))

# X.a. Testing approaches to generate starting values based on data -----------
test_case <- list_mod_dfs[[1]][[4]][!is.na(HEIGHT_TOTA),]

test_start <- start.fun(dat = test_case,
                        height = "HEIGHT_TOTA",
                        dbh = "DBH")
model_form <- formula(HEIGHT_TOTA ~ 1.3+a*((1-exp(-b*DBH))^c))
# fix(test_start)
library(nls2)
ht_proj <- nls2(model_form,
                data = test_case,
                control = nls.control(maxiter = 2000),      # max iterations
                algorithm = "brute-force",                  # method for looking for starting values
                start = do.call(data.frame, test_start))    # starting values from brute force
ht_proj
# PLOT
ggplot(test_case, aes(DBH, HEIGHT_TOTA)) + 
  geom_point() + 
  geom_smooth(method = "nls",                               # Smooth line based on best fit.
              formula = y ~ 1.3 + a * ((1 - exp(-b * x))^c),
              method.args = list(ht_proj$m$getPars()),
              se = FALSE) +
  # geom_smooth(method = "lm",
  #             se = FALSE,
  #             color = "steelblue",
  #             linetype = "dashed") +
  lims(x = c(0, max(test_case$DBH)),
       y = c(0, max(test_case$HEIGHT_TOTA))) +
  theme_bw()


# Using nls2 coefs as starting values for nls.
brute_start <- ht_proj$m$getPars()
# fix(brute_start)
# a <- lm(HEIGHT_TOTA~DBH, data = test_case)

ht_proj2 <- nls(HEIGHT_TOTA ~ 1.3 + a*((1-exp(-b*DBH))^c), 
                data = test_case,
                start = brute_start) 
ht_proj2

ggplot(test_case, aes(DBH, HEIGHT_TOTA)) + 
  geom_point() + 
  geom_smooth(method = "nls",
              formula = y ~ 1.3 + a * ((1 - exp(-b * x))^c),
              method.args = list(ht_proj2$m$getPars()),
              se = FALSE,
              color = "red3") +
  lims(x = c(0, max(test_case$DBH)),
       y = c(0, max(test_case$HEIGHT_TOTA))) +
  theme_bw()

# Using minpack.lm to fit nls.
library(minpack.lm)
ht_proj3 <- nlsLM(formula = model_form,
                  data = test_case,
                  start = brute_start,
                  control = nls.lm.control(maxiter = 5000, maxfev = 5000))
ht_proj3

ggplot(test_case, aes(DBH, HEIGHT_TOTA)) + 
  geom_point() + 
  geom_smooth(method = "nls",
              formula = y ~ 1.3 + a * ((1 - exp(-b * x))^c),
              method.args = list(ht_proj3$m$getPars()),
              se = FALSE,
              color = "green3") +
  lims(x = c(0, max(test_case$DBH)),
       y = c(0, max(test_case$HEIGHT_TOTA))) +
  theme_bw()

rm(ht_proj, ht_proj2, ht_proj3)

# STEPS OUTLINED
# Create a function that takes the model_data, left_joins the groupings and keeps only the data to be model. i.e. line 104 plus filter
# Use that function to create a list of data used in the modelling according to each group.
# Create another function that uses nls2 to get good starting values for each model.
# Create another function that used those starting values in a nls equation.

# Run all models and store then in a list.

# Visualize all curves created to fish for odd things. Talk to Mihai about them.
# Once all is approved and correct, pass on the model on ingrowth trees to get height at DBH threshold.
# With the threshold table ready, pass it on to the LambertUng from CTAE package.
# Recalculate biomass for all available measurements and compare NFI to Lambert Ung.
# If nothing crazy shows up we are done with biomass calculations.
# Now summarize the difference between NFI and recalculations.
# Summarize large tree data as before in bins. 
# Summarize total AGB change for each plot and component. Stacked horizontal bars.
# Visualize it spatially.

# X.b. INTUITION ON NLS - Nonlinear least squares -------------------------------------------
test_case <- list_mod_dfs[[1]][[4]][!is.na(HEIGHT_TOTA),]
ht_proj <- nls(HEIGHT_TOTA ~ 1.3 + a*((1-exp(-b*DBH))^c), 
               data = test_case,
               start = list(a = 70, b = 0.01, c = 1.5))

# Can't find good starting values, trying brute-force on nls2
# It still requires a range of starting values. 
# Based on Ministry of Forestry pamphlets (Catolog of Curves for Curve Fitting), I could found what each parameter does.

# Based on chatGPT let's try to make functions that find good starting values.
# A is easy and is close to max height. Always a bit smaller so something between 0.9 of max height and max Height works
# B and C are harder and is what we are trying to figure it out.

library(nls2)
ht_proj <- nls2(HEIGHT_TOTA ~ 1.3 + a*((1-exp(-b*DBH))^c),     # Chapman-Richards eq.
                data = test_case,
                control = nls.control(maxiter = 1500),   # max iterations
                algorithm = "brute-force",               # method for looking for starting values
                start = data.frame(a = c(0.75*max(test_case$HEIGHT_TOTA),max(test_case$HEIGHT_TOTA)),         # Height at which dbh growth flattens
                                   b = c(0.01, 0.5),      # Growth rate i.e. steepness of the curve
                                   c = c(.01,7)))          # Growth rate at early stages. High value means slow growth when young
ht_proj # Awesome it worked.
# Using nls2 coefs as starting values for nls.
brute_start <- ht_proj$m$getPars()

ht_proj2 <- nls(HEIGHT_TOTA ~ 1.3 + a*((1-exp(-b*DBH))^c), 
                data = test_case,
                start = brute_start) # starting values from brute force
ht_proj2 # Slight improvement (898.9 vs 894.4. ~4 units of MRSE)

# Practical understanding of what which parameter does.
a <- 69.94811
b2 <- 0.01835  # Growth rate i.e. steepness of the curve
c2 <- 1.41508  # Growth rate at early stages. High value means slow growth when young

# Adding predicted values for plotting
test_case <- test_case %>% 
  mutate(FIT_MOD = 1.3 + a*((1-exp(-b2*DBH))^c2))

ggplot(test_case, aes(DBH, HEIGHT_TOTA)) + 
  geom_point() + 
  geom_smooth(method = "nls",                             # Smooth line based on best fit
              formula = y ~ 1.3 + a * ((1 - exp(-b * x))^c),
              method.args = list(ht_proj2$m$getPars()),
              se = FALSE) +
  geom_line(aes(y = FIT_MOD, x = DBH)) +                  # Modified line to check behavior
  lims(x = c(0, max(test_case$DBH)),
       y = c(0, max(test_case$HEIGHT_TOTA))) +
  theme_bw()

# X.b.animate -------------------------------------------------------------
test_case <- list_mod_dfs[[8]][[3]][!is.na(HEIGHT_TOTA),]
mod_group <- as.character(unique(test_case[,12]))

# Model it
test_start <- start.fun(dat = test_case,
                        height = "HEIGHT_TOTA",
                        dbh = "DBH")
model_form <- formula(HEIGHT_TOTA ~ 1.3+a*((1-exp(-b*DBH))^c))
# fix(test_start)
library(nls2)
ht_proj <- nls2(model_form,
                data = test_case,
                control = nls.control(maxiter = 2000),      # max iterations
                algorithm = "brute-force",                  # method for looking for starting values
                start = do.call(data.frame, test_start))    # starting values from brute force
ht_proj

# Fine tune
brute_start <- ht_proj$m$getPars()
# fix(brute_start)
# a <- lm(HEIGHT_TOTA~DBH, data = test_case)
ht_proj2 <- nls(HEIGHT_TOTA ~ 1.3 + a*((1-exp(-b*DBH))^c), 
                data = test_case,
                start = brute_start) 
ht_proj2

# Practical understanding of what which parameter does.
b <- ht_proj2$m$getAllPars()["b"]
c <- ht_proj2$m$getAllPars()["c"]
lenout <- 10
static <- (lenout*2)+((lenout-1)*2)
# 24 steps to cover entire range from origin and back to origin.
b2 <- c(seq(b, b*2, length.out = lenout),                # Up
        rev(seq(b, b*2, length.out = lenout)[-lenout]),  # Reverse up to starting point  
        seq(b, b/2, length.out = lenout)[-1],            # Down
        rev(seq(b, b/2, length.out = lenout)),           # Reverse down to starting point
        rep(b, static))                                        # Keep this static so it can vary c parameter

c2 <- c(rep(c, static),                                        # Keep this static so it can vary c parameter
        seq(c, c*2, length.out = lenout),                # Up
        rev(seq(c, c*2, length.out = lenout)[-lenout]),  # Reverse up to starting point  
        seq(c, c/2, length.out = lenout)[-1],            # Down
        rev(seq(c, c/2, length.out = lenout)))           # Reverse down to starting point

a <- rep(ht_proj2$m$getAllPars()["a"], length(c2))

# Vary both at the same time
# b2 <- c(seq(0.01835, 0.01835*3, length.out = 5))
# b2 <- c(c2, rev(c2[-5]), seq(0.01835, 0.01835/3, length.out = 5)[-1]) # Growth rate at early stages. High value means slow growth when young
# 
# c2 <- c(seq(1.41508, 1.41508*3, length.out = 5))
# c2 <- c(c2, rev(c2[-5]), seq(1.41508, 1.41508/3, length.out = 5)[-1]) # Growth rate at early stages. High value means slow growth when young

start_list <- data.frame(a,b2,c2) %>% split(f = 1:length(c2))

# Adding predicted values for plotting
test_case <- test_case %>% 
  mutate(FIT_MOD = 1.3 + start_list$`1`$a*((1-exp(-start_list$`1`$b2*DBH))^start_list$`1`$c2),
         TIME_STEP = 1,
         start_list$`1`)
test_case3 <- test_case
for (i in 2:length(c2)){
  start <- start_list[[i]]
  test_case2 <- test_case %>% 
    mutate(FIT_MOD = 1.3 + start$a*((1-exp(-start$b2*DBH))^start$c2),
           TIME_STEP = i,
           start)
  test_case3 <- rbind(test_case3, test_case2)
}

test_case4 <- test_case3 %>% mutate(LABEL = paste0("a = ", round(a, 3), "\nb = ", round(b2, 4), "\nc = ", round(c2, 4)))
lab_height <- max(test_case4$HEIGHT_TOTA)*.9
library(gganimate)
ggplot(test_case4, aes(DBH, HEIGHT_TOTA)) + 
  geom_point() + 
  geom_smooth(method = "nls",                             # Smooth line based on best fit
              formula = y ~ 1.3 + a * ((1 - exp(-b * x))^c),
              method.args = list(ht_proj2$m$getPars()),
              se = FALSE) +
  geom_line(aes(y = FIT_MOD, x = DBH)) +                  # Modified line to check behavior
  geom_text(data = unique(select(test_case4, TIME_STEP, LABEL)),
            aes(x = 0, y = lab_height, label = LABEL), hjust = 0) +
  lims(x = c(0, max(test_case4$DBH)),
       y = c(0, max(test_case4$HEIGHT_TOTA))) +
  transition_time(TIME_STEP) +
  labs(title = paste0("Chapman-Richards Curve for NFI model group: ", mod_group), subtitle = "Frame {frame} of {nframes}",
       x = "Diameter at breast height (cm)", y = "Height (m)") +
  ease_aes("linear") +
  theme_bw()

anim_save(filename = paste0("8_animated_params_", mod_group, ".gif"), animation = last_animation())

# X.c. Model lm() HEIGHT from DBH for the entire database --------------------
# Model
mod_lm <- function(df) {
  lm(HEIGHT_TOTA~DBH, data = df)
}
# separate each site, location and growth component.
ht_est_m <- ht_est %>% 
  group_by(nfi_plot, ecozone, SPECIES) %>% 
  select(DBH, HEIGHT_TOTA) %>% 
  nest() %>% 
  filter(map(data, nrow)>=10) %>% 
  mutate(DBH_MOD = map(data, mod_lm))

map_dbl(ht_est_m$data, min) %>% min()
# Extract coefficients 
b_fun <- function(mod){
  beta <- coefficients(mod)[[2]]
}

msre_fun <- function(mod){
  msre <- mean(residuals(mod)^2)
}  
# separate into a new table
ht_dbh_coef <- transmute(ht_est_m, 
                         beta = map_dbl(DBH_MOD, b_fun),
                         msre = map_dbl(DBH_MOD, msre_fun)) %>% ungroup

# Predicting height at dbh N cm
N = 15
predict(ht_est_m$DBH_MOD[[1]], newdata = data.frame(DBH = N))
ht_est_m <- ht_est_m %>% mutate(HEIGHT_AT_NCM_DBH = map_dbl(DBH_MOD, predict, newdata = data.frame(DBH = N)))

# Visualizing model and investigating rare values
ht_est$lgtree_genus %>% table %>% sort
nfiid <- 1452016
test <- filter(ht_est, nfi_plot == nfiid)

test %>% 
  ggplot(aes(DBH, HEIGHT_TOTA)) + 
  geom_point(alpha = 0.1) + 
  
  # Predicted value
  geom_point(data = ht_est_m[ht_est_m$nfi_plot == nfiid,],
             aes(x = N, y = HEIGHT_AT_NCM_DBH), size = 3) +
  
  geom_smooth(se = F, method = "lm") +
  facet_grid(ecozone~SPECIES) +
  theme_minimal()

# Linear model overestimate height at smaller dbh.
# High RMSE are generally associated with very logarithmic like curves
# Small RMSE are generally associated with calculated values (as opposed to measured)
# Overall should be enough to continue.

# 2.b. Creating height models based on groups and starting parameters ---------
model_form <- formula(HEIGHT_TOTA ~ 1.3 + a*((1-exp(-b*DBH))^c))
start <- lapply(list_mod_dfs, function(x) lapply(x, start.fun, "HEIGHT_TOTA", "DBH"))

library(nls2)

model_fit_list <- list()
for (i in seq_along(list_mod_dfs)){
  model_fit_list[[i]] <- list()
  for (j in seq_along(list_mod_dfs[[i]])){
    
    dat <- list_mod_dfs[[i]][[j]]
    model_nls2 <- nls2(model_form,
                       data = dat,
                       control = nls.control(maxiter = 2000),
                       algorithm = "brute-force",
                       start = do.call(data.frame, start[[i]][[j]]),
                       )
    
    brute_start <- model_nls2$m$getPars()
    
    # THIS FUNCTIONS REQUIRES A TRY CATCH BLOCK.
    # TO BETTER PREDICT STARTING POINTS A GROWTH AND YIELD MODELLER MAY BE REQUIRED.
    
    ht_proj <- nls(formula = model_form, 
                   data = dat,
                   start = brute_start) 
    model_fit_list[[i]][[j]] <- ht_proj
    print(names(list_mod_dfs[i]))
    print(names(list_mod_dfs[[i]][j]))
  }
}

#
#
#
#
# Stopped here
#
#
#
#
#

# 2.c. Now build threshold table -------------------------------------------------
trshd <- ht_est_m %>% ungroup %>% select(nfi_plot:SPECIES, HEIGHT_AT_NCM_DBH)
# 1066 groups modeled for height

# Let's see if there are ingrowth cases that could not have height modeled.
trshd <- ht_est %>% 
  select(SPECIES, ecozone, GROWTH_COMP, nfi_plot) %>% 
  unique() %>% 
  filter(GROWTH_COMP %in% c("INGROWTH", "INGROWTH_DEATH")) %>% 
  left_join(trshd)

# 3,422 combinations of nfi_plot, ecozone and SPECIES
# from those 1178 have an INGROWTH component
is.na(trshd$HEIGHT_AT_NCM_DBH) %>% sum()
1178-640

# from those I could model height for 538 combinations
#
#
#
#
#
#
#
#

# Think about how to model those 640 sites later
#
#
#
#
#
#
#
#

# 3.   Now onto the AGB --------------------------------------------

# There are two ways to do that:
# (1) Recalculate volume using volume to biomass formulas (complicated)
# (2) Calculate biomass directly from DBH and HEIGHT using standardized national equations from Lambert and Ung (simpler and similar to NFI)

# Using package from Piotr and Juha for biomass estimations

# Check which species are not in Lambert-Ung table
common <- intersect(unique(trshd$SPECIES), unique(parameters_LambertUng$species)) %>% sort()

ltp_tree <- ltp_tree %>% 
  mutate(SPECIES = paste0(lgtree_genus, ".",lgtree_species),
         CTAE_HT = ifelse(stem_cond == "B", height_prj, height),
         CTAE_AD = ifelse(is.na(vol_total), 1, vol_total)/ifelse(is.na(vol_prj), vol_total, vol_prj))

x <- ltp_tree[ltp_tree$SPECIES %in% common,] %>% 
  left_join(trshd) %>% 
  filter(!is.na(HEIGHT_AT_NCM_DBH))

x_CTAE <- tibble()
start <- Sys.time()
for (i in 1:nrow(x)){
  t1 <- Sys.time()
  
  # Recalculating biomass based on Lambert 2005 and Ung 2006.
  agb <- CTAE::AGB_LambertUngDBHHT(DBH = x$dbh[i],
                                   height = x$CTAE_HT[i],
                                   species = x$SPECIES[i]) %>% 
    reshape2::melt(value.name = "LAMBERTUNG") %>% 
    mutate(LAMBERTUNG_ADJ = LAMBERTUNG*x$CTAE_AD[i]) # Adjusting for broken top
  
  # Calculating biomass at threshold based on Lambert 2005 and Ung 2006.
  agb_trshd <- CTAE::AGB_LambertUngDBHHT(DBH = N,
                                         height = x$HEIGHT_AT_NCM_DBH[i],
                                         species = x$SPECIES[i]) %>% 
    reshape2::melt(value.name = "LAMBERTUNG_TRSHD") %>% 
    mutate(LAMBERTUNG_TRSHD_ADJ = LAMBERTUNG_TRSHD*x$CTAE_AD[i]) # Adjusting for broken top
  
  # Binding both calculations together
  agb_all <- cbind(x[i,1:6],merge(agb, agb_trshd)) %>% 
    rename("BIO_COMP" = "L1") %>% 
    # Setting crown biomass to zero is dead.
    mutate(LAMBERTUNG_ADJ = ifelse(GROWTH_COMP == "INGROWTH_DEATH" & 
                                     BIO_COMP %in% c("Bcrown", "Bbranch", "Bfoliage"), 0, LAMBERTUNG_ADJ),
           LAMBERTUNG_TRSHD_ADJ = ifelse(GROWTH_COMP == "INGROWTH_DEATH" & 
                                           BIO_COMP %in% c("Bcrown", "Bbranch", "Bfoliage"), 0, LAMBERTUNG_TRSHD_ADJ),
           GR_SINCE_TRSHD = LAMBERTUNG_ADJ - LAMBERTUNG_TRSHD_ADJ,
           BIOMASS_NFI = x$biomass_total[i],
           dbh = x$dbh[i],
           height = x$CTAE_HT[i],
           dbh_increase = x$dbh[i] - N,
           height_increase = x$CTAE_HT[i] - x$HEIGHT_AT_NCM_DBH[i],
           meas_est_height = x$meas_est_height[i],
           stem_cond = x$stem_cond[i],
           SPECIES = x$SPECIES[i])
  # Adding to the master table
  x_CTAE <- rbind(x_CTAE, agb_all)
  t2 <- Sys.time()
  # Just a print to track the loop.
  print(paste("nfi_id =", x[i,"nfi_plot"], " | ", "tree_num = ", x[i,"tree_num"], "|", t2-t1))
}
end <- Sys.time()
end-start

ggplot(filter(x_CTAE, BIO_COMP == "Btotal")) + 
  geom_density(aes(x = GR_SINCE_TRSHD, fill = GROWTH_COMP), alpha = 0.4) +
  # scale_x_continuous(limits = c(0,50)) +
  # geom_ridgeline(aes(x = GR_SINCE_TRSHD, y = BIO_COMP, height = 0.5, fill = GROWTH_COMP), alpha = 0.4) +
  # stat_smooth(aes(biomass_total, Btotal, color = CTAE_SP), method = "lm") + 
  theme_minimal()

x_CTAE$meas_est_height <- factor(x_CTAE$meas_est_height, 
                                 levels = c("M", "C", "E", "S"))
ggplot(filter(x_CTAE, BIO_COMP == "Btotal"), 
       aes(x = dbh_increase, y = height_increase, color = SPECIES)) + 
  geom_point(size = 2, alpha = 0.4) + 
  stat_smooth(method = "loess") + 
  theme_minimal()


ggplot(filter(x_CTAE, BIO_COMP == "Btotal"), aes(x = BIOMASS_NFI, y = LAMBERTUNG_ADJ, color = meas_est_height)) + 
  geom_point() + 
  stat_smooth(method = "loess")

comp <- lm(LAMBERTUNG_ADJ ~ BIOMASS_NFI, data = filter(x_CTAE, BIO_COMP == "Btotal"))
plot(comp$model$LAMBERTUNG_ADJ ~ comp$model$BIOMASS_NFI)

plot(residuals(comp))
abline(h = 0)

msre <- mean(residuals(comp)^2)
