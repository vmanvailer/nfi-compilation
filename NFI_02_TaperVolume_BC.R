#' Calculates merchantable tree-level volume using Kozak (1997) and Smalian's equations and tree DBH and height as input.
#'
#' @description Calculates tree-level volume [m3] using taper equation 2 described in Kozak (1988).
#' @param province Province where the tree is located.
#' @param ecozone Ecozone where the tree is located.
#' @param beczone (Optional) If province is BC, which BEC zone it is located.
#' @param taper code for the taper equation to use either "BC" or "QCI".
#' @param species Species name folowing NFI standandars e.g. "PINU.CON".
#' @param dbh Measured Tree DBH in [cm].
#' @param height Tree Height in [m].
#' @return A single numerical value representing total merchantable tree volume in [m3]. 
#' 
#' @references
#' Kozak, A. 1997. Effects of multicollinearity and autocorrelation on the variable-exponent taper functions.
#' Can. J. For. Res. 27:619-629. https://doi.org/10.1139/x97-011
#' 
#' @examples
#' 
#' taper_vol_kozak(province = "BC",
#'                 ecozone  = "13",
#'                 beczone  = "SBS",
#'                 taper    = "BC",
#'                 species  = "PINU.CON", 
#'                 dbh      = 30.1,
#'                 height   = 25.9)
#' 
#' @export


# Example is:
# Lodgepole Pine - PL
# BEC Sub-Boreal Spruce (SBS)
# nfi_plot 1328086
# tree_num 101 (Live Standing no broken top)
# NFI calculated volume is 0.8135 m3

# Example parameters
a0 <- 1.011506
a1 <- 1.05571
a2 <- 0.996782
b0 <- 93.9356
b1 <- -244.593
b2 <- 234.778
b3 <- -83.3503
b4 <- -50.4453
b5 <- -0.3566
b6 <- 0.000373
DBH <- 30.1 # cm
H <- 25.9 # m
p <- 0.01

taper_vol_kozak <- function(province, ecozone, beczone, taper, species, dbh, height){

# REQUIRES THE CREATION OF A LOOKUP OF COEFFICIENTS  
  
# Helper functions
calculate_xi <- function(hi, H) {
  (1 - (hi/H)^0.5)/(1 - p^0.5)
}

calculate_qi <- function(hi, H) {
  1 -(hi/H)^0.5
}

# Diameter inside the bark based on Eq 2 of Kozak 1997
calculate_di <- function(DBH, hi, H) {
  zi <- hi/H
  xi <- calculate_xi(hi, H)
  qi <- calculate_qi(hi, H)
  a0 * (DBH^a1) * (a2^DBH) * xi^(b0 + b1 * (zi^0.25) + b2 * (zi^(1/3)) + b3 * (zi^0.5) + b4 * asin(qi) + b5 * (1 / (DBH/H + zi)) + b6 * H)
}

# Smalian's formula
calculate_volume_slice <- function(d_top, d_bottom, length) {
  (pi * (d_bottom^2 + d_top^2)/20000) * length
}

# Step 1. Calculate the diameter inside bark (Di) at stump height (0.3 m)
stump_height <- 0.3 # m
stump_diameter <- calculate_di(DBH, stump_height, H)

# Step 2. Calculate the merchantable height (MH) which is the tree height at which the top diameter (Di) is 10.0 cm.
find_MH <- function(hi) {
  zi <- hi/H
  Xi <- (1 - (hi/H)^0.5) / (1 - 0.01^0.5)
  Qi <- 1 - (hi/H)^0.5
  Di <- a0 * (DBH^a1) * (a2^DBH) * Xi^(b0 + b1 * (zi^0.25) + b2 * (zi^(0.33333)) + b3 * (zi^0.5) + b4 * asin(Qi) + b5 * (1/(DBH/H + zi)) + b6 * H)
  return(Di - 10)
}

MH <- uniroot(find_MH, c(0.3, H))$root

# Step 3. Calculate the volume for stump which is assumed to be a cylinder with diameter equal to the stump diameter.
stump_volume <- (pi * (stump_diameter/2)^2/10000) * stump_height

# Step 4. Split the merchantable length, which is the length between the stump and 
# the merchantable height, into 10 cm slices and calculate the top diameters for each slice.
merchantable_length <- MH - stump_height
number_of_slices <- merchantable_length/0.1
top_diameters <- numeric(number_of_slices)

for (i in 1:number_of_slices) {
  hi <- stump_height + i * 0.1
  top_diameters[i] <- calculate_di(DBH, hi, H)
}

# Step 5: Calculate the volume for each slice using Smalian’s formula
bottom_diameters <- c(stump_diameter, top_diameters[1:(number_of_slices - 1)])
slices_volumes <- numeric(number_of_slices)
for (i in 1:number_of_slices) {
  slices_volumes[i] <- calculate_volume_slice(top_diameters[i], bottom_diameters[i], 0.1)
}


# Step 6. Split the Tree Top Length, which the length between the merchantable 
# height and Total Tree Height (H) into 10 cm slices and calculate the volume of
# each slice using the Smalian's formula.
height_last_slice <- (stump_height + floor(number_of_slices) * 0.1)
tree_top_length <- H - height_last_slice
number_of_slices_tree_top <- tree_top_length/0.1
top_diameters_tree_top <- numeric(number_of_slices_tree_top)

for (i in 1:number_of_slices_tree_top) {
  hi <- height_last_slice + (i) * 0.1
  top_diameters_tree_top[i] <- calculate_di(DBH, hi, H)
}

bottom_diameters_tree_top <- c(top_diameters[number_of_slices], top_diameters_tree_top[1:(number_of_slices_tree_top - 1)])
tree_top_volume <- 0
for (i in 1:number_of_slices_tree_top) {
  tree_top_volume <- tree_top_volume + calculate_volume_slice(top_diameters_tree_top[i], bottom_diameters_tree_top[i], 0.1)
}


# Step 7. Calculate the total volume of the merchantable length by summing up the volumes of each slice.
merchantable_volume <- sum(slices_volumes)


# Step 7. Calculate the total volume of the tree by adding the volume of the stump to the volume of the merchantable length.
total_volume <- stump_volume + merchantable_volume + tree_top_volume
}
