library(tidyverse)

# There are inconsistencies with species name e.g. Alnus as softwood when it should be hardwood
# Info at Genus level is enough to figure that out
sp_name_list <- read_csv("00_data/treesplistV45-added names.csv")

# Let's see how many inconsistencies
genus_unique <- unique(sp_name_list[,c("GROUP", "CODE_GENU")])
duplicated(genus_unique$CODE_GENU) %>% sum 
# 21 instances where genera are classified as both soft and hardwood.

# Let's fix this manually based on https://powo.science.kew.org/
# For those 21 instances let's create a csv file (03_genus_to_fix.csv) so we can
# edit based on what the Kew website says about which group a given genera belongs to.
dups <- genus_unique[duplicated(genus_unique$CODE_GENU),]$CODE_GENU
sp_name_list[sp_name_list$CODE_GENU %in% dups,] %>% arrange(CODE_GENU) %>% write_csv("00_data-fix/00_data-fix_genus_to_fix.csv")
# It seems that all exotic species have also been marked as their opposite group.

# Checked everything in the 03_genus_to_fix.csv against the Kew website 
# and created the file 03_genus_to_fixed.csv
genus_fixed <- read_csv("00_data-fix/00_data-fix_genus_fixed.csv", 
                        col_names = c("GROUP_FIXED", "CODE_GENU"))

# Now substitute old group calssification by new fixed one.
sp_name_list_fixed <- sp_name_list %>% 
  left_join(genus_fixed, by = join_by("CODE_GENU")) %>% 
  mutate(GROUP = ifelse(!is.na(GROUP_FIXED), GROUP_FIXED, GROUP)) %>% 
  select(-GROUP_FIXED)

write_csv(sp_name_list_fixed, "00_data-fix/00_data-fix_sp_name_list_fixed.csv")
