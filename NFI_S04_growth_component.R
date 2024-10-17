library(tidyverse)

ltp_tree <- read_csv("00_data/all_gp_trees/all_gp_trees/all_gp_ltp_tree.csv")
ltp_tree$lgtree_status %>% table

# Flagging function
fun_ingrowth <- function(status_bef, status_aft = NULL){
  
  case_when(
    # Condition 1 - Ingrowth (trees that crossed the threshold in between measurements and survived)
    is.na(status_bef) & status_aft %in% c("LS", "LF")                        ~ "INGROWTH",
    # Condition 2 - Crossed the threshold and died in between measurements
    is.na(status_bef) & status_aft == "DS"                                   ~ "INGROWTH_DEATH",
    # Condition 3 - Survived
    status_bef %in% c("LS", "LF", "M", "DS") & status_aft %in% c("LS", "LF") ~ "SURVIVED",
    # Condition 4 - Dead in both census (assumed dead is data is missing on second census)
    status_bef == "DS" & status_aft %in% c("DS", "M")                        ~ "DEAD",
    # Condition 5 - Live on first census but died before second
    status_bef %in% c("LS", "LF") & status_aft == "DS"                       ~ "DIED",
    .default = NA
  )
  
}

# Flagging growth component of trees

treen <- ltp_tree %>% 
  select(nfi_plot, loc_id, meas_num, tree_num, lgtree_status) %>% 
  pivot_wider(names_from = "meas_num",
              names_prefix = "status_",
              values_from = "lgtree_status") %>% 
  mutate(GROWTH_COMP_0 = ifelse(
    is.na(status_0), "NOT_PRESENT", ifelse(
      status_0 %in% c("LS", "LF", "M", "DS") & status_1 %in% c("LS", "LF"), "SURVIVED", ifelse(
        status_0 %in% c("LS", "LF"), "SURVIVED", ifelse(  
          status_0 == "DS", "DEAD", NA)
      ))),
    GROWTH_COMP_1 = fun_ingrowth(status_0, status_1),
    GROWTH_COMP_2 = fun_ingrowth(status_1, status_2))


treen
treen$GROWTH_COMP_0 %>% table
treen$GROWTH_COMP_0 %>% table %>% sum
treen %>% filter(is.na(GROWTH_COMP_0)) %>% print(n=41)

# 619 plots flagged from a total of 641 plots that were remeasured. Reconcile later...
treen %>% filter(!is.na(GROWTH_COMP_1)) %>% select(nfi_plot) %>% unique

treen$GROWTH_COMP_1 %>% table %>% sum # 23331 trees flagged

# From the total of *plots measured before and after*,  what is the percentage of each flag
treen$GROWTH_COMP_1 %>% table %>% as_tibble() %>% mutate(PERC = (n/sum(n))*100)

# Adding growth component table to original table
ltp_tree_grwcomp <- treen %>% 
  pivot_longer(cols = status_0:ncol(treen),
               names_to = c(".value", "meas_num"),
               names_pattern = "(.*)_(\\d)",
               names_transform = list(meas_num = as.numeric)) %>% 
  right_join(ltp_tree)

write_csv(ltp_tree_grwcomp, "NFI_S04_ltp_tree_grwcomp.csv")
