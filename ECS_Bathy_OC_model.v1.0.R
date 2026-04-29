# ==============================================================================
# Script Name: ECS_Bathy_model.R
# Description: Mechanistic modeling of microbial-carbon dynamics in deep-sea 
#              sediments. This model simulates the physiological trade-offs 
#              driving archaeal dominance and carbon turnover in the deep subsurface.
# manuscript: "Physiological trade-offs drive the archaeal 
#              dominance and carbon turnover in deep subsurface"
# Code Author: Jialin Hou (houjialin6@sjtu.edu.cn)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment Setup & Data Loading
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(deSolve) 
  library(tidyverse) 
  library(readr) 
  library(ggplot2) 
  library(patchwork)
  library(scales)
  library(ggpubr)
  library(ggsci)
})

# Sedimentation constants and unit conversions
# - sedimentation: 0.61 cm/yr (Yang et al., 2023)
# - bulk dry density: 1.1 g/cm^3 (Yang et al., 2023)
# - cell dry weight: 24 fg/cell (Bar-on et al., 2018)
# Cell carbon conversion: 24 fg C/cell -> ugC/cm3 = cells/g * 24 fg * 1.1 g/cm3 * 1e-9

# Load biological observations
biom <- read.table("./Estimated_cell_num.txt", sep = "\t", header = TRUE)
observed.biomass.df <- as_tibble(biom) %>%
  mutate(
    time = Depth / 0.61,   
    Arc_cell   = Arc * 24 * 1.1 * 1e-9,
    Bac_cell   = Bac * 24 * 1.1 * 1e-9,
    Bathy_cell = Bathy * 24 * 1.1 * 1e-9,
    Total_cell = Arc_cell + Bac_cell
  )

# Load Total Organic Carbon (TOC) observations
toc <- read.table("./TOC.txt", header = TRUE, sep = "\t")
observed.POC.df <- as_tibble(toc) %>%
  mutate(
    time = Depth / 0.61,  
    observed_POC = TOC * 1.1 * 1e4
  )

# ------------------------------------------------------------------------------
# 2. Model Parameters & Initial Conditions
# ------------------------------------------------------------------------------
# Ecological Strategies:
# B1: Copiotroph (surface) - High V_max, High mq, Low YG, High alpha
# B2: Generalist (mid-depth) - Slower growth, prefers POC2
# A1: Mid-depth Archaea - Conservative, moderate V_max and maintenance
# A2: Oligotroph (deep Bathyarchaeia) - Extremely slow growth, high efficiency, 
#     minimal maintenance, high survival rate.

parms <- list(
  # Maximum growth rates (kyr-1, Bradley et al., 2019)
  V_max_B1 = 60,   V_max_A1 = 50,
  V_max_B2 = 4,    V_max_A2 = 2,
  
  # Half-saturation constant for growth (ugC/cm3)
  K_V_B1_L = 6300, K_V_A1_L = 6800,
  K_V_B2_R = 8000, K_V_A2_R = 4000,
  
  # True growth yield (substrate conversion efficiency)
  YG_B1 = 0.006,   YG_A1 = 0.003,
  YG_B2 = 0.012,   YG_A2 = 0.002, 
  
  # Mortality rates (per kyr, Bradley et al., 2019)
  alpha_B1 = 7e-2, alpha_A1 = 1e-2,
  alpha_B2 = 2e-2, alpha_A2 = 8e-4,
  
  # Threshold POC for maintenance energy provenance (ug C/cm3)
  K_M_B1 = 600,    K_M_A1 = 400, 
  K_M_B2 = 3000,   K_M_A2 = 2000,  
  
  # Maintenance demands (per 1000 yr, Bradley et al., 2019)
  mq_B1 = 0.3,     mq_A1 = 0.3, 
  mq_B2 = 0.1,     mq_A2 = 0.1, 
  
  # Steepness of maintenance power dependency
  st_M_B1 = 0.1,   st_M_A1 = 0.1,  
  st_M_B2 = 0.1,   st_M_A2 = 0.1, 
  
  # Necromass recycling parameters
  death_to_POC = 0.4,
  death_to_POC1_frac = 0.9
)

# Initial states based on surface sediment observations
start <- c(
  B1 = 0.5, B2 = 0.1, A1 = 0.08, A2 = 0.08,           # Biomass
  POC1 = 3000, POC2 = 5800,                           # POC pools
  
  # Cumulative trackers
  dc_POC_growth_total = 0, dc_POC_Growth_B1 = 0, dc_POC_Growth_B2 = 0, dc_POC_Growth_A1 = 0, dc_POC_Growth_A2 = 0,
  dc_M_total_POC = 0, dc_M_POC_B1 = 0, dc_M_POC_B2 = 0, dc_M_POC_A1 = 0, dc_M_POC_A2 = 0,
  dc_M_total_bio = 0, dc_M_bio_B1 = 0, dc_M_bio_B2 = 0, dc_M_bio_A1 = 0, dc_M_bio_A2 = 0,
  dc_Death_total = 0, dc_Death_B1 = 0, dc_Death_B2 = 0, dc_Death_A1 = 0, dc_Death_A2 = 0,
  Cons_POC_total = 0
)

# ------------------------------------------------------------------------------
# 3. ODE Model Formulation
# ------------------------------------------------------------------------------
model <- function(time, y, p) {
  with(as.list(c(y, p)), {
    
    # Monod equation for substrate uptake
    monod <- function(B, P, V, K) V * B * P / (P + K + 1e-12)
    
    # Maintenance gating: pool-specific availability 
    Theta_M_B1 <- 1 / (exp((-POC1 + K_M_B1) / (st_M_B1 * K_M_B1)) + 1)
    Theta_M_B2 <- 1 / (exp((-POC2 + K_M_B2) / (st_M_B2 * K_M_B2)) + 1)
    Theta_M_A1 <- 1 / (exp((-POC1 + K_M_A1) / (st_M_A1 * K_M_A1)) + 1)
    Theta_M_A2 <- 1 / (exp((-POC2 + K_M_A2) / (st_M_A2 * K_M_A2)) + 1)
    
    # Uptake (strict diagonal allocation; leak = 0)
    leak <- 0
    UB1 <- monod(B1, POC1 + leak * POC2, V_max_B1, K_V_B1_L)
    UA1 <- monod(A1, POC1 + leak * POC2, V_max_A1, K_V_A1_L)
    UB2 <- monod(B2, POC2 + leak * POC1, V_max_B2, K_V_B2_R)
    UA2 <- monod(A2, POC2 + leak * POC1, V_max_A2, K_V_A2_R)
    
    # Growth
    GB1 <- YG_B1 * UB1; GB2 <- YG_B2 * UB2
    GA1 <- YG_A1 * UA1; GA2 <- YG_A2 * UA2
    
    # Maintenance (Mx = exogenous from POC; Mn = endogenous from biomass)
    MxB1 <- mq_B1 * B1 * Theta_M_B1; MnB1 <- mq_B1 * B1 * (1 - Theta_M_B1)
    MxB2 <- mq_B2 * B2 * Theta_M_B2; MnB2 <- mq_B2 * B2 * (1 - Theta_M_B2)
    MxA1 <- mq_A1 * A1 * Theta_M_A1; MnA1 <- mq_A1 * A1 * (1 - Theta_M_A1)
    MxA2 <- mq_A2 * A2 * Theta_M_A2; MnA2 <- mq_A2 * A2 * (1 - Theta_M_A2)
    
    # Death
    DB1 <- alpha_B1 * B1; DB2 <- alpha_B2 * B2
    DA1 <- alpha_A1 * A1; DA2 <- alpha_A2 * A2
    Dtot <- DB1 + DB2 + DA1 + DA2
    
    # Necromass recycling to POC pools
    D_to_POC <- death_to_POC * Dtot
    D_to_POC1 <- death_to_POC1_frac * D_to_POC
    D_to_POC2 <- (1 - death_to_POC1_frac) * D_to_POC
    
    # State Derivatives (Substrate & Biomass balances)
    dPOC1 <- D_to_POC1 - (UB1 + UA1) - (MxB1 + MxA1)
    dPOC2 <- D_to_POC2 - (UB2 + UA2) - (MxB2 + MxA2)
    
    dB1 <- GB1 - MnB1 - DB1
    dB2 <- GB2 - MnB2 - DB2
    dA1 <- GA1 - MnA1 - DA1
    dA2 <- GA2 - MnA2 - DA2
    
    # Flux Accounting & Cumulative trackers
    Uptake_total <- UB1 + UB2 + UA1 + UA2
    derivs <- c(
      B1=dB1, B2=dB2, A1=dA1, A2=dA2,
      POC1=dPOC1, POC2=dPOC2,
      dc_POC_growth_total=Uptake_total,
      dc_POC_Growth_B1=UB1, dc_POC_Growth_B2=UB2, dc_POC_Growth_A1=UA1, dc_POC_Growth_A2=UA2,
      dc_M_total_POC=(MxB1+MxB2+MxA1+MxA2),
      dc_M_POC_B1=MxB1, dc_M_POC_B2=MxB2, dc_M_POC_A1=MxA1, dc_M_POC_A2=MxA2,
      dc_M_total_bio=(MnB1+MnB2+MnA1+MnA2),
      dc_M_bio_B1=MnB1, dc_M_bio_B2=MnB2, dc_M_bio_A1=MnA1, dc_M_bio_A2=MnA2,
      dc_Death_total=Dtot,
      dc_Death_B1=DB1, dc_Death_B2=DB2, dc_Death_A1=DA1, dc_Death_A2=DA2,
      Cons_POC_total=(Uptake_total + MxB1+MxB2+MxA1+MxA2)
    )
    list(derivs)
  })
}

# ------------------------------------------------------------------------------
# 4. Model Execution
# ------------------------------------------------------------------------------
times <- seq(0, 1000, length.out = 500)

out <- ode(
  y = start, times = times, func = model, parms = parms,
  method = "lsoda", rtol = 1e-8, atol = 1e-10
)

pred <- as_tibble(out) %>%
  mutate(
    time = as.numeric(time),
    across(c(B1, B2, A1, A2), as.numeric)
  )

write_csv(pred, file = "pred.csv")

# ------------------------------------------------------------------------------
# 5. Data Post-Processing (Flux & Fate Calculations)
# ------------------------------------------------------------------------------

# 5.1 Cumulative POC Budget Analysis
pred_POC_fate <- pred %>%
  mutate(
    cum_POC1_growth = dc_POC_Growth_B1 + dc_POC_Growth_A1,
    cum_POC1_maint  = dc_M_POC_B1 + dc_M_POC_A1,
    cum_POC1_necro  = parms$death_to_POC * parms$death_to_POC1_frac * dc_Death_total,
    tot_cum_POC1    = cum_POC1_growth + cum_POC1_maint - cum_POC1_necro,
    
    cum_POC2_growth = dc_POC_Growth_B2 + dc_POC_Growth_A2,
    cum_POC2_maint  = dc_M_POC_B2      + dc_M_POC_A2,
    cum_POC2_necro  = parms$death_to_POC * (1 - parms$death_to_POC1_frac) * dc_Death_total,
    tot_cum_POC2    = cum_POC2_growth + cum_POC2_maint - cum_POC2_necro,
    
    cum_POC_growth_total = cum_POC1_growth + cum_POC2_growth,
    cum_POC_maint_total  = cum_POC1_maint  + cum_POC2_maint,
    cum_POC_necro_total  = cum_POC1_necro  + cum_POC2_necro,
    tot_cum_POC          = cum_POC_growth_total + cum_POC_maint_total - cum_POC_necro_total
  )

# 5.2 Instantaneous Flux Calculations
pred_flux <- pred %>% arrange(time)
cum_cols <- grep("^dc_", names(pred_flux), value = TRUE)
flux_mat <- sapply(cum_cols, function(col) {
  vals <- pred_flux[[col]]
  c(NA, diff(vals) / diff(pred_flux$time)) 
})
colnames(flux_mat) <- paste0("flux_", cum_cols)
pred_flux <- bind_cols(pred_flux, as_tibble(flux_mat))

pred_POC_flux <- pred_flux %>%
  mutate(
    flux_POC1_growth = flux_dc_POC_Growth_B1 + flux_dc_POC_Growth_A1,
    flux_POC1_maint  = flux_dc_M_POC_B1      + flux_dc_M_POC_A1,
    flux_POC1_necro  = parms$death_to_POC * parms$death_to_POC1_frac * flux_dc_Death_total,
    
    flux_POC2_growth = flux_dc_POC_Growth_B2 + flux_dc_POC_Growth_A2,
    flux_POC2_maint  = flux_dc_M_POC_B2      + flux_dc_M_POC_A2,
    flux_POC2_necro  = parms$death_to_POC * (1 - parms$death_to_POC1_frac) * flux_dc_Death_total,
    
    flux_POC_growth_total = flux_POC1_growth + flux_POC2_growth,
    flux_POC_maint_total  = flux_POC1_maint  + flux_POC2_maint,
    flux_POC_necro_total  = flux_POC1_necro  + flux_POC2_necro,
    
    flux_POC1_net   = (flux_POC1_growth + flux_POC1_maint) - flux_POC1_necro,
    flux_POC2_net   = (flux_POC2_growth + flux_POC2_maint) - flux_POC2_necro,
    flux_POCtot_net = (flux_POC_growth_total + flux_POC_maint_total) - flux_POC_necro_total
  )

# Helper functions for fraction transformation
GUILD_LEVELS <- c("B1","B2","A1","A2", "No flux", "No activity")
guild_level <- c("B1", "B2", "A1", "A2")

make_POC_fraction_df <- function(dat, cols_vec) {
  dat %>% select(time, all_of(cols_vec)) %>%
    pivot_longer(cols = -time, names_to = "group", values_to = "value") %>%
    mutate(value = pmax(value, 0)) %>%
    group_by(time) %>%
    mutate(total = sum(value), frac = ifelse(total > 0, value / total, 0)) %>%
    ungroup()
}

make_flux_fraction_df <- function(dat, cols_vec) {
  long_df <- dat %>% select(time, all_of(cols_vec)) %>%
    pivot_longer(cols = -time, names_to = "Guild", values_to = "value") %>%
    mutate(value = pmax(value, 0))
  total_df <- long_df %>% group_by(time) %>% summarise(total = sum(value), .groups = "drop")
  
  df_nonzero <- long_df %>% inner_join(total_df, by = "time") %>%
    filter(total > 0) %>% mutate(frac = value / total)
  
  df_noflux <- total_df %>% filter(total == 0) %>%
    transmute(time = time, Guild = "No flux", frac = 1)
  
  bind_rows(df_nonzero %>% select(time, Guild, frac), df_noflux)
}

make_pathway_fraction_df <- function(dat, cols_vec) {
  long_df <- dat %>% select(time, all_of(cols_vec)) %>%
    pivot_longer(cols = -time, names_to = "Guild", values_to = "value") %>%
    mutate(value = pmax(as.numeric(ifelse(is.na(value), 0, value)), 0))
  
  total_df <- long_df %>% group_by(time) %>% summarise(total = sum(value), .groups = "drop")
  df_nonzero <- long_df %>% inner_join(total_df, by = "time") %>% filter(total > 0) %>% mutate(frac = value / total)
  df_no <- total_df %>% filter(total == 0) %>% transmute(time = time, Guild = "No activity", frac = 1)
  
  bind_rows(df_nonzero %>% select(time, Guild, frac), df_no)
}

# Data subsets for visualizations
POC_tot_group <- pred_POC_fate %>%
  mutate(
    cum_tot_B1 = dc_POC_Growth_B1 + dc_M_POC_B1,
    cum_tot_B2 = dc_POC_Growth_B2 + dc_M_POC_B2,
    cum_tot_A1 = dc_POC_Growth_A1 + dc_M_POC_A1,
    cum_tot_A2 = dc_POC_Growth_A2 + dc_M_POC_A2
  ) %>%
  make_POC_fraction_df(c("cum_tot_B1", "cum_tot_B2", "cum_tot_A1", "cum_tot_A2")) %>%
  mutate(
    group = recode(group, cum_tot_B1="B1", cum_tot_B2="B2", cum_tot_A1="A1", cum_tot_A2="A2"),
    group = factor(group, levels = guild_level)
  )

POC2_flux_group <- pred_flux %>%
  transmute(time, B2 = flux_dc_POC_Growth_B2 + flux_dc_M_POC_B2, A2 = flux_dc_POC_Growth_A2 + flux_dc_M_POC_A2) %>%
  make_flux_fraction_df(c("B2", "A2")) %>% mutate(Guild = factor(Guild, levels = GUILD_LEVELS))

exo_flux <- pred_flux %>%
  transmute(time, B1 = flux_dc_M_POC_B1, B2 = flux_dc_M_POC_B2, A1 = flux_dc_M_POC_A1, A2 = flux_dc_M_POC_A2) %>%
  make_pathway_fraction_df(c("B1","B2","A1","A2")) %>% mutate(Guild = factor(Guild, levels = GUILD_LEVELS))

death_flux <- pred_flux %>%
  transmute(time, B1 = flux_dc_Death_B1, B2 = flux_dc_Death_B2, A1 = flux_dc_Death_A1, A2 = flux_dc_Death_A2) %>%
  make_pathway_fraction_df(c("B1","B2","A1","A2")) %>% mutate(Guild = factor(Guild, levels = GUILD_LEVELS))


# ------------------------------------------------------------------------------
# 6. Core Visualizations (Main & Supplementary Panels)
# ------------------------------------------------------------------------------
# Shared plotting functions
plot_pathway_fraction <- function(df, title, y_lab = "Percentage (%)") {
  ggplot(df, aes(x = time, y = frac, fill = Guild)) +
    geom_area(position = "stack", alpha = 0.9, color = "white", size = 0.5) +
    coord_flip() + scale_x_reverse(expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = NULL), expand = c(0, 0)) +
    scale_fill_atlassian() +
    labs(x = "Sedimentary age (yr)", y = y_lab, title = title, fill = NULL) +
    theme_pubr(border = TRUE) + theme(plot.title = element_text(hjust = 0.5), legend.position = "top")
}

plot_POC_flux_guild <- function(df, title, y_lab = "Carbon flux (%)") {
  ggplot(df, aes(x = time, y = frac, fill = Guild)) +
    geom_area(position = "stack", alpha = 0.9) + coord_flip() + scale_x_reverse(expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    scale_fill_atlassian() + labs(x = "Sedimentary age (yr)", y = y_lab, title = title, fill = NULL) +
    theme_pubr(border = TRUE) + theme(plot.title = element_text(hjust = 0.5), legend.position = "top")
}

# Build Individual Panels
POC_budget.p <- ggplot(pred, aes(y = time)) +
  geom_ribbon(aes(xmin = 0, xmax = POC2, fill = "POC2")) +
  geom_ribbon(aes(xmin = POC2, xmax = POC1 + POC2, fill = "POC1")) +
  geom_path(aes(x = POC2), linewidth = 1, linetype = 2) +
  geom_path(aes(x = POC1 + POC2), linewidth = 1, color = "#dc4527", linetype = 2) +
  geom_point(data = observed.POC.df, aes(x = observed_POC, y = time), inherit.aes = FALSE, shape = 21, size = 4, alpha = 0.8, fill = "#182268", color = "white") +
  scale_y_reverse(limits = c(1000, 0), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 9000), expand = c(0, 0), name = "µgC/cm³") +
  scale_fill_manual(values = c("POC2" = "#667156", "POC1" = "#b7bca5"), name = NULL) +
  labs(x = "µgC/cm³", y = "Sedimentary age (yr)", title = "") +
  theme_pubr(border = TRUE) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

Bac.p <- ggplot(pred) +
  geom_path(aes(x = B1 + B2, y = time), color = "blue", linewidth = 1.5) +
  geom_point(data = observed.biomass.df, aes(x = Bac_cell, y = time), shape = 21, color = "white", size = 5, fill = "#8CA9FF") +
  scale_y_reverse(limits = c(1000, 0), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 2.6), expand = c(0, 0), name = "µgC/cm³") +
  labs(y = "Sedimentary age (yr)", x = " µg C/cm³", title = "") +
  theme_pubr(border = TRUE)

Arc.p <- ggplot(pred) +
  geom_path(aes(x = A1 + A2, y = time), color = "firebrick", linewidth = 1.5) +
  geom_point(data = observed.biomass.df, aes(x = Arc_cell, y = time), shape = 21, color = "white", fill = "#f48d8a", size = 5) +
  scale_y_reverse(limits = c(1000, 0), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 0.75), expand = c(0, 0), name = "µgC/cm³") +
  labs(y = "Sedimentary age (yr)", x = " µg C/cm³", title = "") +
  theme_pubr(border = TRUE)

Bathy.p <- ggplot(pred, aes(y = time)) +
  geom_ribbon(aes(xmin = 0, xmax = A2), fill = "#f2cd40", alpha = 0.8) +
  geom_path(aes(x = A2, y = time), linewidth = 1, linetype = 2, color = "#e00035") +
  geom_point(data = observed.biomass.df, aes(x = Bathy_cell, y = time), shape = 21, size = 5, alpha = 0.8, fill = "#182268", color = "white") +
  scale_y_reverse(limits = c(1000, 0), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 0.5), expand = c(0, 0), name = "µgC/cm³") +
  labs(y = "Sedimentary age (yr)", x = " µg C/cm³", title = "") +
  theme_pubr(border = TRUE)

p_exo_flux <- plot_pathway_fraction(exo_flux, "Exogenous maintenance")
p_death_flux <- plot_pathway_fraction(death_flux, "Death")
p_POC2_flux_guild <- plot_POC_flux_guild(POC2_flux_group, title = "")

POC_budget_guild.p <- ggplot(POC_tot_group, aes(x = time, y = frac, fill = group)) +
  geom_area(position = "stack", alpha = 0.85) +
  coord_flip() + scale_x_reverse(expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = ""), expand = c(0, 0)) +
  scale_fill_atlassian() +
  labs(x = "Sedimentary age (yr)", y = "Carbon budget", title = "Carbon budget") +
  theme_pubr(border = TRUE) + theme(plot.title = element_text(hjust = 0.5), legend.position = "top")

# ==============================================================================
# 7. Final Manuscript Output Figure Array
# ==============================================================================
final_figure <- ggarrange(
  ncol = 4, nrow = 2, 
  POC_budget.p, Bac.p, Arc.p, Bathy.p, 
  p_exo_flux, p_death_flux, POC_budget_guild.p, p_POC2_flux_guild,
  legend = "none", labels = "auto"
)

# Render the assembled main figure
print(final_figure)
