# ==============================================================================
# ===== Galleria Staphylococcus Phage: Phylogenetic MCMCglmms ==================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# This script analyses experimental Galleria mellonella infection data, using 
# non-linear least squares models to summarise variation in the non-linear
# dynamics to a small number of parameters, then analysing variation in these
# parameters using phylogenetic MCMCglmms.

# In vitro data from Walsh et al., (2023) The host phylogeny determines viral
# infectivity and replication across Staphylococcus host species.

# ----- 0.2. Dependencies ------------------------------------------------------

library(tidyverse)
library(MCMCglmm)
library(ape)
library(minpack.lm)
library(here)
library(viridisLite)
library(scales)
library(patchwork)
library(ggtree)

here::i_am("Scripts/01_models.R") # MacOS dynamic path fix

# ----- 0.3. Load Data ---------------------------------------------------------

data <- read_csv(here("Data", "Galleria_main.csv"))

data_weights <- read_csv(here("data", "Galleria_weights.csv"))

tree <- read.tree(here("data", "Bacteria_phylogeny.nwk"))

data_invitro <- read_csv(here("data", "In_vitro.csv"))

# ------------------------------------------------------------------------------
# ----- 1. Data Wrangling ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Format Data -------------------------------------------------------

data <- filter(data, !is.na(Bacteria)) # Remove controls

data$ID <- factor(data$ID)
data$Condition <- factor(data$Condition)
data$Bacteria <- factor(data$Bacteria)
data$animal <- data$Bacteria
data$Plate <- factor(data$Plate)
data$Well <- factor(data$Well)

# Infer weights from pixels
model_weight <- lm(Weight ~ Area, data = data_weights)

data$Weight <- predict(model_weight, newdata = data)

# ----- 1.2. Aggregate by plate ------------------------------------------------

data_plate <- data %>% group_by(Bacteria, Condition, Time, Block) %>%
  filter(!is.na(Melanisation)) %>%
  summarise(u_Alive = mean(Alive),
            u_Melanisation = mean(Melanisation),
            u_Weight = mean(Weight),
            n = n())

data_plate$u_Dead <- 1 - data_plate$u_Alive

data_weights <- data_plate %>% ungroup %>% group_by(Bacteria, Condition, Block) %>%
  summarise(u_Weight = mean(u_Weight))

# ----- 1.3. Fit Logistic Curves to Mortality ----------------------------------

mortality_24h <- data_plate %>% ungroup() %>%
  filter(Time == 24) %>%
  dplyr::select(Bacteria, Condition, Block, U24 = u_Dead)

parameters_mortality_cure <- data_plate %>% ungroup() %>%
  left_join(mortality_24h, by = c("Bacteria", "Condition", "Block")) %>%
  group_by(Bacteria, Condition, Block) %>%
  summarise(fit = list(tryCatch(
    nlsLM(u_Dead ~ U24 / (1 + exp(-k * (Time - t50))),
          start = list(t50 = median(Time), k = 0.5),
          lower = c(t50 = 0, k = 0),
          upper = c(t50 = 48, k = 1),
          control = nls.lm.control(maxiter = 1024)),
    error = function(e) NULL)),
    U = unique(U24),
    .groups = "drop") %>%
  mutate(coef = map(fit, function(x) {
    if (!is.null(x)) coef(x) else c(t50 = NA_real_, k = NA_real_)
  })) %>%
  unnest_wider(coef) %>%
  select(-fit)

parameters_mortality_unobserved <- data_plate %>% ungroup() %>%
  group_by(Bacteria, Condition, Block) %>%
  summarise(fit = list(tryCatch(
    nlsLM(u_Dead ~ U / (1 + exp(-k * (Time - t50))),
          start = list(U = max(u_Dead, na.rm = TRUE), t50 = median(Time), k = 0.5),
          lower = c(U = 0, t50 = 0, k = 0),
          upper = c(U = 1, t50 = 48, k = 1),
          control = nls.lm.control(maxiter = 1024)),
    error = function(e) NULL)),
    .groups = "drop") %>%
  mutate(coef = map(fit, function(x) {
    if (!is.null(x)) coef(x) else c(U = NA_real_, t50 = NA_real_, k = NA_real_)
  })) %>%
  unnest_wider(coef) %>%
  select(-fit)


# ----- 1.4. Fitted Plots (Mortality Curves - Cured) ---------------------------

time_grid <- seq(min(data_plate$Time),
                 max(data_plate$Time),
                 length.out = 200)

predicted_mortality_cure <- merge(parameters_mortality_cure, expand.grid(Time = time_grid),all = TRUE)

predicted_mortality_cure$u_Dead <- with(predicted_mortality_cure, U / (1 + exp(-k * (Time - t50))))

predicted_mortality_cure$u_Dead <- ifelse(is.na(predicted_mortality_cure$u_Dead), 0, predicted_mortality_cure$u_Dead)

ggplot(data_plate, aes(Time, u_Dead, colour = Condition)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_line(data = predicted_mortality_cure, aes(group = interaction(Condition, Block)), alpha = 0.5) +
  facet_wrap(~Bacteria) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(aspect.ratio = 1/2,
        panel.grid = element_blank())

predicted_mortality_cure_pointwise <- left_join(data_plate, parameters_mortality_cure) %>%
  mutate(pred_u_Dead = U / (1 + exp(-k * (Time - t50))))

predicted_mortality_cure_pointwise$pred_u_Dead <- ifelse(is.na(predicted_mortality_cure_pointwise$pred_u_Dead), 0, predicted_mortality_cure_pointwise$pred_u_Dead)

ggplot(predicted_mortality_cure_pointwise, aes(pred_u_Dead, u_Dead)) +
  geom_point(alpha = 0.05) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw() +
  coord_fixed()

# ----- 1.5. Fitted Plots (Mortality Curves - Unobserved) ----------------------

predicted_mortality_unobserved <- merge(parameters_mortality_unobserved, expand.grid(Time = time_grid),all = TRUE)

predicted_mortality_unobserved$u_Dead <- with(predicted_mortality_unobserved, U / (1 + exp(-k * (Time - t50))))

predicted_mortality_unobserved$u_Dead <- ifelse(is.na(predicted_mortality_unobserved$u_Dead), 0, predicted_mortality_unobserved$u_Dead)

ggplot(data_plate, aes(Time, u_Dead, colour = Condition)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_line(data = predicted_mortality_unobserved, aes(group = interaction(Condition, Block)), alpha = 0.5) +
  facet_wrap(~Bacteria) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(aspect.ratio = 1/2,
        panel.grid = element_blank())

predicted_mortality_unobserved_pointwise <- left_join(data_plate, parameters_mortality_unobserved) %>%
  mutate(pred_u_Dead = U / (1 + exp(-k * (Time - t50))))

predicted_mortality_unobserved_pointwise$pred_u_Dead <- ifelse(is.na(predicted_mortality_unobserved_pointwise$pred_u_Dead), 0, predicted_mortality_unobserved_pointwise$pred_u_Dead)

ggplot(predicted_mortality_unobserved_pointwise, aes(pred_u_Dead, u_Dead)) +
  geom_point(alpha = 0.05) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw() +
  coord_fixed()

# ----- 1.6. Fit Michaelis-Menten Curves to Melanisation -----------------------

parameters_melanisation <- data_plate %>% ungroup() %>%
  filter(!Time %in% c(0, 24)) %>%
  group_by(Bacteria, Condition, Block) %>%
  summarise(fit = list(tryCatch(
    nlsLM(u_Melanisation ~ Vmax * Time / (K + Time),
          start = list(Vmax = 1, K = 4),
          lower = c(Vmax = 0, K = 0),
          upper = c(Vmax = 1.2, K = 30),
          control = nls.lm.control(maxiter = 1024)),
    error = function(e) NULL)),
    .groups = "drop") %>%
  mutate(coef = map(fit, function(x) {
    if (!is.null(x)) coef(x) else c(Vmax = NA_real_, K = NA_real_)
  })) %>%
  unnest_wider(coef) %>%
  dplyr::select(-fit)

# ----- 1.8 Fitted Plots (Melanisation Curves) ---------------------------------

predicted_melanisation <- merge(parameters_melanisation, expand.grid(Time = time_grid),all = TRUE)

predicted_melanisation$u_Melanisation <- with(predicted_melanisation, Vmax * Time / (K + Time))

ggplot(data_plate, aes(Time, u_Melanisation, colour = Condition)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_line(data = predicted_melanisation, aes(group = interaction(Condition, Block)), linewidth = 1) +
  facet_wrap(~Bacteria) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(aspect.ratio = 1/2,
        panel.grid = element_blank())

predicted_melanisation_pointwise <- left_join(data_plate, parameters_melanisation) %>%
  mutate(pred_u_Melanisation = Vmax * Time / (K + Time))

ggplot(filter(predicted_melanisation_pointwise, Time != 0), aes(pred_u_Melanisation, u_Melanisation)) +
  geom_point(alpha = 0.05) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw() +
  coord_fixed()

# ----- 1.9 Join Datasets ------------------------------------------------------

colnames(parameters_mortality_unobserved)[4:6] <- c("mor_unobserved_U", "mor_unobserved_t50", "mor_unobserved_k")
colnames(parameters_mortality_cure)[4:6] <- c("mor_cure_U", "mor_cure_t50", "mor_cure_k")
colnames(parameters_melanisation)[4:5] <- c("mel_Vmax", "mel_K")

parameters_all <- parameters_mortality_unobserved %>%
  left_join(parameters_mortality_cure) %>%
  left_join(parameters_melanisation) %>%
  left_join(data_weights)

# ----- 1.10 Create Phage Effect (Delta) Dataset -------------------------------

parameters_all_nophage <- filter(parameters_all, Condition == "NoPhage")
parameters_all_phage <- filter(parameters_all, Condition == "Phage")

parameters_all_delta <- parameters_all_nophage[,1:3]

for (i in c(4:12)){
  parameters_all_delta$placeholder <- parameters_all_phage[[i]] - parameters_all_nophage[[i]]
  colnames(parameters_all_delta)[i] <- colnames(parameters_all_phage)[i]
}

parameters_all_delta$Condition <- "Delta"

# ----- 1.11. Recenter Weights -------------------------------------------------

parameters_all$u_Weight <- scale(parameters_all$u_Weight, scale = FALSE)
parameters_all_delta$u_Weight <- scale(parameters_all_delta$u_Weight, scale = FALSE)

# ----- 1.12. Wrangle in-vitro OD data ------------------------------------------

data_invitro$animal <- ifelse(data_invitro$animal == "8325-4", "83254", data_invitro$animal)

# ------------------------------------------------------------------------------
# ----- w. MCMCglmms -----------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. Run Iterations ----------------------------------------------------

itt <- 10 # 10 for quick exploration, 100 for publication

# ----- 3.2. Priors ------------------------------------------------------------

prior.uni.rep <- list(
  R = list(V = diag(1), nu = 0.02),
  G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000),
           G2 = list(V = diag(1), nu = 0.02)))

prior.uni.her <- list(
  R = list(V = diag(1), nu = 0.02),
  G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000),
           G2 = list(V = diag(1), nu = 0.02),
           G3 = list(V = diag(1), nu = 0.02)))

prior.bi.rep <- list(
  R = list(V = diag(2), nu = 0.02),
  G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
           G3 = list(V = diag(1), nu = 0.02)))


prior.bi.her <- list(
  R = list(V = diag(2), nu = 0.02),
  G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
           G2 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
           G3 = list(V = diag(1), nu = 0.02),
           G4 = list(V = diag(1), nu = 0.02)))

# ----- 3.3. MCMCglmms (Mortality) ---------------------------------------------

parameters_all$animal <- parameters_all$Bacteria
parameters_all_delta$animal <- parameters_all_delta$Bacteria

parameters_all <- as.data.frame(parameters_all)
parameters_all_delta <- as.data.frame(parameters_all_delta)

estimates_fixed <- data.frame(condition = character(),
                              assumption = character(),
                              parameter = character(),
                              mean = numeric(),
                              low = numeric(),
                              high = numeric())

estimates_rep <- data.frame(condition = character(),
                            assumption = character(),
                            parameter = character(),
                            mean = numeric(),
                            low = numeric(),
                            high = numeric())

estimates_her <- data.frame(condition = character(),
                            assumption = character(),
                            parameter = character(),
                            mean = numeric(),
                            low = numeric(),
                            high = numeric())

assumptions <- c("unobserved", "cure")
parameters <- c("U", "t50", "k")

for (assumption in assumptions){
  for (parameter in parameters){
    
    file <- sprintf("fixed_%s_%s.Rdata", assumption, parameter)
    
    par <- sprintf("mor_%s_%s", assumption, parameter)
    
    form <- as.formula(paste(par, "~ Condition + u_weight"))
    
    if (!file.exists(here("models", file))) {
      
      model <- MCMCglmm(fixed = form, random = ~ animal + Block,
                        rcov = ~units, family = "gaussian", data = parameters_all, pedigree = tree,
                        prior = prior.uni.rep, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
      
      save(model, file = here("models", file))
      
    } else {load(here("models", file))}
    
    fixed_nophage <- model$Sol[,1]
    fixed_phage <- model$Sol[,1] + model$Sol[,2]
    fixed_weight <- model$Sol[,3]
    fixed_contrast <- model$Sol[,2]
    
    estimates_fixed <- rbind(estimates_fixed,
                             data.frame(condition = "NoPhage",
                                        assumption = assumption,
                                        parameter = parameter,
                                        mean = round(mean(fixed_nophage), 4),
                                        low = round(HPDinterval(fixed_nophage)[1], 4),
                                        high = round(HPDinterval(fixed_nophage)[2], 4)))
    
    estimates_fixed <- rbind(estimates_fixed,
                             data.frame(condition = "Phage",
                                        assumption = assumption,
                                        parameter = parameter,
                                        mean = round(mean(fixed_phage), 4),
                                        low = round(HPDinterval(fixed_phage)[1], 4),
                                        high = round(HPDinterval(fixed_phage)[2], 4)))
    
    estimates_fixed <- rbind(estimates_fixed,
                             data.frame(condition = "Contrast",
                                        assumption = assumption,
                                        parameter = parameter,
                                        mean = round(mean(fixed_contrast), 4),
                                        low = round(HPDinterval(fixed_contrast)[1], 4),
                                        high = round(HPDinterval(fixed_contrast)[2], 4)))
    
    estimates_fixed <- rbind(estimates_fixed,
                             data.frame(condition = "Weight",
                                        assumption = assumption,
                                        parameter = parameter,
                                        mean = round(mean(fixed_weight), 4),
                                        low = round(HPDinterval(fixed_weight)[1], 4),
                                        high = round(HPDinterval(fixed_weight)[2], 4)))
    
    file <- sprintf("rep_mor_%s_%s_%s.Rdata", "bivariate", assumption, parameter)
    
    par <- sprintf("mor_%s_%s", assumption, parameter)
    
    form <- as.formula(paste(par, "~ Condition + u_weight"))
    
    if (!file.exists(here("models", file))) {
      
      model <- MCMCglmm(fixed = form, random = ~ us(Condition):animal + Condition:Block,
                        rcov = ~idh(Condition):units, family = "gaussian", data = parameters_all, pedigree = tree,
                        prior = prior.bi.rep, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
      
      save(model, file = here("models", file))
      
    } else {load(here("models", file))}
    
    rep_nophage <- model$VCV[,1] / (model$VCV[,1] + model$VCV[,6])
    rep_phage <- model$VCV[,4] / (model$VCV[,4] + model$VCV[,7])
    
    estimates_rep <- rbind(estimates_rep,
                           data.frame(condition = "NoPhage",
                                      assumption = assumption,
                                      parameter = parameter,
                                      mean = round(mean(rep_nophage), 4),
                                      low = round(HPDinterval(rep_nophage)[1], 4),
                                      high = round(HPDinterval(rep_nophage)[2], 4)))
    
    estimates_rep <- rbind(estimates_rep,
                           data.frame(condition = "Phage",
                                      assumption = assumption,
                                      parameter = parameter,
                                      mean = round(mean(rep_phage), 4),
                                      low = round(HPDinterval(rep_phage)[1], 4),
                                      high = round(HPDinterval(rep_phage)[2], 4)))
    
    file <- sprintf("her_mor_%s_%s_%s.Rdata", "bivariate", assumption, parameter)
    
    if (!file.exists(here("models", file))) {
      
      model <- MCMCglmm(fixed = form, random = ~ us(Condition):animal + idh(Condition):Bacteria + Bacteria:Block + Condition:Block,
                        rcov = ~idh(Condition):units, family = "gaussian", data = parameters_all, pedigree = tree,
                        prior = prior.bi.her, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
      
      save(model, file = here("models", file))
      
    } else {load(here("models", file))}
    
    her_nophage <- model$VCV[,1] / (model$VCV[,1] + model$VCV[,5])
    her_phage <- model$VCV[,4] / (model$VCV[,4] + model$VCV[,6])
    
    estimates_her <- rbind(estimates_her,
                           data.frame(condition = "NoPhage",
                                      assumption = assumption,
                                      parameter = parameter,
                                      mean = round(mean(her_nophage), 4),
                                      low = round(HPDinterval(her_nophage)[1], 4),
                                      high = round(HPDinterval(her_nophage)[2], 4)))
    
    estimates_her <- rbind(estimates_her,
                           data.frame(condition = "Phage",
                                      assumption = assumption,
                                      parameter = parameter,
                                      mean = round(mean(her_phage), 4),
                                      low = round(HPDinterval(her_phage)[1], 4),
                                      high = round(HPDinterval(her_phage)[2], 4)))
    
  }
}

# ----- 3.4. MCMCglmms (Melanisation) ------------------------------------------

parameters <- c("Vmax", "K")

for (parameter in parameters){
  
  file <- sprintf("fixed_%s.Rdata", parameter)
  
  par <- sprintf("mel_%s", parameter)
  
  form <- as.formula(paste(par, "~ Condition + u_weight"))
  
  if (!file.exists(here("models", file))) {
    
    model <- MCMCglmm(fixed = form, random = ~ animal + Block,
                      rcov = ~units, family = "gaussian", data = parameters_all, pedigree = tree,
                      prior = prior.uni.rep, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
    
    save(model, file = here("models", file))
    
  } else {load(here("models", file))}
  
  fixed_nophage <- model$Sol[,1]
  fixed_phage <- model$Sol[,1] + model$Sol[,2]
  fixed_weight <- model$Sol[,3]
  fixed_contrast <- model$Sol[,2]
  
  estimates_fixed <- rbind(estimates_fixed,
                           data.frame(condition = "NoPhage",
                                      assumption = NA,
                                      parameter = parameter,
                                      mean = round(mean(fixed_nophage), 4),
                                      low = round(HPDinterval(fixed_nophage)[1], 4),
                                      high = round(HPDinterval(fixed_nophage)[2], 4)))
  
  estimates_fixed <- rbind(estimates_fixed,
                           data.frame(condition = "Phage",
                                      assumption = NA,
                                      parameter = parameter,
                                      mean = round(mean(fixed_phage), 4),
                                      low = round(HPDinterval(fixed_phage)[1], 4),
                                      high = round(HPDinterval(fixed_phage)[2], 4)))
  
  estimates_fixed <- rbind(estimates_fixed,
                           data.frame(condition = "Contrast",
                                      assumption = NA,
                                      parameter = parameter,
                                      mean = round(mean(fixed_contrast), 4),
                                      low = round(HPDinterval(fixed_contrast)[1], 4),
                                      high = round(HPDinterval(fixed_contrast)[2], 4)))
  
  estimates_fixed <- rbind(estimates_fixed,
                           data.frame(condition = "Weight",
                                      assumption = NA,
                                      parameter = parameter,
                                      mean = round(mean(fixed_weight), 4),
                                      low = round(HPDinterval(fixed_weight)[1], 4),
                                      high = round(HPDinterval(fixed_weight)[2], 4)))
  
  file <- sprintf("rep_mel_%s_%s.Rdata", "bivariate", parameter)
  
  par <- sprintf("mel_%s", parameter)
  
  form <- as.formula(paste(par, "~ Condition + u_weight"))
  
  if (!file.exists(here("models", file))) {
    
    model <- MCMCglmm(fixed = form, random = ~ us(Condition):animal + Condition:Block,
                      rcov = ~idh(Condition):units, family = "gaussian", data = parameters_all, pedigree = tree,
                      prior = prior.bi.rep, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
    
    save(model, file = here("models", file))
    
  } else {load(here("models", file))}
  
  if (parameter == "Vmax"){
    rep_nophage <- model$VCV[,1] / (model$VCV[,1] + model$VCV[,4])
    rep_phage <- model$VCV[,2] / (model$VCV[,2] + model$VCV[,5])
  } else {
    rep_nophage <- model$VCV[,1] / (model$VCV[,1] + model$VCV[,6])
    rep_phage <- model$VCV[,4] / (model$VCV[,4] + model$VCV[,7])
  }
  
  
  estimates_rep <- rbind(estimates_rep,
                         data.frame(condition = "NoPhage",
                                    assumption = NA,
                                    parameter = parameter,
                                    mean = round(mean(rep_nophage), 4),
                                    low = round(HPDinterval(rep_nophage)[1], 4),
                                    high = round(HPDinterval(rep_nophage)[2], 4)))
  
  estimates_rep <- rbind(estimates_rep,
                         data.frame(condition = "Phage",
                                    assumption = NA,
                                    parameter = parameter,
                                    mean = round(mean(rep_phage), 4),
                                    low = round(HPDinterval(rep_phage)[1], 4),
                                    high = round(HPDinterval(rep_phage)[2], 4)))
  
  
  file <- sprintf("her_mel_%s_%s.Rdata", "bivariate", parameter)
  
  if (!file.exists(here("models", file))) {
    
    model <- MCMCglmm(fixed = form, random = ~ idh(Condition):animal + idh(Condition):Bacteria + Bacteria:Block + Condition:Block,
                      rcov = ~idh(Condition):units, family = "gaussian", data = parameters_all, pedigree = tree,
                      prior = prior.bi.her, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
    
    save(model, file = here("models", file))
    
  } else {load(here("models", file))}
  
  if (parameter == "K"){
    her_nophage <- model$VCV[,1] / (model$VCV[,1] + model$VCV[,3])
    her_phage <- model$VCV[,2] / (model$VCV[,2] + model$VCV[,4])
    
  } else {
    her_nophage <- model$VCV[,1] / (model$VCV[,1] + model$VCV[,5])
    her_phage <- model$VCV[,4] / (model$VCV[,4] + model$VCV[,6])
  }
  
  
  estimates_her <- rbind(estimates_her,
                         data.frame(condition = "NoPhage",
                                    assumption = NA,
                                    parameter = parameter,
                                    mean = round(mean(her_nophage), 4),
                                    low = round(HPDinterval(her_nophage)[1], 4),
                                    high = round(HPDinterval(her_nophage)[2], 4)))
  
  estimates_her <- rbind(estimates_her,
                         data.frame(condition = "Phage",
                                    assumption = NA,
                                    parameter = parameter,
                                    mean = round(mean(her_phage), 4),
                                    low = round(HPDinterval(her_phage)[1], 4),
                                    high = round(HPDinterval(her_phage)[2], 4)))
  
}

# ----- 3.5. MCMCglmms (Mortality Deltas) --------------------------------------

assumptions <- c("unobserved", "cure")
parameters <- c("U", "t50", "k")

prior.uni.rep <- list(
  R = list(V = diag(1), nu = 0.02),
  G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000)))

for (parameter in parameters){
  for (assumption in assumptions){
    
    file <- sprintf("fixedDelta_%s_%s.Rdata", assumption, parameter)
    
    par <- sprintf("mor_%s_%s", assumption, parameter)
    
    form <- as.formula(paste(par, "~ u_weight"))
    
    if (!file.exists(here("models", file))) {
      
      model <- MCMCglmm(fixed = form, random = ~ animal,
                        rcov = ~units, family = "gaussian", data = parameters_all_delta, pedigree = tree,
                        prior = prior.uni.rep, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
      
      save(model, file = here("models", file))
      
    } else {load(here("models", file))}
    
    fixed_delta <- model$Sol[,1]
    fixed_weight <- model$Sol[,2]
    
    estimates_fixed <- rbind(estimates_fixed,
                             data.frame(condition = "Delta",
                                        assumption = assumption,
                                        parameter = parameter,
                                        mean = round(mean(fixed_delta), 4),
                                        low = round(HPDinterval(fixed_delta)[1], 4),
                                        high = round(HPDinterval(fixed_delta)[2], 4)))
    
    estimates_fixed <- rbind(estimates_fixed,
                             data.frame(condition = "Delta Weight",
                                        assumption = assumption,
                                        parameter = parameter,
                                        mean = round(mean(fixed_weight), 4),
                                        low = round(HPDinterval(fixed_weight)[1], 4),
                                        high = round(HPDinterval(fixed_weight)[2], 4)))
    
    file <- sprintf("rep_mor_%s_%s_%s.Rdata", "bivariateDelta", assumption, parameter)
    
    par <- sprintf("mor_%s_%s", assumption, parameter)
    
    form <- as.formula(paste(par, "~ u_weight"))
    
    if (!file.exists(here("models", file))) {
      
      model <- MCMCglmm(fixed = form, random = ~ animal + Bacteria:Block,
                        rcov = ~units, family = "gaussian", data = parameters_all_delta, pedigree = tree,
                        prior = prior.uni.rep, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
      
      save(model, file = here("models", file))
      
    } else {load(here("models", file))}
    
    rep <- model$VCV[,1] / (model$VCV[,1] + model$VCV[,3])
    
    estimates_rep <- rbind(estimates_rep,
                           data.frame(condition = "Delta",
                                      assumption = assumption,
                                      parameter = parameter,
                                      mean = round(mean(rep), 4),
                                      low = round(HPDinterval(rep)[1], 4),
                                      high = round(HPDinterval(rep)[2], 4)))
    
    file <- sprintf("her_mor_%s_%s_%s.Rdata", "bivariateDelta", assumption, parameter)
    
    if (!file.exists(here("models", file))) {
      
      model <- MCMCglmm(fixed = form, random = ~ animal + Bacteria + Bacteria:Block,
                        rcov = ~units, family = "gaussian", data = parameters_all_delta, pedigree = tree,
                        prior = prior.uni.her, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
      
      save(model, file = here("models", file))
      
    } else {load(here("models", file))}
    
    her <- model$VCV[,1] / (model$VCV[,1] + model$VCV[,2])
    
    estimates_her <- rbind(estimates_her,
                           data.frame(condition = "Delta",
                                      assumption = assumption,
                                      parameter = parameter,
                                      mean = round(mean(her), 4),
                                      low = round(HPDinterval(her)[1], 4),
                                      high = round(HPDinterval(her)[2], 4)))
    
  }
}

# ----- 3.6. MCMCglmms (Melanisation Deltas) -----------------------------------

parameters <- c("Vmax", "K")

for (parameter in parameters){
  
  file <- sprintf("fixedDelta_%s.Rdata", parameter)
  
  par <- sprintf("mel_%s", parameter)
  
  form <- as.formula(paste(par, "~ u_weight"))
  
  if (!file.exists(here("models", file))) {
    
    model <- MCMCglmm(fixed = form, random = ~ animal,
                      rcov = ~units, family = "gaussian", data = parameters_all_delta, pedigree = tree,
                      prior = prior.uni.rep, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
    
    save(model, file = here("models", file))
    
  } else {load(here("models", file))}
  
  fixed_delta <- model$Sol[,1]
  fixed_weight <- model$Sol[,2]
  
  estimates_fixed <- rbind(estimates_fixed,
                           data.frame(condition = "Delta",
                                      assumption = assumption,
                                      parameter = parameter,
                                      mean = round(mean(fixed_delta), 4),
                                      low = round(HPDinterval(fixed_delta)[1], 4),
                                      high = round(HPDinterval(fixed_delta)[2], 4)))
  
  estimates_fixed <- rbind(estimates_fixed,
                           data.frame(condition = "Delta Weight",
                                      assumption = assumption,
                                      parameter = parameter,
                                      mean = round(mean(fixed_weight), 4),
                                      low = round(HPDinterval(fixed_weight)[1], 4),
                                      high = round(HPDinterval(fixed_weight)[2], 4)))
  
  
  file <- sprintf("rep_mel_%s_%s.Rdata", "bivariateDelta", parameter)
  
  par <- sprintf("mel_%s", parameter)
  
  form <- as.formula(paste(par, "~ u_weight"))
  
  if (!file.exists(here("models", file))) {
    
    model <- MCMCglmm(fixed = form, random = ~ animal + Bacteria:Block,
                      rcov = ~units, family = "gaussian", data = parameters_all_delta, pedigree = tree,
                      prior = prior.uni.rep, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
    
    save(model, file = here("models", file))
    
  } else {load(here("models", file))}
  
  rep <- model$VCV[,1] / (model$VCV[,1] + model$VCV[,3])
  
  estimates_rep <- rbind(estimates_rep,
                         data.frame(condition = "Delta",
                                    assumption = NA,
                                    parameter = parameter,
                                    mean = round(mean(rep), 4),
                                    low = round(HPDinterval(rep)[1], 4),
                                    high = round(HPDinterval(rep)[2], 4)))
  
  file <- sprintf("her_mel_%s_%s.Rdata", "bivariateDelta", parameter)
  
  if (!file.exists(here("models", file))) {
    
    model <- MCMCglmm(fixed = form, random = ~ animal + Bacteria + Bacteria:Block,
                      rcov = ~units, family = "gaussian", data = parameters_all_delta, pedigree = tree,
                      prior = prior.uni.her, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
    
    save(model, file = here("models", file))
    
  } else {load(here("models", file))}
  
  her <- model$VCV[,1] / (model$VCV[,1] + model$VCV[,2])
  
  estimates_her <- rbind(estimates_her,
                         data.frame(condition = "Delta",
                                    assumption = NA,
                                    parameter = parameter,
                                    mean = round(mean(her), 4),
                                    low = round(HPDinterval(her)[1], 4),
                                    high = round(HPDinterval(her)[2], 4)))
  
}

estimates_rep$sig <- estimates_rep$low > 0.05
estimates_her$sig <- estimates_her$low > 0.05

estimates_fixed$sig <- sign(estimates_fixed$low) == sign(estimates_fixed$high)

estimates_rep
estimates_her

estimates_fixed

# ------------------------------------------------------------------------------
# ----- 4. Correlations --------------------------------------------------------
# ------------------------------------------------------------------------------

prior.corr <- list(
  R = list(V = diag(2), nu = 2),
  G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
           G3 = list(V = diag(2), nu = 0.02)))

# ----- 4.1. Condition Correlations --------------------------------------------

corr_conditions <- data.frame(Condition = character(),
                              x = character(),
                              y = character(),
                              Fx_mean = numeric(),
                              Fx_low  = numeric(),
                              Fx_high = numeric(),
                              Fy_mean = numeric(),
                              Fy_low  = numeric(),
                              Fy_high = numeric(),
                              Bp_mean = numeric(),
                              Bp_low = numeric(),
                              Bp_high = numeric(),
                              Rp_mean = numeric(),
                              Rp_low = numeric(),
                              Rp_high = numeric())


parameters <- c("U", "t50", "Vmax", "K")

for (parameter in parameters){
  
  form <- as.formula(paste0("cbind(", "NoPhage", ", ", "Phage", ") ~ trait - 1"))
  
  par <- ifelse(parameter %in% c("U", "t50", "k"), sprintf("mor_cure_%s", parameter), sprintf("mel_%s", parameter))
  
  parameter_current <- dplyr::select(parameters_all, Bacteria, Condition, Block, par)
  
  parameter_phage <- filter(parameter_current, Condition == "Phage")
  parameter_nophage <- filter(parameter_current, Condition != "Phage")
  
  parameter_phage$Condition <- NULL
  parameter_nophage$Condition <- NULL
  
  colnames(parameter_phage)[3] <- "Phage"
  colnames(parameter_nophage)[3] <- "NoPhage"
  
  parameter_model <- left_join(parameter_phage, parameter_nophage)
  parameter_model$animal <- parameter_model$Bacteria
  parameter_model <- as.data.frame(parameter_model)
  
  file <- sprintf("cor_condition_%s.Rdata", parameter)
  
  if (!file.exists(here("models", file))) {
    
    model <- MCMCglmm(fixed = form, random = ~ us(trait):animal + idh(trait):Bacteria:Block,
                      rcov = ~idh(trait):units, family = c("gaussian", "gaussian"), data = parameter_model, pedigree = tree,
                      prior = prior.corr, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
    
    save(model, file = here("models", file))
    
  } else {load(here("models", file))}
  
  var_p_x <- model$VCV[,1]
  
  var_p_y <- model$VCV[,4]
  
  cov_p <- model$VCV[,2]
  
  Bp <- as.mcmc(cov_p / var_p_x)
  
  Rp <- as.mcmc(cov_p /(sqrt(var_p_x * var_p_y)))
  
  Fx <- as.mcmc(model$Sol[,1])
  Fy <- as.mcmc(model$Sol[,2])
  
  corr_conditions <- rbind(corr_conditions, data.frame(Condition = parameter,
                                                       x = "NoPhage",
                                                       y = "Phage",
                                                       
                                                       Fx_mean = mean(Fx),
                                                       Fx_low  = HPDinterval(Fx)[1],
                                                       Fx_high = HPDinterval(Fx)[2],
                                                       Fy_mean = mean(Fy),
                                                       Fy_low  = HPDinterval(Fy)[1],
                                                       Fy_high = HPDinterval(Fy)[2],
                                                       
                                                       Bp_mean = mean(Bp),
                                                       Bp_low = HPDinterval(Bp)[1],
                                                       Bp_high = HPDinterval(Bp)[2],
                                                       Rp_mean = mean(Rp),
                                                       Rp_low = HPDinterval(Rp)[1],
                                                       Rp_high = HPDinterval(Rp)[2]))
}

# ----- 4.2. Parameter Correlations --------------------------------------------

correlations <- data.frame(Condition = character(),
                           Assumption = character(),
                           x = character(),
                           y = character(),
                           Fx_mean = numeric(),
                           Fx_low  = numeric(),
                           Fx_high = numeric(),
                           Fy_mean = numeric(),
                           Fy_low  = numeric(),
                           Fy_high = numeric(),
                           Bp_mean = numeric(),
                           Bp_low = numeric(),
                           Bp_high = numeric(),
                           Rp_mean = numeric(),
                           Rp_low = numeric(),
                           Rp_high = numeric())

assumptions <- c("unobserved", "cure")

parameters <- c("U", "t50", "Vmax", "K")

for (assumption in assumptions){
  
  params <- sprintf("mor_%s_%s", assumption, parameters[1:2])
  params <- c(params, sprintf("mel_%s", parameters[3:4]))
  
  pairs <- as.data.frame(t(combn(params, 2)))
  
  for(i in c(1:nrow(pairs))){
    
    form <- as.formula(paste0("cbind(", pairs[i, 1], ", ", pairs[i, 2], ") ~ trait - 1"))
    
    for (condition in c("NoPhage", "Phage")){
      
      current_data <- parameters_all |>
        filter(Condition == condition) |>
        dplyr::select(all_of(as.character(unlist(pairs[i, ]))), animal, Bacteria, Block) |>
        na.omit()
      
      model_name <- sprintf("cor_%s_%s_%s.Rdata", condition, pairs[i, 1], pairs[i, 2])
      
      if (!file.exists(here("models", model_name))) {
        
        model <- MCMCglmm(fixed = form, random = ~ us(trait):animal + idh(trait):Bacteria:Block,
                          rcov = ~idh(trait):units, family = c("gaussian", "gaussian"), data = current_data, pedigree = tree,
                          prior = prior.corr, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
        
        save(model, file = here("models", model_name))
        
      } else {load(here("models", model_name))}
      
      var_p_x <- model$VCV[,1]
      
      var_p_y <- model$VCV[,4]
      
      cov_p <- model$VCV[,2]
      
      Bp <- as.mcmc(cov_p / var_p_x)
      
      Rp <- as.mcmc(cov_p /(sqrt(var_p_x * var_p_y)))
      
      Fx <- as.mcmc(model$Sol[,1])
      Fy <- as.mcmc(model$Sol[,2])
      
      correlations <- rbind(correlations, data.frame(Condition = condition,
                                                     Assumption = assumption,
                                                     x = pairs[i, 1],
                                                     y = pairs[i, 2],
                                                     
                                                     Fx_mean = mean(Fx),
                                                     Fx_low  = HPDinterval(Fx)[1],
                                                     Fx_high = HPDinterval(Fx)[2],
                                                     Fy_mean = mean(Fy),
                                                     Fy_low  = HPDinterval(Fy)[1],
                                                     Fy_high = HPDinterval(Fy)[2],
                                                     
                                                     Bp_mean = mean(Bp),
                                                     Bp_low = HPDinterval(Bp)[1],
                                                     Bp_high = HPDinterval(Bp)[2],
                                                     Rp_mean = mean(Rp),
                                                     Rp_low = HPDinterval(Rp)[1],
                                                     Rp_high = HPDinterval(Rp)[2]))
      
    }
    
  }
  
}

# ----- 4.3. Correlations in-vivo:in-vitro -------------------------------------

parameters_all_delta <- left_join(parameters_all_delta, data_invitro)

correlations_invitro <- data.frame(Condition = character(),
                                   Assumption = character(),
                                   x = character(),
                                   y = character(),
                                   Fx_mean = numeric(),
                                   Fx_low  = numeric(),
                                   Fx_high = numeric(),
                                   Fy_mean = numeric(),
                                   Fy_low  = numeric(),
                                   Fy_high = numeric(),
                                   Bp_mean = numeric(),
                                   Bp_low = numeric(),
                                   Bp_high = numeric(),
                                   Rp_mean = numeric(),
                                   Rp_low = numeric(),
                                   Rp_high = numeric())

assumptions <- c("unobserved", "cure")

parameters <- c("U", "t50", "k")

invitros <- c("OD", "qPCR")

for (assumption in assumptions){
  for (parameter in parameters){
    for (iv in invitros){
      
      file <- sprintf("cor_in_vitro_%s_%s_%s_%s.Rdata", "bivariate", assumption, parameter, iv)
      
      par <- sprintf("mor_%s_%s", assumption, parameter)
      
      form <- as.formula(paste0("cbind(", par, ", ", iv, ") ~ trait - 1"))
      
      if (!file.exists(here("models", file))) {
        
        model <- MCMCglmm(fixed = form, random = ~ us(trait):animal + us(trait):Bacteria:Block,
                          rcov = ~us(trait):units, family = c("gaussian", "gaussian"), data = parameters_all_delta, pedigree = tree,
                          prior = prior.corr, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
        
        save(model, file = here("models", file))
        
      } else {load(here("models", file))}
      
      var_p_x <- model$VCV[,1]
      
      var_p_y <- model$VCV[,4]
      
      cov_p <- model$VCV[,2]
      
      Bp <- as.mcmc(cov_p / var_p_x)
      
      Rp <- as.mcmc(cov_p /(sqrt(var_p_x * var_p_y)))
      
      Fx <- as.mcmc(model$Sol[,1])
      Fy <- as.mcmc(model$Sol[,2])
      
      correlations_invitro <- rbind(correlations_invitro, data.frame(Condition = condition,
                                                                     Assumption = assumption,
                                                                     x = parameter,
                                                                     y = iv,
                                                                     
                                                                     Fx_mean = mean(Fx),
                                                                     Fx_low  = HPDinterval(Fx)[1],
                                                                     Fx_high = HPDinterval(Fx)[2],
                                                                     Fy_mean = mean(Fy),
                                                                     Fy_low  = HPDinterval(Fy)[1],
                                                                     Fy_high = HPDinterval(Fy)[2],
                                                                     
                                                                     Bp_mean = mean(Bp),
                                                                     Bp_low = HPDinterval(Bp)[1],
                                                                     Bp_high = HPDinterval(Bp)[2],
                                                                     Rp_mean = mean(Rp),
                                                                     Rp_low = HPDinterval(Rp)[1],
                                                                     Rp_high = HPDinterval(Rp)[2]))
      
    }
    
  }
  
}

parameters <- c("Vmax", "K")

for (parameter in parameters){
  for (iv in invitros){
    
    file <- sprintf("cor_in_vitro_%s_%s_%s.Rdata", "bivariate", parameter, iv)
    
    par <- sprintf("mel_%s", parameter)
    
    form <- as.formula(paste0("cbind(", par, ", ", iv, ") ~ trait - 1"))
    
    if (!file.exists(here("models", file))) {
      
      model <- MCMCglmm(fixed = form, random = ~ us(trait):animal + us(trait):Bacteria:Block,
                        rcov = ~us(trait):units, family = c("gaussian", "gaussian"), data = parameters_all_delta, pedigree = tree,
                        prior = prior.corr, nitt = 130000*itt, thin = 50*itt, burnin = 30000*itt)
      
      save(model, file = here("models", file))
      
    } else {load(here("models", file))}
    
    var_p_x <- model$VCV[,1]
    
    var_p_y <- model$VCV[,4]
    
    cov_p <- model$VCV[,2]
    
    Bp <- as.mcmc(cov_p / var_p_x)
    
    Rp <- as.mcmc(cov_p /(sqrt(var_p_x * var_p_y)))
    
    Fx <- as.mcmc(model$Sol[,1])
    Fy <- as.mcmc(model$Sol[,2])
    
    correlations_invitro <- rbind(correlations_invitro, data.frame(Condition = condition,
                                                                   Assumption = assumption,
                                                                   x = parameter,
                                                                   y = iv,
                                                                   
                                                                   Fx_mean = mean(Fx),
                                                                   Fx_low  = HPDinterval(Fx)[1],
                                                                   Fx_high = HPDinterval(Fx)[2],
                                                                   Fy_mean = mean(Fy),
                                                                   Fy_low  = HPDinterval(Fy)[1],
                                                                   Fy_high = HPDinterval(Fy)[2],
                                                                   
                                                                   Bp_mean = mean(Bp),
                                                                   Bp_low = HPDinterval(Bp)[1],
                                                                   Bp_high = HPDinterval(Bp)[2],
                                                                   Rp_mean = mean(Rp),
                                                                   Rp_low = HPDinterval(Rp)[1],
                                                                   Rp_high = HPDinterval(Rp)[2]))
    
  }
  
}

