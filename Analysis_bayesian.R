
### Bayesian analysis of beach cohort survey data 

pacman::p_load(
  Matrix,
  tidyverse,  
  rstatix,
  janitor,
  brms,
  tidybayes, 
  bayesplot,
  marginaleffects,
  cmdstanr,
  rstan,
  performance
)


## Respiratory illness outcome 
# Any water contact

# Set weakly informative priors and run model- starting with baseline model 

priors <- c(set_prior("normal(0,1.5)", class= "b"),
            set_prior("exponential(1)", class = "sd"))

m_resp <- brm(respiratory3 ~ 0 + Intercept + water_contact + 
                (1 | beach/recruit_date/house_id),
            family = bernoulli, data = data_follow, prior = priors,
            iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 5764, control = list(adapt_delta = 0.9),
            backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_resp, robust = TRUE)
get_variables(m_resp)
plot(m_resp)
pp_check(m_resp, ndraws=100)
pp_check(m_resp, type = "stat", stat = "mean")

conditional_effects(m_resp, effects = "water_contact") 

# Examine PAV-adjusted calibration plot diagnostic

y <- data_follow |> mutate(respiratory4 = if_else(respiratory3=="No", 0, 1)) 
y <- y$respiratory4
yrep <- posterior_predict(m_resp, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

# Model has poor predictive accuracy/fit
# Add full model-covariates and FIB interaction

m_resp2 <- brm(respiratory3 ~ 0 + Intercept + water_contact*log_e_coli_max_s + age5 + gender + education2 + ethnicity2 + cond_resp + 
               cond_immune + cond_allergy + other_rec_act + beach_exp_food + sand_contact + 
               (1 | beach/recruit_date/house_id),
             family = bernoulli, data = data_follow, prior = priors,
             iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 3273, control = list(adapt_delta = 0.95),
             backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_resp2, robust = TRUE)
get_variables(m_resp2)
plot(m_resp2)
pp_check(m_resp2, ndraws=100)
pp_check(m_resp2, type = "stat", stat = "mean")

conditional_effects(m_resp2, effects = "log_e_coli_max_s:water_contact")
conditional_effects(m_resp2, effects = "water_contact") -> fit
fit$water_contact

# Examine PAV-adjusted calibration plot diagnostic

y <- data_follow |> mutate(respiratory4 = if_else(respiratory3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_e_coli_max_s", "age5", "gender", "education2", "ethnicity2", "cond_resp", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$respiratory4
yrep <- posterior_predict(m_resp2, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

# Model has poor predictive accuracy - drop random effect for house and add house group variable instead

m_resp3 <- brm(respiratory3 ~ 0 + Intercept + water_contact*log_e_coli_max_s + age5 + gender + education2 + ethnicity2 + cond_resp + 
                 cond_immune + cond_allergy + other_rec_act + beach_exp_food + sand_contact + household_group +
                 (1 | site/beach/recruit_date),
               family = bernoulli, data = data_follow, prior = priors,
               iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 9696, control = list(adapt_delta = 0.95),
               backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_resp3, robust = TRUE)
get_variables(m_resp3)
plot(m_resp3)
pp_check(m_resp3, ndraws=100)
pp_check(m_resp3, type = "stat", stat = "mean")

conditional_effects(m_resp3, effects = "log_e_coli_max_s:water_contact")
conditional_effects(m_resp3, effects = "water_contact") -> fit
fit$water_contact

yrep <- posterior_predict(m_resp3, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])


# Body immersion exposure

m_resp3.1 <- brm(respiratory3 ~ 0 + Intercept + water_exp_body*log_e_coli_max_s + age5 + gender + education2 + ethnicity2 + cond_resp + 
                cond_immune + cond_allergy + other_rec_act + beach_exp_food + sand_contact + household_group +
                (1 | site/beach/recruit_date),
              family = bernoulli, data = data_follow, prior = priors,
              iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 1228, control = list(adapt_delta = 0.99),
              backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_resp3.1, robust = TRUE)
plot(m_resp3.1)
pp_check(m_resp3.1, ndraws=100)
pp_check(m_resp3.1, type = "stat", stat = "mean")

conditional_effects(m_resp3.1, effects = "log_e_coli_max_s:water_exp_body")
conditional_effects(m_resp3.1, effects = "water_exp_body") -> fit
fit$water_exp_body

yrep <- posterior_predict(m_resp3.1, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

# Head immersion exposure

m_resp3.2 <- brm(respiratory3 ~ 0 + Intercept + water_exp_head*log_e_coli_max_s + age5 + gender + education2 + ethnicity2 + cond_resp + 
                   cond_immune + cond_allergy + other_rec_act + beach_exp_food + sand_contact + household_group +
                   (1 | site/beach/recruit_date),
                 family = bernoulli, data = data_follow, prior = priors,
                 iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 7762, control = list(adapt_delta = 0.99),
                 backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_resp3.2, robust = TRUE)
get_variables(m_resp3.2)
plot(m_resp3.2)
pp_check(m_resp3.2, ndraws=100)
pp_check(m_resp3.2, type = "stat", stat = "mean")

conditional_effects(m_resp3.2, effects = "log_e_coli_max_s:water_exp_head")
conditional_effects(m_resp3.2, effects = "water_exp_head")

yrep <- posterior_predict(m_resp3.2, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

loo(m_resp3.1, m_resp3.2)


# Examine time in water as predictor

m_resp3.3 <- brm(respiratory3 ~ 0 + Intercept + water_time_s*log_e_coli_max_s + age5 + gender + education2 + ethnicity2 + cond_resp + 
                   cond_immune + cond_allergy + other_rec_act + beach_exp_food + sand_contact + household_group +
                   (1 | site/beach/recruit_date),
                 family = bernoulli, data = data_follow, prior = priors,
                 iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 3804, control = list(adapt_delta = 0.95),
                 backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_resp3.3, robust = TRUE)
get_variables(m_resp3.3)
plot(m_resp3.3)
pp_check(m_resp3.3, ndraws=100)
pp_check(m_resp3.3, type = "stat", stat = "mean")

conditional_effects(m_resp3.3, effects = "log_e_coli_max_s:water_time_s")
conditional_effects(m_resp3.3, effects = "water_time_s")

yrep <- posterior_predict(m_resp3.3, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

loo(m_resp3.1, m_resp3.2, m_resp3.3)



### Skin illness outcome ###
# Any water contact

m_skin <- brm(skin_infection3 ~ 0 + Intercept + water_contact*log_e_coli_max_s + age5 + gender + education2 + ethnicity2 +  
                cond_skin + cond_immune + cond_allergy + other_rec_act + beach_exp_sunscreen + beach_exp_repellent +
                sand_contact + household_group +
                (1 | site/beach/recruit_date),
              family = bernoulli, data = data_follow, prior = priors,
              iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 1014, control = list(adapt_delta = 0.95),
              backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_skin, robust = TRUE)
get_variables(m_skin)
plot(m_skin)
pp_check(m_skin, ndraws=100)
pp_check(m_skin, type = "stat", stat = "mean")

conditional_effects(m_skin, effects = "log_e_coli_max_s:water_contact")
conditional_effects(m_skin, effects = "water_contact") -> fit
fit$water_contact

y <- data_follow |> mutate(skin_infection4 = if_else(skin_infection3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_e_coli_max_s", "age5", "gender", "education2", "ethnicity2", "cond_skin", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$skin_infection4
yrep <- posterior_predict(m_skin, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

# Body immersion exposure

m_skin1 <- brm(skin_infection3 ~ 0 + Intercept + water_exp_body*log_e_coli_max_s + age5 + gender + education2 + ethnicity2 +  
                cond_skin + cond_immune + cond_allergy + other_rec_act + beach_exp_sunscreen + beach_exp_repellent +
                sand_contact + household_group +
                (1 | site/beach/recruit_date),
              family = bernoulli, data = data_follow, prior = priors,
              iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 9180, control = list(adapt_delta = 0.95),
              backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_skin1, robust = TRUE)
plot(m_skin1)
pp_check(m_skin1, ndraws=100)
pp_check(m_skin1, type = "stat", stat = "mean")

conditional_effects(m_skin1, effects = "log_e_coli_max_s:water_exp_body")
conditional_effects(m_skin1, effects = "water_exp_body") -> fit
fit$water_exp_body

yrep <- posterior_predict(m_skin1, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])


m_skin1.2 <- brm(skin_infection3 ~ 0 + Intercept + water_exp_head*log_e_coli_max_s + age5 + gender + education2 + ethnicity2 +  
                 cond_skin + cond_immune + cond_allergy + other_rec_act + beach_exp_sunscreen + beach_exp_repellent +
                 sand_contact + household_group +
                 (1 | site/beach/recruit_date),
               family = bernoulli, data = data_follow, prior = priors,
               iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 5502, control = list(adapt_delta = 0.99),
               backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_skin1.2, robust = TRUE)
plot(m_skin1.2)
pp_check(m_skin1.2, ndraws=100)
pp_check(m_skin1.2, type = "stat", stat = "mean")

conditional_effects(m_skin1.2, effects = "log_e_coli_max_s:water_exp_head")
conditional_effects(m_skin1.2, effects = "water_exp_head") -> fit
fit$water_exp_body

loo(m_skin1, m_skin1.2)


# Examine time in water as predictor

m_skin2 <- brm(skin_infection3 ~ 0 + Intercept + water_time_s*log_e_coli_max_s + age5 + gender + education2 + ethnicity2 +  
                 cond_skin + cond_immune + cond_allergy + other_rec_act + beach_exp_sunscreen + beach_exp_repellent +
                 sand_contact + household_group +
                 (1 | site/beach/recruit_date),
                 family = bernoulli, data = data_follow, prior = priors,
                 iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 5173, control = list(adapt_delta = 0.95),
                 backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_skin2, robust = TRUE)
get_variables(m_skin2)
plot(m_skin2)
pp_check(m_skin2, ndraws=100)
pp_check(m_skin2, type = "stat", stat = "mean")

conditional_effects(m_skin2, effects = "log_e_coli_max_s:water_time_s")
conditional_effects(m_skin2, effects = "water_time_s")

yrep <- posterior_predict(m_skin2, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

loo(m_skin1, m_skin2)

# Re-run body immersion model with other FIB

m_skin_entero <- brm(skin_infection3 ~ 0 + Intercept + water_exp_body*log_entero_max_s + age5 + gender + education2 + ethnicity2 +  
                 cond_skin + cond_immune + cond_allergy + other_rec_act + beach_exp_sunscreen + beach_exp_repellent +
                 sand_contact + household_group +
                 (1 | site/beach/recruit_date),
                 family = bernoulli, data = data_follow, prior = priors,
                 iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 6782, control = list(adapt_delta = 0.99),
                 backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_skin_entero , robust = TRUE)
get_variables(m_skin_entero)
plot(m_skin_entero)
pp_check(m_skin_entero, ndraws=100)
pp_check(m_skin_entero, type = "stat", stat = "mean")

conditional_effects(m_skin_entero, effects = "log_entero_max_s:water_exp_body")
conditional_effects(m_skin_entero, effects = "water_exp_body")

y <- data_follow |> mutate(skin_infection3 = if_else(skin_infection3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_entero_max_s", "age5", "gender", "education2", "ethnicity2", "cond_resp", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$skin_infection3
yrep <- posterior_predict(m_skin_entero, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

# Examine with turbidity as water quality indicator

m_skin_turb <- brm(skin_infection3 ~ 0 + Intercept + water_exp_body*log_turbidity_s + age5 + gender + education2 + ethnicity2 +  
                 cond_skin + cond_immune + cond_allergy + other_rec_act + beach_exp_sunscreen + beach_exp_repellent +
                 sand_contact + household_group +
                 (1 | site/beach/recruit_date),
               family = bernoulli, data = data_follow, prior = priors,
               iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 4591, control = list(adapt_delta = 0.99),
               backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_skin_turb, robust = TRUE)
get_variables(m_skin_turb)
plot(m_skin_turb)
pp_check(m_skin_turb, ndraws=100)
pp_check(m_skin_turb, type = "stat", stat = "mean")

conditional_effects(m_skin_turb, effects = "log_turbidity_s:water_exp_body")
conditional_effects(m_skin_turb, effects = "water_exp_body")

y <- data_follow |> mutate(skin_infection3 = if_else(skin_infection3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_turbidity_s", "age5", "gender", "education2", "ethnicity2", "cond_resp", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$skin_infection3
yrep <- posterior_predict(m_skin_turb, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])

# Examine MST markers

m_skin_human <- brm(skin_infection3 ~ 0 + Intercept + water_exp_body*log_mst_human_max_s + age5 + gender + education2 + ethnicity2 +  
                 cond_skin + cond_immune + cond_allergy + other_rec_act + beach_exp_sunscreen + beach_exp_repellent +
                 sand_contact + household_group +
                 (1 | site/beach/recruit_date),
               family = bernoulli, data = data_follow, prior = priors,
               iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 3042, control = list(adapt_delta = 0.95),
               backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_skin_human, robust = TRUE)
get_variables(m_skin_human)
plot(m_skin_human)
pp_check(m_skin_human, ndraws=100)
pp_check(m_skin_human, type = "stat", stat = "mean")

conditional_effects(m_skin_human, effects = "log_mst_human_max_s:water_exp_body")
conditional_effects(m_skin_human, effects = "water_exp_body")

y <- data_follow |> mutate(skin_infection3 = if_else(skin_infection3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_mst_human_max_s", "age5", "gender", "education2", "ethnicity2", "cond_resp", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$skin_infection3
yrep <- posterior_predict(m_skin_human, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])


m_skin_human_mt <- brm(skin_infection3 ~ 0 + Intercept + water_exp_body*log_mst_human_mt_max_s + age5 + gender + education2 + ethnicity2 +  
                 cond_skin + cond_immune + cond_allergy + other_rec_act + beach_exp_sunscreen + beach_exp_repellent +
                 sand_contact + household_group +
                 (1 | site/beach/recruit_date),
               family = bernoulli, data = data_follow, prior = priors,
               iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 5817, control = list(adapt_delta = 0.95),
               backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_skin_human_mt, robust = TRUE)
get_variables(m_skin_human_mt)
plot(m_skin_human_mt)
pp_check(m_skin_human_mt, ndraws=100)
pp_check(m_skin_human_mt, type = "stat", stat = "mean")

conditional_effects(m_skin_human_mt, effects = "log_mst_human_mt_max_s:water_exp_body")
conditional_effects(m_skin_human_mt, effects = "water_exp_body")

y <- data_follow |> mutate(skin_infection3 = if_else(skin_infection3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_mst_human_mt_max_s", "age5", "gender", "education2", "ethnicity2", "cond_resp", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$skin_infection3
yrep <- posterior_predict(m_skin_human_mt, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])


m_skin_gull <- brm(skin_infection3 ~ 0 + Intercept + water_exp_body*log_mst_gull_max_s + age5 + gender + education2 + ethnicity2 +  
                 cond_skin + cond_immune + cond_allergy + other_rec_act + beach_exp_sunscreen + beach_exp_repellent +
                 sand_contact + household_group +
                 (1 | site/beach/recruit_date),
               family = bernoulli, data = data_follow, prior = priors,
               iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 2400, control = list(adapt_delta = 0.95),
               backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_skin_gull, robust = TRUE)
get_variables(m_skin_gull)
plot(m_skin_gull)
pp_check(m_skin_gull, ndraws=100)
pp_check(m_skin_gull, type = "stat", stat = "mean")

conditional_effects(m_skin_gull, effects = "log_mst_gull_max_s:water_exp_body")
conditional_effects(m_skin_gull, effects = "water_exp_body")

y <- data_follow |> mutate(skin_infection3 = if_else(skin_infection3=="No", 0, 1)) |> 
  drop_na(any_of(c("log_mst_gull_max_s", "age5", "gender", "education2", "ethnicity2", "cond_resp", 
                   "other_rec_act", "beach_exp_food", "sand_contact")))
y <- y$skin_infection3
yrep <- posterior_predict(m_skin_gull, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])



# Sensitivity analysis

m_skin_3day <- brm(skin_infection_3day ~ 0 + Intercept + water_exp_body*log_e_coli_max_s + age5 + gender + education2 + ethnicity2 +  
                 cond_skin + cond_immune + cond_allergy + other_rec_act + beach_exp_sunscreen + beach_exp_repellent +
                 sand_contact + household_group +
                 (1 | site/beach/recruit_date),
               family = bernoulli, data = data_follow, prior = priors,
               iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 7749, control = list(adapt_delta = 0.95),
               backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m_skin_3day, robust = TRUE)
plot(m_skin_3day)
pp_check(m_skin_3day, ndraws=100)
pp_check(m_skin_3day, type = "stat", stat = "mean")

conditional_effects(m_skin_3day, effects = "log_e_coli_max_s:water_exp_body")
conditional_effects(m_skin_3day, effects = "water_exp_body") -> fit
fit$water_exp_body

yrep <- posterior_predict(m_skin_3day, draws = 500)
ppc_calibration_pava(y, 0.95, yrep[1:100,])


### Negative control - skin infection

data_neg_control <- data_follow |> filter(water_contact3 == "No contact")

m.nc <- brm(skin_infection3 ~ 0 + Intercept + log_e_coli_max_s + age5 + gender + education2 + ethnicity2 +  
              cond_skin + cond_immune + cond_allergy + other_rec_act + beach_exp_sunscreen + beach_exp_repellent +
              sand_contact + household_group +
              (1 | site/beach/recruit_date),
              family = bernoulli, data = data_neg_control, prior = priors,
              iter = 2000, chains = 4, cores = 4, warmup = 1000, seed = 1868, control = list(adapt_delta = 0.99),
              backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1")))

summary(m.nc, robust = TRUE)
plot(m.nc)
plot(m.nc, variable = "b_log_e_coli_s")
plot(m.nc, variable = "b_log_e_coli_s", combo = c("dens", "rank_overlay"))
pp_check(m.nc, ndraws=100)
pp_check(m.nc, type = "stat", stat = "mean")

conditional_effects(m.nc, effects = "log_e_coli_max_s")



