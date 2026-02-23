
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
  modelr,
  patchwork,
  rstan,
  viridis
)

### Conditional and marginal effects with 'marginaleffects' R package
# Predictions ignore cluster-level variables (re_formula = NA) to get overall averages

### Skin Infections 

# Time in the water exposure

quantile(data_follow$water_time, na.rm = TRUE)
quantile(data_follow$water_time_s, na.rm = TRUE)

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_time_s = seq(-0.6142013, 9.2603016, by = 0.4),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity2 = "White", education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

avg_comparisons(m_skin2, re_formula = NA, variables = list(water_time_s = "iqr"), newdata = nd)
avg_comparisons(m_skin2, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd)

pred <- predictions(m_skin2, re_formula = NA, type = "response", newdata = nd) |> get_draws()

pred <- pred |> 
  mutate(water_time = water_time_s*sd(data_follow$water_time, na.rm=TRUE) + mean(data_follow$water_time, na.rm=TRUE))

ggplot(pred, aes(x = water_time, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Amount of Time in the Water (Min)",
       y = "Predicted Probability of Illness",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") 


data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_time_s = quantile(water_time_s, probs = c(0.25, 0.50, 0.75, 0.95), na.rm = TRUE),
            log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2), 
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity2 = "White", education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

avg_slopes(m_skin2, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_time_s")

pred <- predictions(m_skin2, re_formula = NA, type = "response", newdata = nd) |> get_draws()

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) |> 
  mutate(water_time = water_time_s*sd(data_follow$water_time, na.rm=TRUE) + mean(data_follow$water_time, na.rm=TRUE))

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of Illness",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_time)


# 3-day outcome variable

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_exp_body = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity2 = "White", education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

avg_comparisons(m_skin_3day, re_formula = NA, variables = "water_exp_body", newdata = nd)

mfx <- comparisons(m_skin_3day, re_formula = NA, variables = "water_exp_body", by = "water_exp_body", 
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Yes - No" = "Body Immersion"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Risk Difference per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#440154") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(-25, 40) 

avg_comparisons(m_skin_3day, re_formula = NA, variables = "water_exp_body", newdata = nd,
                comparison = "lnratioavg", transform = "exp")

mfx <- comparisons(m_skin_3day, re_formula = NA, comparison = "lnratio", transform = "exp", 
                   variables = "water_exp_body",   
                   by = "water_exp_body", newdata = nd) |> posterior_draws()

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Yes) / mean(No))" = "Body Immersion"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Risk Ratio", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#440154") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(0,4)












### Respiratory illness outcomes ###

# Time in the water 

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

exp(list$log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_exp_body = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity2 = "White", education2 = "bachelors", cond_resp = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            sand_contact = "No", household_group = "Yes") 

predictions(m_resp3.1, re_formula = NA, by = "water_exp_body", type = "response", newdata = nd)

pred <- predictions(m_resp3.1, re_formula = NA, by = "water_exp_body", type = "response", newdata = nd) |> get_draws()

pred <- pred |> mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = water_exp_body, fill = water_exp_body)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Predicted Respiratory Illness Incident Risk per 1000 Beachgoers", y = "Body Immersion Status",
       subtitle = "Posterior Predictions", fill = "Water contact") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0, 75)  



# 3-day outcome variable

avg_comparisons(m_resp3.1, re_formula = NA, variables = "water_exp_body", newdata = nd)

mfx <- comparisons(m_resp3.1, re_formula = NA, variables = "water_exp_body", by = "water_exp_body", 
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Yes - No" = "Body Immersion"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Risk Difference per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#21918c") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(-25, 40) 


avg_comparisons(m_resp3.1, re_formula = NA, variables = "water_exp_body", newdata = nd,
                comparison = "lnratioavg", transform = "exp")

mfx <- comparisons(m_resp3.1, re_formula = NA, comparison = "lnratio", transform = "exp", 
                   variables = "water_exp_body",   
                   by = "water_exp_body", newdata = nd) |> posterior_draws()

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Yes) / mean(No))" = "Body Immersion"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Risk Ratio", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#21918c") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(0,4) 



