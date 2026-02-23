
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

### Conditional and marginal effects with marginaleffects R package
# Predictions ignore cluster-level variables (re_formula = NA) to get overall averages

### Skin Infections and Body Immersion

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
            ethnicity2 = "White", education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

predictions(m_skin1, re_formula = NA, by = "water_exp_body", type = "response", newdata = nd)

pred <- predictions(m_skin1, re_formula = NA, by = "water_exp_body", type = "response", newdata = nd) |> get_draws()

pred <- pred |> mutate(draw = draw*1000)

ggplot(pred, aes(x = draw, y = water_exp_body, fill = water_exp_body)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "Predicted Skin Infection Incident Risk per 1000 Beachgoers", y = "Body Immersion Status",
       subtitle = "Posterior Predictions", fill = "Water contact") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  xlim(0, 75)  

# Examine marginal effects/contrast of water contact exposure effect - probability scale

avg_comparisons(m_skin1, re_formula = NA, variables = "water_exp_body", newdata = nd)

mfx <- comparisons(m_skin1, re_formula = NA, variables = "water_exp_body", by = "water_exp_body", 
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
  xlim(-25, 40) -> Skin_body_RD

# Check proportion of posterior that is greater than 0 and other values

mfx |> group_by(contrast) |> 
  summarize(proportion_0 = mean(draw > 0),
            proportion_1 = mean(draw > 1),
            proportion_5 = mean(draw > 5),
            proportion_10 = mean(draw > 10),
            proportion_20 = mean(draw > 20))

# Population-averaged (marginal) adjusted risk ratios

avg_comparisons(m_skin1, re_formula = NA, variables = "water_exp_body", newdata = nd,
                comparison = "lnratioavg", transform = "exp")

mfx <- comparisons(m_skin1, re_formula = NA, comparison = "lnratio", transform = "exp", 
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
  xlim(0,4) -> Skin_body


# Gender specific estimates 

avg_comparisons(m_skin1, re_formula = NA, variables = "water_exp_body", newdata = nd, by = "gender")

mfx <- comparisons(m_skin1, re_formula = NA, variables = "water_exp_body", by = "gender",
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(gender = recode(gender, "man/boy" = "Man/boy", "woman/girl" = "Woman/girl",
                         "fluid/trans" = "Fluid/trans"))

ggplot(mfx, aes(x = draw, y = gender, fill = gender)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect of Water Contact on Skin Infection Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-5, 50) +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast)

# Age specific estimates 

avg_comparisons(m_skin1, re_formula = NA, variables = "water_exp_body", newdata = nd, by = "age5")

mfx <- comparisons(m_skin1, re_formula = NA, variables = "water_exp_body", by = "age5",
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(age5 = fct_relevel(age5, "0-9", "10-19", "20+"))

ggplot(mfx, aes(x = draw, y = age5, fill = age5)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect of Water Contact on Skin Infection Incident Risk per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE, option = "turbo") +
  facet_wrap(~ contrast) +
  xlim(-5, 50) 

# Compare to other exposure measures

nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_contact = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity2 = "White", education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

mfx <- comparisons(m_skin, re_formula = NA, variables = "water_contact", by = "water_contact", 
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Yes - No" = "Any Water Contact"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Risk Difference per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#440154") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(-25, 40) -> Skin_any_RD

mfx <- comparisons(m_skin, re_formula = NA, comparison = "lnratio", transform = "exp", 
                   variables = "water_contact",   
                   by = "water_contact", newdata = nd) |> posterior_draws()

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Yes) / mean(No))" = "Any Water Contact"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Risk Ratio", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#440154") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(0,4) -> Skin_any


nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_exp_head = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity2 = "White", education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

mfx <- comparisons(m_skin1.2, re_formula = NA, variables = "water_exp_head", by = "water_exp_head", 
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Yes - No" = "Head Immersion"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Risk Difference per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#440154") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(-25, 40) -> Skin_head_RD


mfx <- comparisons(m_skin1.2, re_formula = NA, comparison = "lnratio", transform = "exp", 
                   variables = "water_exp_head",   
                   by = "water_exp_head", newdata = nd) |> posterior_draws()

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Yes) / mean(No))" = "Head immersion"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Risk Ratio", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#440154") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(0,4) -> Skin_head


Fig_skin_RR <- Skin_any  + Skin_body + Skin_head
Fig_skin_RR + plot_layout(ncol = 1, axes = 'collect')

remove(Skin_any, Skin_body, Skin_head)



### Respiratory illness outcomes ###


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

# Examine marginal effects/contrast of water contact exposure effect - probability scale

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
  xlim(-25, 40) -> Resp_body_RD

# Check proportion of posterior that is greater than 0 and other values

mfx |> group_by(contrast) |> 
  summarize(proportion_0 = mean(draw > 0),
            proportion_1 = mean(draw > 1),
            proportion_5 = mean(draw > 5),
            proportion_10 = mean(draw > 10),
            proportion_20 = mean(draw > 20))

# Population-averaged (marginal) adjusted risk ratios

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
  xlim(0,4) -> Resp_body

# Compare to other exposure measures

nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_contact = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity2 = "White", education2 = "bachelors", cond_resp = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            sand_contact = "No", household_group = "Yes") 

mfx <- comparisons(m_resp3, re_formula = NA, variables = "water_contact", by = "water_contact", 
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Yes - No" = "Any Water Contact"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Risk Difference per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#21918c") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(-25, 40) -> Resp_any_RD


mfx <- comparisons(m_resp3, re_formula = NA, comparison = "lnratio", transform = "exp", 
                   variables = "water_contact",   
                   by = "water_contact", newdata = nd) |> posterior_draws()

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Yes) / mean(No))" = "Any Water Contact"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Risk Ratio", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#21918c") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(0,4) -> Resp_any


nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_exp_head = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity2 = "White", education2 = "bachelors", cond_resp = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            sand_contact = "No", household_group = "Yes") 


mfx <- comparisons(m_resp3.2, re_formula = NA, variables = "water_exp_head", by = "water_exp_head", 
                   newdata = nd) |> posterior_draws()

mfx <- mfx |> mutate(draw = draw*1000)

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "Yes - No" = "Head Immersion"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Risk Difference per 1000 Beachgoers", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#21918c") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(-25, 40) -> Resp_head_RD


mfx <- comparisons(m_resp3.2, re_formula = NA, comparison = "lnratio", transform = "exp", 
                   variables = "water_exp_head",   
                   by = "water_exp_head", newdata = nd) |> posterior_draws()

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Yes) / mean(No))" = "Head immersion"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Risk Ratio", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#21918c") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(0,4) -> Resp_head

Fig_resp_RR <- Resp_any  + Resp_body + Resp_head
Fig_resp_RR + plot_layout(ncol = 1, axes = 'collect')


### Combine RR plots ###

Skin_any <- Skin_any + ggtitle("Skin Infection")
Resp_any <- Resp_any + ggtitle("Respiratory Illness")

Fig_RR <- Skin_any + Skin_body + Skin_head + Resp_any + Resp_body + Resp_head
Fig_RR + plot_layout(ncol = 1, axes = 'collect')

remove(Resp_any, Resp_body, Resp_head)
remove(Skin_any, Skin_body, Skin_head)
remove(Fig_resp_RR, Fig_skin_RR)


### Combine RD plots ###

Skin_any_RD <- Skin_any_RD + ggtitle("Skin Infection")
Resp_any_RD <- Resp_any_RD + ggtitle("Respiratory Illness")

Fig_RD <- Skin_any_RD + Skin_body_RD + Skin_head_RD + Resp_any_RD + Resp_body_RD + Resp_head_RD
Fig_RD + plot_layout(ncol = 1, axes = 'collect')

remove(Resp_any_RD, Resp_body_RD, Resp_head_RD)
remove(Skin_any_RD, Skin_body_RD, Skin_head_RD)
remove(Fig_resp_RD, Fig_skin_RD)





