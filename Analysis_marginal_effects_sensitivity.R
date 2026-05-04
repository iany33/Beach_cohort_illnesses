
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
            water_time_s = seq(-0.6249203, 9.2603016, by = 0.4),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_skin = "No", cond_immune = "No",
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
  theme(legend.position = "bottom") -> time_skin


data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_time_s = quantile(water_time_s, probs = c(0.25, 0.50, 0.75, 0.95), na.rm = TRUE),
            log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2), 
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

avg_slopes(m_skin2, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_time_s")

pred <- predictions(m_skin2, re_formula = NA, type = "response", newdata = nd) |> get_draws()

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) |> 
  mutate(water_time = water_time_s*sd(data_follow$water_time, na.rm=TRUE) + mean(data_follow$water_time, na.rm=TRUE)) |> 
  mutate(water_time = round(water_time, digits = 1)) |> 
  mutate(water_time = ifelse(water_time == -0.5, 0, water_time))

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of Illness",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "none") +
  facet_wrap(~ water_time) -> time_ecoli_skin


# 3-day outcome variable comparison

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_exp_body = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

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
  xlim(-20, 30) -> RD_skin_3day_body

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
  xlim(0,3) -> RR_skin_3day_body


nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_contact = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

mfx <- comparisons(m_skin_3day_any, re_formula = NA, variables = "water_contact", by = "water_contact", 
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
  xlim(-20, 30) -> RD_skin_3day_any

mfx <- comparisons(m_skin_3day_any, re_formula = NA, comparison = "lnratio", transform = "exp", 
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
  xlim(0,3) -> RR_skin_3day_any


nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_exp_head = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

mfx <- comparisons(m_skin_3day_head, re_formula = NA, variables = "water_exp_head", by = "water_exp_head", 
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
  xlim(-20, 30) -> RD_skin_3day_head


mfx <- comparisons(m_skin_3day_head, re_formula = NA, comparison = "lnratio", transform = "exp", 
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
  xlim(0,3) -> RR_skin_3day_head



### Respiratory illness outcomes ###

# Time in the water 

quantile(data_follow$water_time, na.rm = TRUE)
quantile(data_follow$water_time_s, na.rm = TRUE)

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_time_s = seq(-0.6249203, 9.2603016, by = 0.4),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_resp = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            sand_contact = "No", household_group = "Yes") 

avg_comparisons(m_resp3.3, re_formula = NA, variables = list(water_time_s = "iqr"), newdata = nd)
avg_comparisons(m_resp3.3, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd)

pred <- predictions(m_resp3.3, re_formula = NA, type = "response", newdata = nd) |> get_draws()

pred <- pred |> 
  mutate(water_time = water_time_s*sd(data_follow$water_time, na.rm=TRUE) + mean(data_follow$water_time, na.rm=TRUE))

ggplot(pred, aes(x = water_time, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Amount of Time in the Water (Min)",
       y = "Predicted Probability of Illness",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") -> time_resp


data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(water_time_s = quantile(water_time_s, probs = c(0.25, 0.50, 0.75, 0.95), na.rm = TRUE),
            log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2), 
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_resp = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            sand_contact = "No", household_group = "Yes") 

avg_slopes(m_resp3.3, re_formula = NA, variables = "log_e_coli_max_s", newdata = nd, by = "water_time_s")

pred <- predictions(m_resp3.3, re_formula = NA, type = "response", newdata = nd) |> get_draws()

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) |> 
  mutate(water_time = water_time_s*sd(data_follow$water_time, na.rm=TRUE) + mean(data_follow$water_time, na.rm=TRUE)) |> 
  mutate(water_time = round(water_time, digits = 1)) |> 
  mutate(water_time = ifelse(water_time == -0.5, 0, water_time))

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of Illness",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_time) -> time_ecoli_resp



# 3-day outcome variable

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
            education2 = "bachelors", cond_resp = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            sand_contact = "No", household_group = "Yes") 

mfx <- comparisons(m_resp_3day_body, re_formula = NA, variables = "water_exp_body", by = "water_exp_body", 
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
  xlim(-20, 30) -> RD_resp_3day_body

mfx <- comparisons(m_resp_3day_body, re_formula = NA, comparison = "lnratio", transform = "exp", 
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
  xlim(0,3) -> RR_resp_3day_body


nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_exp_head = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_resp = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            sand_contact = "No", household_group = "Yes") 

mfx <- comparisons(m_resp_3day_head, re_formula = NA, variables = "water_exp_head", by = "water_exp_head", 
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
  xlim(-20, 30) -> RD_resp_3day_head

mfx <- comparisons(m_resp_3day_head, re_formula = NA, comparison = "lnratio", transform = "exp", 
                   variables = "water_exp_head",   
                   by = "water_exp_head", newdata = nd) |> posterior_draws()

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Yes) / mean(No))" = "Head Immersion"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Risk Ratio", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#21918c") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(0,3) -> RR_resp_3day_head


nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_contact = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_resp = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            sand_contact = "No", household_group = "Yes") 

mfx <- comparisons(m_resp_3day_any, re_formula = NA, variables = "water_contact", by = "water_contact", 
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
  xlim(-20, 30) -> RD_resp_3day_any

mfx <- comparisons(m_resp_3day_any, re_formula = NA, comparison = "lnratio", transform = "exp", 
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
  xlim(0,3) -> RR_resp_3day_any


### Combine RD plots ###

RD_skin_3day_any <- RD_skin_3day_any + ggtitle("Skin Infection")
RD_resp_3day_any <- RD_resp_3day_any + ggtitle("Respiratory Illness")

Fig_RD <- RD_skin_3day_any + RD_skin_3day_body + RD_skin_3day_head + 
  RD_resp_3day_any + RD_resp_3day_body + RD_resp_3day_head
Fig_RD + plot_layout(ncol = 1, axes = 'collect')

remove(RD_resp_3day_any, RD_resp_3day_head, RD_resp_3day_body, 
       RD_skin_3day_any, RD_skin_3day_body, RD_skin_3day_head)

### Combine RR plots ###

RR_skin_3day_any <- RR_skin_3day_any + ggtitle("Skin Infection")
RR_resp_3day_any <- RR_resp_3day_any + ggtitle("Respiratory Illness")

Fig_RR_3day <- RR_skin_3day_any + RR_skin_3day_body + RR_skin_3day_head + 
  RR_resp_3day_any + RR_resp_3day_body + RR_resp_3day_head
Fig_RR_3day + plot_layout(ncol = 1, axes = 'collect')

remove(RR_resp_3day_any, RR_resp_3day_head, RR_resp_3day_body, 
       RR_skin_3day_any, RR_skin_3day_body, RR_skin_3day_head)


## Combine time in water plots

time_skin <- time_skin + ggtitle("Skin Infection")
time_resp <- time_resp + ggtitle("Respiratory Illness")

Fig_time <- time_skin + time_resp
Fig_time + plot_layout(ncol = 1, axes = 'collect')

remove(time_skin, time_resp)


time_ecoli_skin <- time_ecoli_skin + ggtitle("Skin Infection")
time_ecoli_resp <- time_ecoli_resp + ggtitle("Respiratory Illness")

Fig_ecoli_time <- time_ecoli_skin + time_ecoli_resp
Fig_ecoli_time + plot_layout(ncol = 1, axes = 'collect')

remove(time_ecoli_skin, time_ecoli_resp)





### Respiratory illness outcomes ###

# One participant per household model comparison

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
            education2 = "bachelors", cond_resp = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            sand_contact = "No", household_group = "Yes") 

mfx <- comparisons(m_resp3.1.house, re_formula = NA, variables = "water_exp_body", by = "water_exp_body", 
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
  xlim(-20, 40) -> RD_resp_house_body

mfx <- comparisons(m_resp3.1.house, re_formula = NA, comparison = "lnratio", transform = "exp", 
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
  xlim(0,3) -> RR_resp_house_body


nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_exp_head = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_resp = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            sand_contact = "No", household_group = "Yes") 

mfx <- comparisons(m_resp3.2.house, re_formula = NA, variables = "water_exp_head", by = "water_exp_head", 
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
  xlim(-20, 40) -> RD_resp_house_head

mfx <- comparisons(m_resp3.2.house, re_formula = NA, comparison = "lnratio", transform = "exp", 
                   variables = "water_exp_head",   
                   by = "water_exp_head", newdata = nd) |> posterior_draws()

mfx <- mfx |> 
  mutate(contrast = recode(contrast, "ln(mean(Yes) / mean(No))" = "Head Immersion"))

ggplot(mfx, aes(x = draw, y = contrast, fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Risk Ratio", y="") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#21918c") +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  xlim(0,3) -> RR_resp_house_head


nd <- data_follow |> 
  data_grid(log_e_coli_max_s = list$log_e_coli_max_s, 
            water_contact = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            ethnicity2 = "White", education2 = "bachelors", cond_resp = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            sand_contact = "No", household_group = "Yes") 

mfx <- comparisons(m_resp3.house, re_formula = NA, variables = "water_contact", by = "water_contact", 
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
  xlim(-20, 40) -> RD_resp_house_any

mfx <- comparisons(m_resp3.house, re_formula = NA, comparison = "lnratio", transform = "exp", 
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
  xlim(0,3) -> RR_resp_house_any


### Combine RD plots ###

RD_resp_house_any <- RD_resp_house_any + ggtitle("Respiratory Illness")

Fig_RD_house <- RD_resp_house_any + RD_resp_house_body + RD_resp_house_head
Fig_RD_house + plot_layout(ncol = 1, axes = 'collect')

remove(RD_resp_house_any, RD_resp_house_head, RD_resp_house_body)

### Combine RR plots ###

RR_resp_house_any <- RR_resp_house_any + ggtitle("Respiratory Illness")

Fig_RR_house <- RR_resp_house_any + RR_resp_house_body + RR_resp_house_head
Fig_RR_house + plot_layout(ncol = 1, axes = 'collect')

remove(RR_resp_house_any, RR_resp_house_head, RR_resp_house_body)

remove(Fig_RD_house, Fig_RR_house)


