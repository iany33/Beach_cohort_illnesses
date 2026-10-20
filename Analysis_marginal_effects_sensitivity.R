
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

avg_comparisons(m_skin2, re_formula = NULL, variables = list(water_time_s = "iqr"))
avg_comparisons(m_skin2, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"))

avg_pred <- avg_predictions(m_skin2, type = "response", re_formula = NULL, variables = list(
  water_time_s = seq(-0.6142013, 9.2603016, by = 0.4)))

pred <- posterior_draws(avg_pred)

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
  reframe(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

quantile(data_follow$water_time_s, probs = c(0.25, 0.50, 0.75, 0.95), na.rm = TRUE)

avg_pred <- avg_predictions(m_skin2, type = "response", re_formula = NULL, variables = list(
  water_time_s = c(-0.6249203, -0.3283636, 0.1658975, 1.7475330),
  log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.4)))

pred <- posterior_draws(avg_pred)

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

avg_comp <- avg_comparisons(m_skin_3day, re_formula = NULL, variables = "water_exp_body")
avg_comp

mfx <- posterior_draws(avg_comp)
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
  xlim(-25, 30) -> RD_skin_3day_body


avg_comp <- avg_comparisons(m_skin_3day, re_formula = NULL, variables = "water_exp_body", 
                            comparison = "lnratioavg", transform = "exp")
avg_comp

mfx <- posterior_draws(avg_comp)

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
  xlim(0,3.5) -> RR_skin_3day_body


avg_comp <- avg_comparisons(m_skin_3day_any, re_formula = NULL, variables = "water_contact")
avg_comp

mfx <- posterior_draws(avg_comp)
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
  xlim(-25, 30) -> RD_skin_3day_any


avg_comp <- avg_comparisons(m_skin_3day_any, re_formula = NULL, 
                            comparison = "lnratio", transform = "exp", variables = "water_contact")
avg_comp

mfx <- posterior_draws(avg_comp)

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
  xlim(0,3.5) -> RR_skin_3day_any


avg_comp <- avg_comparisons(m_skin_3day_head, re_formula = NULL, variables = "water_exp_head")
avg_comp

mfx <- posterior_draws(avg_comp)
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
  xlim(-25, 30) -> RD_skin_3day_head


avg_comp <- avg_comparisons(m_skin_3day_head, re_formula = NULL, 
                            comparison = "lnratio", transform = "exp", variables = "water_exp_head")
avg_comp

mfx <- posterior_draws(avg_comp)

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
  xlim(0,3.5) -> RR_skin_3day_head



### Respiratory illness outcomes ###

# Time in the water 

quantile(data_follow$water_time, na.rm = TRUE)
quantile(data_follow$water_time_s, na.rm = TRUE)

list <- data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = mean(log_e_coli_max_s, na.rm=TRUE))
list <- as.list(list)

avg_comparisons(m_resp3.3, re_formula = NULL, variables = list(water_time_s = "iqr"))
avg_comparisons(m_resp3.3, re_formula = NULL, variables = list(log_e_coli_max_s = "iqr"))

avg_pred <- avg_predictions(m_resp3.3, type = "response", re_formula = NULL, variables = list(
  water_time_s = seq(-0.6142013, 9.2603016, by = 0.4)))

pred <- posterior_draws(avg_pred)

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
  reframe(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

quantile(data_follow$water_time_s, probs = c(0.25, 0.50, 0.75, 0.95), na.rm = TRUE)

avg_pred <- avg_predictions(m_resp3.3, type = "response", re_formula = NULL, variables = list(
  water_time_s = c(-0.6249203, -0.3283636, 0.1658975, 1.7475330),
  log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.4)))

pred <- posterior_draws(avg_pred)

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
  facet_wrap(~ water_time) -> time_ecoli_resp


# 3-day outcome variable

avg_comp <- avg_comparisons(m_resp_3day_body, re_formula = NULL, variables = "water_exp_body")
avg_comp

mfx <- posterior_draws(avg_comp)
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
  xlim(-25, 30) -> RD_resp_3day_body


avg_comp <- avg_comparisons(m_resp_3day_body, re_formula = NULL, variables = "water_exp_body", 
                            comparison = "lnratioavg", transform = "exp")
avg_comp

mfx <- posterior_draws(avg_comp)

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
  xlim(0,3.5) -> RR_resp_3day_body



avg_comp <- avg_comparisons(m_resp_3day_head, re_formula = NULL, variables = "water_exp_head")
avg_comp

mfx <- posterior_draws(avg_comp)
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
  xlim(-25, 30) -> RD_resp_3day_head


avg_comp <- avg_comparisons(m_resp_3day_head, re_formula = NULL, variables = "water_exp_head", 
                            comparison = "lnratioavg", transform = "exp")
avg_comp

mfx <- posterior_draws(avg_comp)

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
  xlim(0,3.5) -> RR_resp_3day_head



avg_comp <- avg_comparisons(m_resp_3day_any, re_formula = NULL, variables = "water_contact")
avg_comp

mfx <- posterior_draws(avg_comp)
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
  xlim(-25, 30) -> RD_resp_3day_any



avg_comp <- avg_comparisons(m_resp_3day_any, re_formula = NULL, 
                            comparison = "lnratio", transform = "exp", variables = "water_contact")
avg_comp

mfx <- posterior_draws(avg_comp)

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
  xlim(0,3.5) -> RR_resp_3day_any


### Combine RD plots ###

RD_skin_3day_any <- RD_skin_3day_any + ggtitle("Skin Symptoms")
RD_resp_3day_any <- RD_resp_3day_any + ggtitle("Respiratory Illness")

Fig_RD <- RD_skin_3day_any + RD_skin_3day_body + RD_skin_3day_head + 
  RD_resp_3day_any + RD_resp_3day_body + RD_resp_3day_head
Fig_RD + plot_layout(ncol = 1, axes = 'collect')

remove(RD_resp_3day_any, RD_resp_3day_head, RD_resp_3day_body, 
       RD_skin_3day_any, RD_skin_3day_body, RD_skin_3day_head)

### Combine RR plots ###

RR_skin_3day_any <- RR_skin_3day_any + ggtitle("Skin Symptoms")
RR_resp_3day_any <- RR_resp_3day_any + ggtitle("Respiratory Illness")

Fig_RR_3day <- RR_skin_3day_any + RR_skin_3day_body + RR_skin_3day_head + 
  RR_resp_3day_any + RR_resp_3day_body + RR_resp_3day_head
Fig_RR_3day + plot_layout(ncol = 1, axes = 'collect')

remove(RR_resp_3day_any, RR_resp_3day_head, RR_resp_3day_body, 
       RR_skin_3day_any, RR_skin_3day_body, RR_skin_3day_head)


## Combine time in water plots

time_skin <- time_skin + ggtitle("Skin Symptoms")
time_resp <- time_resp + ggtitle("Respiratory Illness")

Fig_time <- time_skin + time_resp
Fig_time + plot_layout(ncol = 1, axes = 'collect')

remove(time_skin, time_resp)


time_ecoli_skin <- time_ecoli_skin + ggtitle("Skin Symptoms")
time_ecoli_resp <- time_ecoli_resp + ggtitle("Respiratory Illness")

Fig_ecoli_time <- time_ecoli_skin + time_ecoli_resp
Fig_ecoli_time + plot_layout(ncol = 1, axes = 'collect')

remove(time_ecoli_skin, time_ecoli_resp)





### Respiratory illness outcomes ###

# One participant per household model comparison


avg_comp <- avg_comparisons(m_resp3.1.house, re_formula = NULL, variables = "water_exp_body")
avg_comp

mfx <- posterior_draws(avg_comp)
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
  xlim(-35, 45) -> RD_resp_house_body


avg_comp <- avg_comparisons(m_resp3.1.house, re_formula = NULL, variables = "water_exp_body", 
                            comparison = "lnratioavg", transform = "exp")
avg_comp

mfx <- posterior_draws(avg_comp)

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
  xlim(0,3.5) -> RR_resp_house_body


avg_comp <- avg_comparisons(m_resp3.2.house, re_formula = NULL, variables = "water_exp_head")
avg_comp

mfx <- posterior_draws(avg_comp)
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
  xlim(-35, 45) -> RD_resp_house_head


avg_comp <- avg_comparisons(m_resp3.2.house, re_formula = NULL, variables = "water_exp_head", 
                            comparison = "lnratioavg", transform = "exp")
avg_comp

mfx <- posterior_draws(avg_comp)

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
  xlim(0,3.5) -> RR_resp_house_head



avg_comp <- avg_comparisons(m_resp3.house, re_formula = NULL, variables = "water_contact")
avg_comp

mfx <- posterior_draws(avg_comp)
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
  xlim(-35, 45) -> RD_resp_house_any


avg_comp <- avg_comparisons(m_resp3.house, re_formula = NULL, 
                            comparison = "lnratio", transform = "exp", variables = "water_contact")
avg_comp

mfx <- posterior_draws(avg_comp)

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
  xlim(0,3.5) -> RR_resp_house_any


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


