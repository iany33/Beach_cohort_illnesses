
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

### Skin Infections

# Predicted probabilities of E. coli
# Sequence E. coli by range of logged, standardized and centered variable then back-transform

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2), 
            water_exp_body = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

pred <- predictions(m_skin1, re_formula = NA, by = c("water_exp_body", "log_e_coli_max_s"), 
                    type = "response", newdata = nd) |> get_draws()

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of Skin Infection",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_exp_body)   -> Skin_ecoli

avg_comparisons(m_skin1, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd)
avg_comparisons(m_skin1, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd, by = "water_exp_body")


# qPCR enterococci model 

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_entero_max_s = range(log_entero_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(log_entero_max_s = seq(-2.014249, 3.1953, by = 0.4), 
            water_exp_body = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

pred <- predictions(m_skin_entero, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(entero = exp(log_entero_max_s*sd(data_follow$log_entero_max, na.rm=TRUE) + mean(data_follow$log_entero_max, na.rm=TRUE))) 

pred <- pred |> 
  mutate(log_entero_max = log_entero_max_s*sd(data_follow$log_entero_max, na.rm=TRUE) + mean(data_follow$log_entero_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_entero_max, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Log Enterococci Highest Single Sample",
       y = "Predicted Probability of Skin Infection",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_exp_body)  -> Skin_entero

avg_comparisons(m_skin_entero, re_formula = NA, variables = list(log_entero_max_s = "iqr"), newdata = nd)
avg_comparisons(m_skin_entero, re_formula = NA, variables = list(log_entero_max_s = "iqr"), newdata = nd, by = "water_exp_body")


# MST human marker mt model 

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_mst_human_mt_max_s = range(log_mst_human_mt_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(log_mst_human_mt_max_s = seq(-1.695368, 1.567463, by = 0.4), 
            water_exp_body = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

pred <- predictions(m_skin_human_mt, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(mst_human_mt = exp(log_mst_human_mt_max_s*sd(data_follow$log_mst_human_mt_max, na.rm=TRUE) + mean(data_follow$log_mst_human_mt_max, na.rm=TRUE))) 

pred <- pred |> 
  mutate(log_mst_human_mt = log_mst_human_mt_max_s*sd(data_follow$log_mst_human_mt_max, na.rm=TRUE) + mean(data_follow$log_mst_human_mt_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_mst_human_mt, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Log Human Mitochondrial DNA Highest Single Sample",
       y = "Predicted Probability of Skin Infection",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_exp_body)   -> Skin_human_mt

avg_comparisons(m_skin_human_mt, re_formula = NA, variables = list(log_mst_human_mt_max_s = "iqr"), newdata = nd)
avg_comparisons(m_skin_human_mt, re_formula = NA, variables = list(log_mst_human_mt_max_s = "iqr"), newdata = nd, by = "water_exp_body")


### Marginal effects for MST human sewage biomarker model ###

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_mst_human_max_s = range(log_mst_human_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(log_mst_human_max_s = seq(-0.8906186, 2.3137543, by = 0.4), 
            water_exp_body = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

pred <- predictions(m_skin_human, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(log_mst_human = log_mst_human_max_s*sd(data_follow$log_mst_human_max, na.rm=TRUE) + mean(data_follow$log_mst_human_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_mst_human, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Log Human Sewage Biomarker Highest Single Sample",
       y = "Predicted Probability of Skin Infection",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_exp_body) -> Skin_human

avg_comparisons(m_skin_human, re_formula = NA, variables = list(log_mst_human_max_s = "iqr"), newdata = nd)
avg_comparisons(m_skin_human, re_formula = NA, variables = list(log_mst_human_max_s = "iqr"), newdata = nd, by = "water_exp_body")


### MST seagull marker model

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_mst_gull_max_s = range(log_mst_gull_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(log_mst_gull_max_s = seq(-2.718343, 1.832380, by = 0.4), 
            water_exp_body = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

pred <- predictions(m_skin_gull, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(log_mst_gull_max = log_mst_gull_max_s*sd(data_follow$log_mst_gull_max, na.rm=TRUE) + mean(data_follow$log_mst_gull_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_mst_gull_max, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Log MST Seagull Biomarker Highest Single Sample",
       y = "Predicted Probability of Skin Infection",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_exp_body)  -> Skin_gull

avg_comparisons(m_skin_gull, re_formula = NA, variables = list(log_mst_gull_max_s = "iqr"), newdata = nd)
avg_comparisons(m_skin_gull, re_formula = NA, variables = list(log_mst_gull_max_s = "iqr"), newdata = nd, by = "water_exp_body")


### Turbidity model

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_turbidity_s = range(log_turbidity_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(log_turbidity_s = seq(-1.149963, 2.982043, by = 0.4), 
            water_exp_body = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_skin = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            beach_exp_sunscreen = "Yes", beach_exp_repellent = "No",
            sand_contact = "No", household_group = "Yes") 

pred <- predictions(m_skin_turb, re_formula = NA, type = "response", newdata = nd) |> 
  posterior_draws()

pred <- pred |> 
  mutate(log_turbidity = log_turbidity_s*sd(data_follow$log_turbidity_s, na.rm=TRUE) + mean(data_follow$log_turbidity_s, na.rm=TRUE)) 

ggplot(pred, aes(x = log_turbidity, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Log Turbidity",
       y = "Predicted Probability of Skin Infection",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_exp_body)  -> Skin_turbidity

avg_comparisons(m_skin_turb, re_formula = NA, variables = list(log_turbidity_s = "iqr"), newdata = nd)
avg_comparisons(m_skin_turb, re_formula = NA, variables = list(log_turbidity_s = "iqr"), newdata = nd, by = "water_exp_body")

## Combine FIB plots together

Skin_ecoli <- Skin_ecoli + theme(legend.position = "none")
Skin_human <- Skin_human + theme(legend.position = "none")
Skin_human_mt <- Skin_human_mt + theme(legend.position = "none")
Skin_gull <- Skin_gull + theme(legend.position = "none")

Skin_FIB <- Skin_ecoli + Skin_human + Skin_human_mt + Skin_gull + Skin_entero + Skin_turbidity
Skin_FIB + plot_annotation(tag_levels = 'A') + plot_layout(ncol = 2)

remove(Skin_ecoli, Skin_human, Skin_human_mt, Skin_gull, Skin_entero, Skin_turbidity)



## Evaluate E. coli cut-points
# Cut-points of 25th, 50th, 75th & 95th percentiles

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(quantile = scales::percent(c(0.25, 0.5, 0.75, 0.95)),
            e_coli_max = quantile(e_coli_max, na.rm=TRUE, c(0.25, 0.5, 0.75, 0.95)),
            log_e_coli_max_s = quantile(log_e_coli_max_s, na.rm=TRUE, c(0.25, 0.5, 0.75, 0.95)))

list <- data |> distinct(recruit_date, .keep_all = TRUE) |>
  summarize(log_e_coli_max_s = quantile(log_e_coli_max_s, na.rm=TRUE, c(0.25, 0.5, 0.75, 0.95)))
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

mfx <- comparisons(m_skin1, re_formula = NA, variables = "water_exp_body", by = "log_e_coli_max_s", 
                   newdata = nd) |> get_draws()

mfx <- mfx |> mutate(e_coli = exp(log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE))) |> 
  mutate(e_coli = round(e_coli, digits = 0)) |> 
  mutate(draw = draw*1000)

ggplot(mfx, aes(x = draw, y = factor(e_coli), fill = factor(log_e_coli_max_s))) +
  stat_halfeye(slab_alpha = .5)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Risk Difference per 1000 Beachgoers at Specific E. coli Values", y = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete=TRUE) +
  facet_wrap(~ contrast) +
  xlim(-20, 60)

avg_comparisons(m_skin1, re_formula = NA, variables = "water_exp_body", newdata = nd, by = "log_e_coli_max_s")

avg_comparisons(m_skin1, re_formula = NA, variables = "water_exp_body", by = "log_e_coli_max_s",
                newdata = nd, comparison = "lnratioavg", transform = "exp")


# Check proportion of posterior that is greater than 0 and other values

mfx |> group_by(contrast, log_e_coli_max_s) |> 
  summarize(proportion_0 = mean(draw > 0),
            proportion_1 = mean(draw > 1),
            proportion_5 = mean(draw > 5),
            proportion_10 = mean(draw > 10),
            proportion_20 = mean(draw > 20))





### Respiratory illness

data |> distinct(recruit_date, .keep_all = TRUE) |> 
  summarize(log_e_coli_max_s = range(log_e_coli_max_s, na.rm=TRUE))

nd <- data_follow |> 
  data_grid(log_e_coli_max_s = seq(-2.186539, 2.281874, by = 0.2), 
            water_contact = c("No", "Yes"),
            age5 = c("0-9", "10-19", "20+"),
            gender = c("woman/girl", "man/boy", "fluid/trans"),
            education2 = "bachelors", cond_resp = "No", cond_immune = "No",
            cond_allergy = "No", other_rec_act = "Yes", beach_exp_food = "Yes",
            sand_contact = "No", household_group = "Yes") 

pred <- predictions(m_resp3, re_formula = NA, by = c("water_contact", "log_e_coli_max_s"), 
                    type = "response", newdata = nd) |> get_draws()

pred <- pred |> 
  mutate(log_e_coli = log_e_coli_max_s*sd(data_follow$log_e_coli_max, na.rm=TRUE) + mean(data_follow$log_e_coli_max, na.rm=TRUE)) 

ggplot(pred, aes(x = log_e_coli, y = draw)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Log E. coli Highest Single Sample",
       y = "Predicted Probability of Respiratory Illness",
       fill = "") +
  theme_classic() + 
  theme(legend.position = "bottom") +
  facet_wrap(~ water_contact) 

avg_comparisons(m_resp3, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd)
avg_comparisons(m_resp3, re_formula = NA, variables = list(log_e_coli_max_s = "iqr"), newdata = nd, by = "water_contact")



