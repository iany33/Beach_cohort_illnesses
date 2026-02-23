
pacman::p_load(
  rio,
  here,
  tidyverse,
  gtsummary,
  rstatix, 
  janitor, 
  flextable,
  viridis
)

# Load dataset

data <- import(here("data.xlsx"))

# Create new variables

data <- data |> 
  mutate(water_contact2 = case_when(
    water_contact == "No" ~ "No contact",
    water_exp_mouth == "Yes" ~ "Swallowed water",
    water_exp_body == "Yes" ~ "Body immersion",
    TRUE ~ "Minimal contact")) |>
  mutate(water_contact2 = as.factor(water_contact2)) |> 
  mutate(water_contact2 = fct_relevel(water_contact2, "No contact", "Minimal contact")) 

data <- data |> 
  mutate(water_contact3 = factor(water_contact2, ordered = T, 
                                 levels = c("No contact", "Minimal contact", "Body immersion", "Swallowed water")))

data <- data |> 
  mutate(water_contact4 = case_when(
    water_contact == "No" ~ "No contact",
    water_exp_body == "Yes" ~ "Body immersion",
    TRUE ~ "Minimal contact")) |>
  mutate(water_contact4 = as.factor(water_contact4)) |> 
  mutate(water_contact4 = fct_relevel(water_contact4, "No contact", "Minimal contact")) 

data<- data |> 
  mutate(water_contact5 = factor(water_contact4, ordered = T, 
                                 levels = c("No contact", "Minimal contact", "Body immersion")))

data <- data |> 
  mutate(age5 = case_when(
    (age4 == "0-4" | age4 == "5-9")  ~ "0-9",
    (age4 == "10-14" | age4 == "15-19") ~ "10-19",
    TRUE ~ "20+")) |>
  mutate(age5 = as.factor(age5)) |> 
  mutate(age5 = fct_relevel(age5, "0-9", "10-19", "20+")) 

data <- data |> 
  mutate(ethnicity2 = case_when(
    ethnicity == "White"  ~ "White",
    TRUE ~ "Non-White")) |>
  mutate(ethnicity2 = as.factor(ethnicity2)) |> 
  mutate(ethnicity2 = fct_relevel(ethnicity2, "White", "Non-White")) 

data <- data |> 
  mutate(water_time = replace_na(water_time, 0)) |> 
  mutate(water_time_s = (water_time - mean(water_time, na.rm = TRUE)) / sd(water_time, na.rm = TRUE))

data |> group_by(date) |> ggplot(aes(x = water_time)) + geom_histogram()
data |> group_by(date) |> ggplot(aes(x = water_time_s)) + geom_histogram()

data_follow <- data |> filter(follow == "Yes") 


# Creating new variables

data |> 
  select(age4, gender, education2, ethnicity, water_contact2, follow) |> 
  tbl_summary(by = follow, digits = list(all_categorical() ~ c(0, 1))) |> 
  add_overall() |>
  as_flex_table() 

# Descriptive tables

data |> 
  select(other_rec_act, beach_exp_food, sand_contact, household_group, follow) |> 
  tbl_summary(by = follow, digits = list(all_categorical() ~ c(0, 1))) |> 
  add_overall() |>
  as_flex_table() 

data_follow |> tabyl(respiratory3)
data_follow |> tabyl(skin_infection3)
data_follow |> tabyl(ear_infection3)
data_follow |> tabyl(eye_infection3)

data_follow |> tabyl(respiratory_3day)
data_follow |> tabyl(skin_infection_3day)
data_follow |> tabyl(ear_infection_3day)
data_follow |> tabyl(eye_infection_3day)


data_follow |> group_by(house_id) |> 
  mutate(n = n(), house_size = case_when(
    n == 1 ~ 1, n == 2 ~ 2, n == 3 ~ 3,
    n == 4 ~ 4, n == 5 ~ 5, n == 6 ~ 6)) |> 
  tabyl(house_size)

# Examine water contact and illness outcome relationships

data_follow |> 
  select(respiratory3, skin_infection3, ear_infection3, eye_infection3, water_contact2) |> 
  tbl_summary(by = water_contact2, digits = list(all_categorical() ~ c(0, 1)),
              type = all_categorical() ~ "categorical")

data_follow |> 
  select(respiratory3, skin_infection3, ear_infection3, eye_infection3, water_contact4) |> 
  tbl_summary(by = water_contact4, digits = list(all_categorical() ~ c(0, 1)),
              type = all_categorical() ~ "categorical")

data_follow |> 
  select(respiratory3, skin_infection3, ear_infection3, eye_infection3, water_exp_body) |> 
  tbl_summary(by = water_exp_body, digits = list(all_categorical() ~ c(0, 1)),
              type = all_categorical() ~ "categorical")

data_follow |> 
  select(respiratory3, skin_infection3, ear_infection3, eye_infection3, water_contact) |> 
  tbl_summary(by = water_contact, digits = list(all_categorical() ~ c(0, 1)),
              type = all_categorical() ~ "categorical")

data_follow |> 
  ggplot(aes(x = respiratory3, y = log_e_coli_max, fill = respiratory3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Respiratory illness",
       y = "Log E. coli highest single sample") +
  facet_grid(~water_contact4)

data_follow |> 
  ggplot(aes(x = skin_infection3, y = log_e_coli_max, fill = skin_infection3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Skin illness",
       y = "Log E. coli highest single sample") +
  facet_grid(~water_contact4)

data_follow |> 
  ggplot(aes(x = ear_infection3, y = log_e_coli_max, fill = ear_infection3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Ear infection",
       y = "Log E. coli highest single sample") +
  facet_grid(~water_contact4)

data_follow |> 
  ggplot(aes(x = eye_infection3, y = log_e_coli_max, fill = eye_infection3)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "cividis") +
  labs(x = "Eye infection",
       y = "Log E. coli highest single sample") +
  facet_grid(~water_contact4)


