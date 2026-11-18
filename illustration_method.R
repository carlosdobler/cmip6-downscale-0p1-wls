library(tidyverse)
library(ggforce)

era <-
  read_csv("intermediate_tables/era5_land_annual_mean.csv")

data <-
  read_csv("intermediate_tables/global_mean_temperatures_old.csv") |>
  arrange(Year) |>
  select(-ERA5_Land)

data <-
  left_join(data, era)


data_long <-
  data |>
  pivot_longer(-Year, names_to = "dataset")


data_long |>
  filter(dataset %in% c("DS_Land", "ERA5_Land"), Year >= 1971, Year <= 2010) |>
  ggplot(aes(x = Year, y = value, color = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  labs(title = "Downscaled data aligns with ERA5-Land") +
  ggview::canvas(8, 5)

data_long |>
  filter(dataset %in% c("DS_Land", "GCM_Land"), Year >= 1961) |>
  ggplot(aes(x = Year, y = value, color = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Diff between downscaled data and driving GCM (land only) seems constant"
  ) +
  ggview::canvas(8, 5)

# bias
diff_ds_gcm <- mean(data$GCM_Land - data$DS_Land, na.rm = TRUE) * -1

data |>
  filter(Year >= 1961) |>
  mutate(diff_ds_gcm = DS_Land - GCM_Land) |>
  ggplot(aes(y = 0, x = diff_ds_gcm, color = Year)) +
  geom_sina(size = 2, alpha = 0.8) +
  geom_vline(
    xintercept = diff_ds_gcm,
    linetype = "dashed",
    color = "grey30"
  ) +
  scale_color_viridis_c(limits = c(1961, 2099), breaks = c(1961, 2099)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  # scale_y_continuous(breaks = NULL) +
  labs(
    y = NULL,
    x = "Difference (K)",
    color = "Year",
    title = str_glue("Mean diff = {round(diff_ds_gcm, 3)}")
  ) +
  ggview::canvas(8, 3)

data <-
  data |>
  mutate(
    DS_Land_ext = if_else(is.na(DS_Land), GCM_Land + diff_ds_gcm, DS_Land)
  )

data |>
  select(Year, DS_Land_ext, GCM_Land) |>
  pivot_longer(-Year, names_to = "dataset") |>
  ggplot(aes(x = Year, y = value, color = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  labs(
    title = "A constant difference allows me to extend\nmy downscaled data (only land) back to 1850"
  ) +
  annotate(
    "rect",
    xmin = 1850,
    xmax = 1961,
    ymin = 280.5,
    ymax = 281.8,
    color = "grey20",
    fill = "transparent",
    linetype = "dashed"
  ) +
  ggview::canvas(8, 5)


data_long |>
  filter(dataset %in% c("GCM_Global", "GCM_Land")) |>
  ggplot(aes(x = Year, y = value, color = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Diff between global and land-only GCM data are not constant",
    subtitle = "Land warms faster"
  ) +
  ggview::canvas(8, 5)


# Fit exponential function
# Center years for numerical stability
year_center <- mean(data$Year)
years_c <- data$Year - year_center

diff_land_global <- data$GCM_Global - data$GCM_Land

fit_exp <- nls(
  delta ~ a * exp(b * yr) + c,
  data = data.frame(delta = diff_land_global, yr = years_c),
  start = list(a = -0.5, b = 0.01, c = -1.5) # adjust if nls() fails to converge
)

# summary(fit_exp)

delta_smooth_exp <- predict(fit_exp)

df_delta_exp <-
  tibble(
    year = data$Year,
    delta_raw = diff_land_global,
    delta_smooth_exp = delta_smooth_exp
  )

ggplot(df_delta_exp, aes(x = year)) +
  geom_point(aes(y = delta_raw), color = "grey70") +
  geom_line(
    aes(y = delta_smooth_exp),
    color = "red",
    linewidth = 1
  ) +
  labs(
    y = "K",
    title = "GCM Global - Land",
    subtitle = "Fitted exponential function"
  ) +
  ggview::canvas(8, 5)


data <-
  data |>
  mutate(DS_Global = DS_Land_ext + delta_smooth_exp)

data |>
  select(Year, DS_Global, GCM_Global, DS_Land_ext, GCM_Land) |>
  pivot_longer(-Year, names_to = "dataset") |>
  ggplot(aes(x = Year, y = value, color = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(
    name = NULL,
    values = c(
      "GCM_Global" = "dodgerblue3",
      "GCM_Land" = "tomato",
      "DS_Land_ext" = "orange",
      "DS_Global" = "deepskyblue1"
    )
  ) +
  labs(
    title = "DS_Global = DS_Land_ext with exponential correction"
  ) +
  ggview::canvas(8, 5)

hadcrut <- read_table(
  "https://crudata.uea.ac.uk/cru/data/temperature/HadCRUT5.1Analysis_gl.txt",
  col_names = F,
  na = "-9.999"
)

hadcrut <-
  hadcrut |>
  select(1, 14) |>
  filter(row_number() %% 2 == 1) |>
  rename_with(~ (c("Year", "anom")))

hadcrut_abs <-
  hadcrut |>
  mutate(HadCRUT = anom + 13.974 + 273.15)

left_join(data, hadcrut_abs) |>
  select(Year, DS_Global, HadCRUT) |>
  pivot_longer(-Year, names_to = "dataset") |>
  ggplot(aes(x = Year, y = value, color = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  labs(title = "Global downscaled data aligns with observational", y = "K") +
  ggview::canvas(8, 5)


pre_ind <-
  data |>
  filter(Year >= 1850, Year <= 1900) |>
  pull(DS_Global) |>
  mean()

data <-
  data |>
  mutate(DS_Global_anom = DS_Global - pre_ind) |>
  mutate(
    DS_Global_anom_w = slider::slide_dbl(
      DS_Global_anom,
      .f = mean,
      .before = 10,
      .after = 10,
      .complete = T
    )
  )


first_year_exceeding <- function(data, threshold) {
  data |>
    filter(DS_Global_anom_w >= threshold) |>
    slice_min(Year, n = 1) |>
    pull(Year)
}

tb_yrs <-
  tibble(wl = seq(0.5, 3, 0.5)) |>
  mutate(year = map_dbl(wl, first_year_exceeding, data = data))


data |>
  select(Year, DS_Global_anom_w) |>
  ggplot(aes(x = Year, y = DS_Global_anom_w)) +
  geom_line() +
  # geom_point() +
  geom_point(
    data = data |> semi_join(tb_yrs, by = c("Year" = "year")),
    size = 3,
    color = "red"
  ) +
  geom_text(
    data = data |> semi_join(tb_yrs, by = c("Year" = "year")),
    aes(label = Year),
    hjust = -0.4,
    vjust = 1.4,
    size = 4,
    color = "red"
  ) +
  labs(y = "K (anomaly from baseline)") +
  ggview::canvas(7, 5)
