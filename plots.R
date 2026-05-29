library(tidyverse)

data <-
  read_csv("global_mean_temperatures.csv") |>
  arrange(Year)


ggplot(data, aes(x = Year)) +
  geom_line(aes(y = GCM_Global), color = "blue") +
  geom_line(aes(y = GCM_Land), color = "red")


ggplot(data |> filter(Year >= 1961), aes(x = Year)) +
  geom_line(aes(y = DS_Land), color = "blue") +
  geom_line(aes(y = GCM_Land), color = "red")

data |>
  filter(Year >= 1961) |>
  mutate(dif = GCM_Land - DS_Land) |>
  ggplot(aes(x = Year, y = dif)) +
  geom_line()

data |>
  filter(Year >= 1961) |>
  mutate(dif = GCM_Global - GCM_Land) |>
  ggplot(aes(x = Year, y = dif)) +
  geom_line() +
  geom_smooth()


C <- data$GCM_Land - data$DS_Land
mean(C, na.rm = T)
sd(C, na.rm = T)


# Center years for numerical stability (important for high-degree polynomials)
year_center <- mean(data$Year)
years_c <- data$Year - year_center

# Fit polynomials of increasing degree and compare
delta <- data$GCM_Global - data$GCM_Land
fits <- lapply(1:5, function(d) lm(delta ~ poly(years_c, d, raw = TRUE)))

# Check R² and AIC for each degree
model_summary <- data.frame(
  degree = 1:5,
  r2 = sapply(fits, function(m) summary(m)$r.squared),
  AIC = sapply(fits, AIC),
  adj_r2 = sapply(fits, function(m) summary(m)$adj.r.squared)
)
print(model_summary)

# Example: use degree 3 (adjust based on model_summary above)
best_degree <- 2
fit_poly <- lm(delta ~ poly(years_c, best_degree, raw = TRUE))

# Generate smoothed delta
delta_smooth <- predict(fit_poly, newdata = data.frame(years_c = years_c))

# Build a data frame with both raw and smoothed delta
df_delta <- data.frame(
  year = data$Year,
  delta_raw = delta,
  delta_smooth = delta_smooth
)

ggplot(df_delta, aes(x = year)) +
  geom_point(aes(y = delta_raw), colour = "grey60", size = 1.2, alpha = 0.6) +
  geom_line(aes(y = delta_smooth), colour = "steelblue", linewidth = 1)


# Good starting values matter a lot for nls() — eyeball your delta plot
# to get rough guesses before running this
fit_exp <- nls(
  delta ~ a * exp(b * years_c) + c,
  data = data.frame(delta = delta, years_c = years_c),
  start = list(a = -0.5, b = 0.01, c = -1.5) # adjust if nls() fails to converge
)

summary(fit_exp)

delta_smooth_exp <- predict(fit_exp)

df_delta_exp <- data.frame(
  year = data$Year,
  delta_raw = delta,
  delta_smooth_exp = delta_smooth_exp
)

ggplot(df_delta_exp, aes(x = year)) +
  geom_point(aes(y = delta_raw, colour = "Raw δ(t)"), size = 1.2, alpha = 0.6) +
  geom_line(
    aes(y = delta_smooth_exp, colour = "Exponential fit"),
    linewidth = 1
  ) +
  scale_colour_manual(
    name = NULL,
    values = c("Raw δ(t)" = "grey50", "Exponential fit" = "tomato")
  )


data$DS_Global <- data$DS_Land + delta_smooth_exp

# Pivot to long format for ggplot
df_long <- pivot_longer(
  data,
  cols = -Year,
  names_to = "series",
  values_to = "temp"
) |>
  filter(Year >= 1961)

# Panel 1: all trajectories
ggplot(
  df_long,
  aes(x = Year, y = temp, colour = series)
) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(
    name = NULL,
    values = c(
      "GCM_Global" = "black",
      "GCM_Land" = "tomato",
      "DS_Land" = "orange",
      "DS_Global" = "steelblue"
    )
  ) +
  theme(legend.position = "bottom")

# Panel 2: residuals of the exponential fit
df_resid <- data.frame(
  year = data$Year,
  resid = delta - delta_smooth_exp
)

ggplot(df_resid, aes(x = year, y = resid)) +
  geom_point(colour = "grey40", size = 1.2, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  geom_smooth(
    method = "loess",
    se = FALSE
  )


gistemp <- read_csv(
  "https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.csv",
  skip = 1,
  na = c("***", "****")
)
