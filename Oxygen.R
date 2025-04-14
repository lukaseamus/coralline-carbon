# 1. Load data ####
# Load output from PreSens Measurement Studio 2
require(tidyverse)
require(here)
files <- here("Oxygen") %>% list.files(pattern = "\\.csv$", full.names = TRUE)

O2 <- files %>%
  map(~ read.csv(.x, skip = 1, header = TRUE) %>%
        drop_na(Value) %>%
        mutate(delta_t = delta_t %>% as.numeric(),
               delta_t_c = delta_t - mean(delta_t), # centre time for linear models
               Date = Date %>% str_c(Time, sep = " ") %>% mdy_hms(),
               Calibration_Date = Calibration_Date %>% # parsing date differently due to varying formats
                 parse_date_time(orders = c("mdy HMS", "mdy HM")))
      ) %>%
  set_names(str_remove(basename(files), "\\.csv$") %>% make.names) %>%
  imap(~ .x %>% mutate(Date_true = .y %>% 
                         str_split(pattern = "_", n = 2) %>% 
                         map_chr(1) %>% str_sub(start = 2) %>%
                         ymd(),
                       ID = .y %>% 
                         str_split(pattern = "_", n = 2) %>% 
                         map_chr(2) %>% fct())
       )

str(O2$X240521_A_1_D) # examples
str(O2$X250121_A_14_L)

# Compare time to last calibration
O2 %>%
  map(~ .x %>% 
        mutate(Time = Calibration_Date %--% Date %>% 
                 as.duration() / ddays())
      ) %>%
  bind_rows() %>%
  group_by(ID) %>%
  summarise(Time = mean(Time)) %>%
  arrange(Time) %>%
  print(n = 136)
# Never more than 45 days between calibration and measurement; 
# mostly done on the same day.

# Compare internal and true dates
O2 %>%
  map(~ .x %>% 
        mutate(Time = Date %--% Date_true %>% 
                 as.duration() / ddays())
  ) %>%
  bind_rows() %>%
  group_by(ID) %>%
  summarise(Time = mean(Time)) %>%
  arrange(Time) %>%
  print(n = 136)
# Internal date never matches true date, due to disconnection from internet

# 2. Raw models ####
# 2.1 Visualise data ####
require(patchwork)
O2 %>%
  imap(~ .x %>%
        ggplot() +
          geom_hline(yintercept = 0) +
          geom_point(aes(delta_t, Value), shape = 16, alpha = 0.2) +
          theme_minimal() +
          theme(panel.grid = element_blank()) +
          ggtitle(.y)
       ) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_data.pdf", path = "Plots", 
         width = 100, height = 50, unit = "cm", device = cairo_pdf)
# First there are some negative values which are impossible. Then there are some 
# completely nonsensical measurement series due to miscalibration. These are
# M_5_D, M_5_L, M_6_D, M_6_L, B_7_D, B_7_L, J_7_D, J_7_L.

# Another way to check this is to filter for values < 0 or > 1000 µM, which
# are (practically) impossible:
O2 %>%
  map(~ .x %>%
        filter(Value > 1e3 | Value < 0)) %>%
  bind_rows() %>%
  group_by(ID) %>%
  summarise(min = min(Value),
            max = max(Value),
            n = length(Value))
# A few other measurement series also get flagged because they had substantial 
# measurement error with some outliers above 1000 µM. M_1_L needs it's negative
# outlier removed but with A_1_L, J_1_L, A_4_L, J_4_L, M_4_L, A_6_L, M_7_L, and A_9_L 
# it has the right shape according to the plot, so can be retained.

# 2.2 Clean data ####
require(magrittr)
O2 %<>%
  map(~ .x %>%
        filter(!ID %in% c("M_5_D", "M_5_L", "M_6_D", "M_6_L", # filter out bad measurements
                          "B_7_D", "B_7_L", "J_7_D", "J_7_L") &
                 Value >= 0)) %>% # filter out negative values
  keep(~ nrow(.x) > 0)

# Visualise clean data
O2 %>%
  imap(~ .x %>%
        ggplot() +
          geom_hline(yintercept = 0) +
          geom_point(aes(delta_t, Value), shape = 16, alpha = 0.2) +
          theme_minimal() +
          theme(panel.grid = element_blank()) +
          ggtitle(.y)
       ) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_data_cleaned.pdf", path = "Plots", 
         width = 100, height = 50, unit = "cm", device = cairo_pdf)

# The biggest problem about this is that two of these measurements are blanks, 
# which are required to correct all other measurements in the same round. To 
# avoid losing the two good samples in round 7, blank must be estimated from 
# the other blanks via multilevel modelling. To do that blanks must be modelled
# separately in a single dataframe rather than a list of dataframes.

O2_B_L <- O2 %>% 
  map(~ .x %>% filter(ID %>% str_detect("B") & # select blanks
                        ID %>% str_detect("L"))) %>% # select light
  keep(~ nrow(.x) > 0)

O2_B_D <- O2 %>% 
  map(~ .x %>% filter(ID %>% str_detect("B") & # select blanks
                        ID %>% str_detect("D"))) %>% # select dark
  keep(~ nrow(.x) > 0)

# 2.3 Model selection ####
# O2 measurements in light often reached saturation (~1000 µM) but I want to retain the 
# full measurement series (setting a limit of linearity is arbitrary). Dark data are on
# the other hand are always linear because we never allowed O2 to reach zero. So different
# models needs to be fit to light and dark data. The easiest way is to split the data.

O2_L <- O2 %>% 
  map(~ .x %>% filter(ID %>% str_detect("L") & # select light
                        !ID %>% str_detect("B"))) %>% # remove blanks
  keep(~ nrow(.x) > 0)
  
O2_D <- O2 %>% 
  map(~ .x %>% filter(ID %>% str_detect("D") & # select dark
                        !ID %>% str_detect("B"))) %>% # remove blanks
  keep(~ nrow(.x) > 0)

# There are a variety of saturating models to choose from that describe an increase towards 
# a plateau with only two parameters: beta (linear slope below the limit of linearity) and 
# O2_max (maximum O2 concentration). 

# The simplest are the rectangular hyperbola (i.e. Michaelis-Menten function)
# O2 = O2_max * beta * t / ( O2_max + beta * c )
# the exponential saturation function (i.e. exponential cumulative distribution function)
# O2 = O2_max * ( 1 - exp( -beta * t / O2_max ) )
# and the hyperbolic tangent
# O2 = O2_max * tanh( beta * t / O2_max )

# However, the water in the incubation medium already has O2, so an intercept term must be added.
# When the intercept (O2_0) is added the interpretation of O2_max becomes the O2 differential 
# between intercept and saturating level, i.e. the maximal O2 addition above baseline. 
# The functions thus become:

# O2 = O2_0 + O2_max * beta * t / ( O2_max + beta * c )
# O2 = O2_0 + O2_max * ( 1 - exp( -beta * t / O2_max ) )
# O2 = O2_0 + O2_max * tanh( beta * t / O2_max )

# I instinctively favour the last function because it "sticks" to the linear phase for longest:

# example parameter values
O2_0 <- 230 # rounded regional empirical mean (228 µM)
O2_max <- 1e3 - O2_0 # the usual supersaturation limit (1000 µM) minus the intercept
beta <- 13 # no usable literature data (assuming saturation in one hour as O2_max / 60)

tibble(t = seq(0, 120)) %>% # each incubation ran for two hours
  mutate(lm = O2_0 + beta * t, # linear model for comparison
         rh = O2_0 + O2_max * beta * t / ( O2_max + beta * t ),
         es = O2_0 + O2_max * ( 1 - exp( -beta * t / O2_max ) ),
         ht = O2_0 + O2_max * tanh( beta * t / O2_max )) %>%
  pivot_longer(cols = -t, names_to = "Function", values_to = "O2") %>%
  ggplot(aes(t, O2, colour = Function)) +
    geom_hline(yintercept = 1e3) +
    geom_line() +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    coord_cartesian(ylim = c(0, 1e3))

# However, the goodness of fit of the three functions must be compared.

# 2.4 Hyperbolic tangent ####
# 2.4.1 Prior simulation ####
tibble(n = 1:1e3,
       # gamma distribution is reparameterised in terms of mean and sd
       O2_max = rgamma(n = 1e3, shape = O2_max^2 / 100^2, rate = O2_max / 100^2), 
       # in this model the slope has to be greater than zero
       beta = rgamma(n = 1e3, shape = beta^2 / 6^2, rate = beta / 6^2),
       O2_0 = rgamma(n = 1e3, shape = O2_0^2 / 20^2, rate = O2_0 / 20^2)) %>%
  expand_grid(t = seq(0, 120)) %>%
  mutate(O2 = O2_0 + O2_max * tanh( beta * t / O2_max )) %>%
  ggplot(aes(t, O2, group = n)) +
    geom_hline(yintercept = c(0, 1e3)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks like plenty of variability.

# 2.4.2 Run model ####
O2_L_ht_stan <- "
data{
  int n;
  vector<lower=0>[n] Value;
  vector<lower=0>[n] delta_t;
}

parameters{
  real<lower=0> O2_max;
  real<lower=0> beta;
  real<lower=0> O2_0;
  real<lower=0> sigma;
}

model{
  // Priors
  O2_0 ~ gamma( 230^2 / 20^2, 230 / 20^2 ); // reparameterised with mean and sd
  beta ~ gamma( 13^2 / 6^2, 13 / 6^2 );
  O2_max ~ gamma( 770^2 / 100^2, 770 / 100^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = O2_0 + O2_max * tanh( beta * delta_t[i] / O2_max );
  }

  // Likelihood
  Value ~ normal( mu , sigma );
}
"

require(cmdstanr)
O2_L_ht_mod <- O2_L_ht_stan %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
O2_L_ht_samples <- O2_L %>%
  map(~ O2_L_ht_mod$sample(
    data = .x %>% 
      select(Value, delta_t) %>%
      compose_data(), 
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4, 
    iter_sampling = 1e4))

# 2.4.3 Model checks ####
# check Rhat, effective sample size and chains
O2_L_ht_samples %>%
  map(~ .x$summary()) %>%
  bind_rows() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

require(bayesplot)
ggsave(
  O2_L_ht_samples %>%
    map(~ .x$draws(format = "df") %>%
          mcmc_rank_overlay()) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_L_ht_rank.pdf", path = "Plots", 
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

O2_L_ht_samples %>%
  map(~ .x$draws(format = "df") %>%
        mcmc_pairs(pars = c("O2_max", "beta", "O2_0"))) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_L_ht_pairs.pdf", path = "Plots", 
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# some correlation between Fmax and beta, indicating some interdependence

# 2.4.4 Prior-posterior comparison ####
source("functions.R")
# sample prior
O2_L_ht_prior <- O2_L %>%
  map(~ prior_samples(model = O2_L_ht_mod,
                      data = .x %>% 
                        select(Value, delta_t) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
kinetics_standard_ht_prior %>%
  map2(kinetics_standard_ht_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA), # no groups so this has to be an empty list or tibble
                               parameters = c("Fmax", "beta", "F0", "sigma"),
                               format = "long")) %>%
  map(~ .x %>% prior_posterior_plot()) %>%
  wrap_plots() +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# some posteriors for F0 broke out of the expected prior probability space

# 2.4.5 Predictions ####
kinetics_standard_ht_predictions <- kinetics_standard_ht_prior %>%
  map2(kinetics_standard_ht_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("Fmax", "beta", "F0", "sigma"),
                               format = "short")) %>%
  map2(kinetics, ~ spread_continuous(.x, .y$standard, predictor_name = "Concentration")) %>%
  map(~ .x %>% mutate(mu = Fmax * tanh( beta * Concentration / Fmax ) + F0,
                      obs = rnorm(n(), mu, sigma)))

kinetics_standard_ht_predictions_summary <- kinetics_standard_ht_predictions %>%
  map(~ .x %>% group_by(distribution, Concentration) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_")) 

kinetics_standard_ht_predictions_summary %>%
  map2(kinetics,
       ~ ggplot() +
            geom_point(data = .y$standard, aes(Concentration, Fluorescence)) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(Concentration, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(Concentration, ymin = obs_ymin, ymax = obs_ymax,
                            alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(Concentration, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
            # + coord_cartesian(xlim = c(0, 1), ylim = c(0, 5e3)) # unhash to check F0
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots()
# fit doesn't look as good as expected
# try exponential saturation function 

# 2.5 Exponential saturation ####



# 2.6 Rectangular hyperbola ####


# 2.7 Dark samples ####


# 2.8 Blanks ####




# 2.7 Data preparation ####

# 3. Conversion ####

# 4. Rate models ####

