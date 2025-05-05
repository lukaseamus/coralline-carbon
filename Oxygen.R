# 1. Load data ####
# Load output from PreSens Measurement Studio 2
require(tidyverse)
require(here)
files <- here("Oxygen") %>% list.files(pattern = "\\.csv$", full.names = TRUE)

O2 <- files %>%
  map(~ .x %>% read.csv(skip = 1, header = TRUE) %>%
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

# Finally, there are sometimes interruptions in the data stream. Rounds 6 and 14 started, 
# then stopped for 3 and 16 min, then restarted and in total ran for longer than 120 min. 
# The first measurements are likely a mistake, i.e. taken before the actual start of the
# incubation and can be dropped. However, in round 7 the interruption happens 8 min in
# and only lasts for 4 min. Even though the incubation runs longer, it is here not possible
# to remove the substantial first part.

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
# the other hand are linear until they lateau towards zero. So different models needs to 
# be fit to light and dark data. The easiest way is to split the data.

O2_L <- O2 %>% 
  map(~ .x %>% filter(ID %>% str_detect("L") & # select light
                        !ID %>% str_detect("B"))) %>% # remove blanks
  keep(~ nrow(.x) > 0)
  
O2_D <- O2 %>% 
  map(~ .x %>% filter(ID %>% str_detect("D") & # select dark
                        !ID %>% str_detect("B"))) %>% # remove blanks
  keep(~ nrow(.x) > 0)

# There are a variety of saturating models to choose from that describe an increase towards 
# a plateau in light incubations with only two parameters: beta (linear slope below the limit 
# of linearity) and O2_max (maximum O2 concentration). 

# The simplest forms from which most others are derived are 
# 1. The piecewise linear function
# O2 = {
#       beta * t, t <= O2_max / beta
#       O2_max, t > O2_max / beta 
#      }
# 2. The rectangular hyperbola (i.e. Michaelis-Menten function)
# O2 = O2_max * beta * t / ( O2_max + beta * t )
# 3. The exponential saturation function (i.e. exponential cumulative distribution function)
# O2 = O2_max * ( 1 - exp( -beta * t / O2_max ) )
# 4. The hyperbolic tangent
# O2 = O2_max * tanh( beta * t / O2_max )
# See e.g. Jassby & Platt (1976) Limnology and Oceanography for details.

# However, the water in the incubation medium already has O2, so an intercept term must be added.
# When the intercept (O2_0) is added the interpretation of O2_max becomes the maximal O2 differential 
# between intercept and saturating level, i.e. the maximal O2 addition above baseline. 
# The functions thus become:

# O2 = {
#       O2_0 + beta * t, t <= O2_max / beta
#       O2_0 + O2_max, t > O2_max / beta 
#      }
# O2 = O2_0 + O2_max * beta * t / ( O2_max + beta * t )
# O2 = O2_0 + O2_max * ( 1 - exp( -beta * t / O2_max ) )
# O2 = O2_0 + O2_max * tanh( beta * t / O2_max )

# I instinctively favour the last function because it "sticks" to the linear phase for longest:

# example parameter values
O2_0 <- 230 # rounded regional empirical mean (228 µM)
O2_max <- 600 # the usual supersaturation limit is below 1000 µM, so O2_max should be below 770 µM
beta <- 8 # no usable literature data (assuming saturation in one hour as O2_max / 60)

tibble(t = seq(0, 120)) %>% # each incubation ran for two hours
  mutate(lm = O2_0 + beta * t, # linear model for comparison
         pl = if_else(t <= O2_max / beta, O2_0 + beta * t, O2_0 + O2_max),
         rh = O2_0 + O2_max * beta * t / ( O2_max + beta * t ),
         es = O2_0 + O2_max * ( 1 - exp( -beta * t / O2_max ) ),
         ht = O2_0 + O2_max * tanh( beta * t / O2_max )) %>%
  pivot_longer(cols = -t, names_to = "Function", values_to = "O2") %>%
  ggplot(aes(t, O2, colour = Function)) +
    geom_hline(yintercept = 1e3) + # generally observed supersaturation limit
    geom_line() +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    coord_cartesian(ylim = c(0, 1e3))

# However, the goodness of fit of the three functions must be compared.

# Now for the dark incubations I'll try the same three functions but modified so
# their slopes are negative and they approach zero. In all cases this means we
# can lose O2_max because the differential between the intercept and plateau is
# always going to be equivalent to the intercept. The modified functions are:

# O2 = {
#       O2_0 - beta * t, t <= O2_0 / beta
#       0, t > O2_0 / beta 
#      }
# O2 = O2_0^2 / ( O2_0 + beta * t )
# O2 = O2_0 * exp( -beta * t / O2_0 )
# O2 = O2_0 * (1 - tanh( beta * t / O2_0 ) )

# the oxygen was practically never fully consumed, so beta must be around O2_0 / 120
# or lower, so around 2
beta_r <- 2 

tibble(t = seq(0, 120)) %>% # each incubation ran for two hours
  mutate(lm = O2_0 - beta_r * t, # linear model for comparison
         pl = if_else(t <= O2_0 / beta_r, O2_0 - beta_r * t, 0),
         rh = O2_0^2 / ( O2_0 + beta_r * t ),
         es = O2_0 * exp( -beta_r * t / O2_0 ),
         ht = O2_0 * ( 1 - tanh( beta_r * t / O2_0 ) )) %>%
  pivot_longer(cols = -t, names_to = "Function", values_to = "O2") %>%
  ggplot(aes(t, O2, colour = Function)) +
    geom_hline(yintercept = 0) +
    geom_line() +
    theme_minimal() +
    theme(panel.grid = element_blank())

# Again, I intuitively think the hyperbolic tangent will be best, but we'll see. 

# 2.4 Photosynthesis: hyperbolic tangent ####
# 2.4.1 Prior simulation ####
tibble(n = 1:1e3,
       # gamma distribution is reparameterised in terms of mean and sd
       O2_max = rgamma(n = 1e3, shape = O2_max^2 / 100^2, rate = O2_max / 100^2),
       # in this model the slope has to be greater than zero
       beta = rgamma(n = 1e3, shape = beta^2 / 5^2, rate = beta / 5^2),
       O2_0 = rgamma(n = 1e3, shape = O2_0^2 / 50^2, rate = O2_0 / 50^2)) %>%
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
  O2_0 ~ gamma( 230^2 / 50^2 , 230 / 50^2 ); // reparameterised with mean and sd
  beta ~ gamma( 8^2 / 5^2 , 8 / 5^2 );
  O2_max ~ gamma( 600^2 / 100^2 , 600 / 100^2 );
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
ggsave(
  O2_L_ht_prior %>%
    map2(O2_L_ht_samples,
         ~ prior_posterior_draws(prior_samples = .x,
                                 posterior_samples = .y,
                                 group = list(NA), # no groups so this has to be an empty list or tibble
                                 parameters = c("O2_max", "beta", "O2_0", "sigma"),
                                 format = "long")) %>%
    imap(~ .x %>% prior_posterior_plot() +
                    ggtitle(.y)) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_L_ht_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# Some posteriors for O2_0 are pretty unconstrained and O2_max sometimes broke out
# of the expected prior probability space, but generally priors are not restrictive.

# 2.4.5 Predictions ####
O2_L_ht_predictions <- O2_L_ht_prior %>%
  map2(O2_L_ht_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("O2_max", "beta", "O2_0", "sigma"),
                               format = "short")) %>%
  map2(O2_L, ~ spread_continuous(.x, .y,
                                 predictor_name = "delta_t",
                                 length = 50)) %>%
  map(~ .x %>% mutate(mu = O2_0 + O2_max * tanh( beta * delta_t / O2_max ),
                      obs = rnorm( n(), mu, sigma )))

O2_L_ht_predictions_summary <- O2_L_ht_predictions %>%
  map(~ .x %>% group_by(distribution, delta_t) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_"))

rm(O2_L_ht_predictions) # remove raw predictions because they are too big

O2_L_ht_predictions_summary %>%
  map2(O2_L,
       ~ ggplot() +
            geom_point(data = .y, aes(delta_t, Value),
                       shape = 16, alpha = 0.05) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(delta_t, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            # geom_ribbon(data = .x %>% filter(distribution == "posterior"), # unhash to check
            #             aes(delta_t, ymin = obs_ymin, ymax = obs_ymax, # predicted observations
            #                 alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_L_ht_prediction.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# The hyperbolic tangent fits well, except for cases where the transition from linear
# increase to plateau is very abrupt, e.g. M_4_L, or where the plateau occurs early 
# and the intercept consequently is not as well constrained as the plateau, e.g. J_6_L.

# 2.5 Photosynthesis: exponential saturation ####
# 2.5.1 Prior simulation ####
tibble(n = 1:1e3,
       O2_max = rgamma(n = 1e3, shape = O2_max^2 / 100^2, rate = O2_max / 100^2),
       beta = rgamma(n = 1e3, shape = beta^2 / 5^2, rate = beta / 5^2),
       O2_0 = rgamma(n = 1e3, shape = O2_0^2 / 50^2, rate = O2_0 / 50^2)) %>%
  expand_grid(t = seq(0, 120)) %>%
  mutate(O2 = O2_0 + O2_max * ( 1 - exp( -beta * t / O2_max ) )) %>%
  ggplot(aes(t, O2, group = n)) +
    geom_hline(yintercept = c(0, 1e3)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.5.2 Run model ####
O2_L_es_stan <- "
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
  O2_0 ~ gamma( 230^2 / 50^2 , 230 / 50^2 ); // reparameterised with mean and sd
  beta ~ gamma( 8^2 / 5^2 , 8 / 5^2 );
  O2_max ~ gamma( 600^2 / 100^2 , 600 / 100^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = O2_0 + O2_max * ( 1 - exp( -beta * delta_t[i] / O2_max ) );
  }

  // Likelihood
  Value ~ normal( mu , sigma );
}
"

O2_L_es_mod <- O2_L_es_stan %>%
  write_stan_file() %>%
  cmdstan_model()

O2_L_es_samples <- O2_L %>%
  map(~ O2_L_es_mod$sample(
    data = .x %>%
      select(Value, delta_t) %>%
      compose_data(),
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4,
    iter_sampling = 1e4))

# 2.5.3 Model checks ####
# check Rhat, effective sample size and chains
O2_L_es_samples %>%
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

ggsave(
  O2_L_es_samples %>%
    map(~ .x$draws(format = "df") %>%
          mcmc_rank_overlay()) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_L_es_rank.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

# 2.5.4 Prior-posterior comparison ####
# sample prior
O2_L_es_prior <- O2_L %>%
  map(~ prior_samples(model = O2_L_es_mod,
                      data = .x %>%
                        select(Value, delta_t) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
ggsave(
  O2_L_es_prior %>%
    map2(O2_L_es_samples,
         ~ prior_posterior_draws(prior_samples = .x,
                                 posterior_samples = .y,
                                 group = list(NA), # no groups so this has to be an empty list or tibble
                                 parameters = c("O2_max", "beta", "O2_0", "sigma"),
                                 format = "long")) %>%
    imap(~ .x %>% prior_posterior_plot() +
                    ggtitle(.y)) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_L_es_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# Some posteriors for O2_0 are pretty unconstrained and O2_max sometimes broke out
# of the expected prior probability space, but generally priors are not restrictive.

# 2.5.5 Predictions ####
O2_L_es_predictions <- O2_L_es_prior %>%
  map2(O2_L_es_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("O2_max", "beta", "O2_0", "sigma"),
                               format = "short")) %>%
  map2(O2_L, ~ spread_continuous(.x, .y,
                                 predictor_name = "delta_t",
                                 length = 50)) %>%
  map(~ .x %>% mutate(mu = O2_0 + O2_max * ( 1 - exp( -beta * delta_t / O2_max ) ),
                      obs = rnorm( n(), mu, sigma )))

O2_L_es_predictions_summary <- O2_L_es_predictions %>%
  map(~ .x %>% group_by(distribution, delta_t) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_"))

rm(O2_L_es_predictions) # remove raw predictions because they are too big

O2_L_es_predictions_summary %>%
  map2(O2_L,
       ~ ggplot() +
            geom_point(data = .y, aes(delta_t, Value),
                       shape = 16, alpha = 0.05) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(delta_t, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            # geom_ribbon(data = .x %>% filter(distribution == "posterior"), # unhash to check
            #             aes(delta_t, ymin = obs_ymin, ymax = obs_ymax, # predicted observations
            #                 alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_L_es_prediction.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# The exponential saturation function produces a worse fit.

# 2.6 Photosynthesis: rectangular hyperbola ####
# 2.6.1 Prior simulation ####
tibble(n = 1:1e3,
       O2_max = rgamma(n = 1e3, shape = O2_max^2 / 100^2, rate = O2_max / 100^2),
       beta = rgamma(n = 1e3, shape = beta^2 / 5^2, rate = beta / 5^2),
       O2_0 = rgamma(n = 1e3, shape = O2_0^2 / 50^2, rate = O2_0 / 50^2)) %>%
  expand_grid(t = seq(0, 120)) %>%
  mutate(O2 = O2_0 + O2_max * beta * t / ( O2_max + beta * t )) %>%
  ggplot(aes(t, O2, group = n)) +
    geom_hline(yintercept = c(0, 1e3)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.6.2 Run model ####
O2_L_rh_stan <- "
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
  O2_0 ~ gamma( 230^2 / 50^2 , 230 / 50^2 ); // reparameterised with mean and sd
  beta ~ gamma( 8^2 / 5^2 , 8 / 5^2 );
  O2_max ~ gamma( 600^2 / 100^2 , 600 / 100^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = O2_0 + O2_max * beta * delta_t[i] / ( O2_max + beta * delta_t[i] );
  }

  // Likelihood
  Value ~ normal( mu , sigma );
}
"

O2_L_rh_mod <- O2_L_rh_stan %>%
  write_stan_file() %>%
  cmdstan_model()

O2_L_rh_samples <- O2_L %>%
  map(~ O2_L_rh_mod$sample(
    data = .x %>%
      select(Value, delta_t) %>%
      compose_data(),
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4,
    iter_sampling = 1e4))

# 2.6.3 Model checks ####
# check Rhat, effective sample size and chains
O2_L_rh_samples %>%
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

ggsave(
  O2_L_rh_samples %>%
    map(~ .x$draws(format = "df") %>%
          mcmc_rank_overlay()) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_L_rh_rank.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

# 2.6.4 Prior-posterior comparison ####
# sample prior
O2_L_rh_prior <- O2_L %>%
  map(~ prior_samples(model = O2_L_rh_mod,
                      data = .x %>%
                        select(Value, delta_t) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
ggsave(
  O2_L_rh_prior %>%
    map2(O2_L_rh_samples,
         ~ prior_posterior_draws(prior_samples = .x,
                                 posterior_samples = .y,
                                 group = list(NA), # no groups so this has to be an empty list or tibble
                                 parameters = c("O2_max", "beta", "O2_0", "sigma"),
                                 format = "long")) %>%
    imap(~ .x %>% prior_posterior_plot() +
                    ggtitle(.y)) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_L_rh_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# Some posteriors for O2_0 are pretty unconstrained and O2_max sometimes broke out
# of the expected prior probability space, but generally priors are not restrictive.

# 2.6.5 Predictions ####
O2_L_rh_predictions <- O2_L_rh_prior %>%
  map2(O2_L_rh_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("O2_max", "beta", "O2_0", "sigma"),
                               format = "short")) %>%
  map2(O2_L, ~ spread_continuous(.x, .y,
                                 predictor_name = "delta_t",
                                 length = 50)) %>%
  map(~ .x %>% mutate(mu = O2_0 + O2_max * beta * delta_t / ( O2_max + beta * delta_t ),
                      obs = rnorm( n(), mu, sigma )))

O2_L_rh_predictions_summary <- O2_L_rh_predictions %>%
  map(~ .x %>% group_by(distribution, delta_t) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_"))

rm(O2_L_rh_predictions) # remove raw predictions because they are too big

O2_L_rh_predictions_summary %>%
  map2(O2_L,
       ~ ggplot() +
            geom_point(data = .y, aes(delta_t, Value),
                       shape = 16, alpha = 0.05) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(delta_t, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            # geom_ribbon(data = .x %>% filter(distribution == "posterior"), # unhash to check
            #             aes(delta_t, ymin = obs_ymin, ymax = obs_ymax, # predicted observations
            #                 alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_L_rh_prediction.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# The rectangular hyperbola produces an even worse fit.

# 2.7 Photosynthesis: piecewise linear ####
# 2.7.1 Prior simulation ####
tibble(n = 1:1e3,
       O2_max = rgamma(n = 1e3, shape = O2_max^2 / 100^2, rate = O2_max / 100^2),
       beta = rgamma(n = 1e3, shape = beta^2 / 5^2, rate = beta / 5^2),
       O2_0 = rgamma(n = 1e3, shape = O2_0^2 / 50^2, rate = O2_0 / 50^2)) %>%
  expand_grid(t = seq(0, 120)) %>%
  mutate(O2 = if_else(t <= O2_max / beta, O2_0 + beta * t, O2_0 + O2_max)) %>%
  ggplot(aes(t, O2, group = n)) +
    geom_hline(yintercept = c(0, 1e3)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.7.2 Run model ####
O2_L_pl_stan <- "
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
  O2_0 ~ gamma( 230^2 / 50^2 , 230 / 50^2 ); // reparameterised with mean and sd
  beta ~ gamma( 8^2 / 5^2 , 8 / 5^2 );
  O2_max ~ gamma( 600^2 / 100^2 , 600 / 100^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    if ( delta_t[i] <= O2_max / beta ) {
      mu[i] = O2_0 + beta * delta_t[i];
    } else {
      mu[i] = O2_0 + O2_max;
    }
  }

  // Likelihood
  Value ~ normal( mu , sigma );
}
"

O2_L_pl_mod <- O2_L_pl_stan %>%
  write_stan_file() %>%
  cmdstan_model()

O2_L_pl_samples <- O2_L %>%
  map(~ O2_L_pl_mod$sample(
    data = .x %>%
      select(Value, delta_t) %>%
      compose_data(),
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4,
    iter_sampling = 1e4))

# 2.7.3 Model checks ####
# check Rhat, effective sample size and chains
O2_L_pl_samples %>%
  map(~ .x$summary()) %>%
  bind_rows() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# some rhat above 1.001, but usually less than 5%
# decent effective sample size

ggsave(
  O2_L_pl_samples %>%
    map(~ .x$draws(format = "df") %>%
          mcmc_rank_overlay()) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_L_pl_rank.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look mostly fine, but clearly not as good as for the other models,
# with a few seriously deviating chains

# 2.7.4 Prior-posterior comparison ####
# sample prior
O2_L_pl_prior <- O2_L %>%
  map(~ prior_samples(model = O2_L_pl_mod,
                      data = .x %>%
                        select(Value, delta_t) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
ggsave(
  O2_L_pl_prior %>%
    map2(O2_L_pl_samples,
         ~ prior_posterior_draws(prior_samples = .x,
                                 posterior_samples = .y,
                                 group = list(NA), # no groups so this has to be an empty list or tibble
                                 parameters = c("O2_max", "beta", "O2_0", "sigma"),
                                 format = "long")) %>%
    imap(~ .x %>% prior_posterior_plot() +
                    ggtitle(.y)) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_L_pl_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# Some posteriors for O2_0 are pretty unconstrained and O2_max posteriors are frequently
# identical with their prior, meaning there is no information on the plateau in the data.
# Also, some posteriors are not uniquely identified, and have two peaks, e.g. M_11_L.

# 2.7.5 Predictions ####
O2_L_pl_predictions <- O2_L_pl_prior %>%
  map2(O2_L_pl_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("O2_max", "beta", "O2_0", "sigma"),
                               format = "short")) %>%
  map2(O2_L, ~ spread_continuous(.x, .y,
                                 predictor_name = "delta_t",
                                 length = 50)) %>%
  map(~ .x %>% mutate(mu = if_else(delta_t <= O2_max / beta, 
                                   O2_0 + beta * delta_t,
                                   O2_0 + O2_max),
                      obs = rnorm( n(), mu, sigma )))

O2_L_pl_predictions_summary <- O2_L_pl_predictions %>%
  map(~ .x %>% group_by(distribution, delta_t) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_"))

rm(O2_L_pl_predictions) # remove raw predictions because they are too big

O2_L_pl_predictions_summary %>%
  map2(O2_L,
       ~ ggplot() +
            geom_point(data = .y, aes(delta_t, Value),
                       shape = 16, alpha = 0.05) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(delta_t, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            # geom_ribbon(data = .x %>% filter(distribution == "posterior"), # unhash to check
            #             aes(delta_t, ymin = obs_ymin, ymax = obs_ymax, # predicted observations
            #                 alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax),
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_L_pl_prediction.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# The piecewise linear sometimes fits well at the lower end of the predictor but clearly struggles
# to find the transition point and therefore often gets the intercept and slope obviously wrong.

# 2.7 Photosynthesis: comparison of fit ####
O2_L_predictions_summary <-
  O2_L_ht_predictions_summary %>%
  map2(O2_L_es_predictions_summary, ~ bind_rows(.x, .y)) %>%
  map2(O2_L_rh_predictions_summary, ~ bind_rows(.x, .y)) %>%
  map2(O2_L_pl_predictions_summary, ~ bind_rows(.x, .y) %>%
         mutate(Function = c(rep("Hyperbolic tangent", nrow(.y)),
                             rep("Exponential saturation", nrow(.y)),
                             rep("Rectangular hyperbola", nrow(.y)),
                             rep("Piecewise linear", nrow(.y))) %>%
                  fct())
  )

ggsave(
  O2_L_predictions_summary %>%
    map2(O2_L,
         ~ ggplot() +
              geom_point(data = .y, aes(delta_t, Value),
                         shape = 16, alpha = 0.05) +
              geom_line(data = .x %>% filter(distribution == "posterior"),
                        aes(delta_t, mu_y, colour = Function)) +
              geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                          aes(delta_t, ymin = mu_ymin, ymax = mu_ymax,
                              alpha = factor(mu_.width), fill = Function)) +
              scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
              theme_minimal() +
              theme(panel.grid = element_blank())
         ) %>%
    imap(~ .x + ggtitle(.y)) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_L_prediction.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# Hyperbolic tangent is the best followed by the exponential saturation function.
# The piecewise linear and rectangular hyperbola are worst.

# Pull out a few example plots for the supplement (M_4_L, J_5_L, A_6_L)
O2_L_predictions_summary_selection <- O2_L_predictions_summary %$%
  bind_rows(X240620_M_4_L, X240708_J_5_L, X240723_A_6_L) %>%
  mutate(Species = c("Metamastophora flabellata",
                     "Jania rosea",
                     "Amphiroa anceps") %>% 
           rep(each = O2_L_predictions_summary$X240620_M_4_L %>% nrow()) %>% 
           fct_relevel("Amphiroa anceps", "Jania rosea"))

O2_L_selection <- O2_L %$%
  bind_rows(X240620_M_4_L, X240708_J_5_L, X240723_A_6_L) %>%
  mutate(Species = case_when(
    ID %>% str_detect("M") ~ "Metamastophora flabellata",
    ID %>% str_detect("J") ~ "Jania rosea",
    ID %>% str_detect("A") ~ "Amphiroa anceps"
  ) %>% fct_relevel("Amphiroa anceps", "Jania rosea"))

# Define custom theme
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = "cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 12, hjust = 0),
                 axis.text = element_text(size = 10, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black", lineend = "square"),
                 legend.key = element_blank(),
                 legend.key.width = unit(.25, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.key.spacing.x = unit(.5, "cm"),
                 legend.key.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.position = "top",
                 legend.justification = 0,
                 legend.text = element_text(size = 12, hjust = 0),
                 legend.title = element_blank(),
                 legend.margin = margin(0, 0, 0, 0, unit = "cm"),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 12, hjust = 0),
                 panel.spacing = unit(0.6, "cm"),
                 text = element_text(family = "Futura"))

Fig_S1a <- 
  O2_L_predictions_summary_selection %>%
    ggplot() +
      geom_point(data = O2_L_selection, 
                 aes(delta_t, Value),
                 shape = 16, alpha = 0.05) +
      geom_line(data = . %>% filter(distribution == "posterior"),
                aes(delta_t, mu_y, colour = Function),
                linewidth = 0.2) +
      geom_ribbon(data = . %>% filter(distribution == "posterior"),
                  aes(delta_t, ymin = mu_ymin, ymax = mu_ymax,
                  alpha = factor(mu_.width), fill = Function)) +
      scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
      scale_colour_manual(values = c("#f5a54a", "#bc90c1", "#5bb5b5", "#2e4a5b")) +
      scale_fill_manual(values = c("#f5a54a", "#bc90c1", "#5bb5b5", "#2e4a5b")) +
      facet_grid(~ Species) +
      labs(x = "Incubation time (min)", y = expression("O"[2]*" (µM)")) +
      coord_cartesian(xlim = c(0, 120), ylim = c(0, 1200),
                      clip = "off", expand = FALSE) +
      mytheme +
      theme(strip.text = element_text(face = "italic"))

Fig_S1a %>%
  ggsave(filename = "Fig_S1a.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# Retain only optimal estimates.
rm(list = setdiff(ls(), c("mytheme", "O2", "O2_B_D", "O2_B_L", "O2_D",
                          "O2_L_ht_samples", "beta_r", "O2_0", 
                          "prior_posterior_draws", "prior_posterior_plot",
                          "prior_samples", "spread_continuous")))

# 2.8 Respiration: hyperbolic tangent ####
# 2.8.1 Prior simulation ####
tibble(n = 1:1e3,
       O2_0 = rgamma(n = 1e3, shape = O2_0^2 / 50^2, rate = O2_0 / 50^2),
       beta = rgamma(n = 1e3, shape = beta_r^2 / 1^2, rate = beta_r / 1^2)) %>%
  expand_grid(t = seq(0, 120)) %>%
  mutate(O2 = O2_0 * ( 1 - tanh( beta * t / O2_0 ) )) %>%
  ggplot(aes(t, O2, group = n)) +
    geom_hline(yintercept = 0) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# 2.8.2 Run model ####
O2_D_ht_stan <- "
data{
  int n;
  vector<lower=0>[n] Value;
  vector<lower=0>[n] delta_t;
}

parameters{
  real<lower=0> O2_0;
  real<lower=0> beta;
  real<lower=0> sigma;
}

model{
  // Priors
  O2_0 ~ gamma( 230^2 / 50^2 , 230 / 50^2 ); // reparameterised with mean and sd
  beta ~ gamma( 2^2 / 1^2 , 2 / 1^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = O2_0 * ( 1 - tanh( beta * delta_t[i] / O2_0 ) );
  }

  // Likelihood
  Value ~ normal( mu , sigma );
}
"

O2_D_ht_mod <- O2_D_ht_stan %>%
  write_stan_file() %>%
  cmdstan_model()

O2_D_ht_samples <- O2_D %>%
  map(~ O2_D_ht_mod$sample(
    data = .x %>%
      select(Value, delta_t) %>%
      compose_data(),
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4,
    iter_sampling = 1e4))

# 2.8.3 Model checks ####
# check Rhat, effective sample size and chains
O2_D_ht_samples %>%
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

ggsave(
  O2_D_ht_samples %>%
    map(~ .x$draws(format = "df") %>%
          mcmc_rank_overlay()) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_D_ht_rank.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look fine, but there is some structure in some
# of them, potentially due to wiggliness in the data

# 2.8.4 Prior-posterior comparison ####
# sample prior
O2_D_ht_prior <- O2_D %>%
  map(~ prior_samples(model = O2_D_ht_mod,
                      data = .x %>%
                        select(Value, delta_t) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
ggsave(
  O2_D_ht_prior %>%
    map2(O2_D_ht_samples,
         ~ prior_posterior_draws(prior_samples = .x,
                                 posterior_samples = .y,
                                 group = list(NA), # no groups so this has to be an empty list or tibble
                                 parameters = c("O2_0", "beta", "sigma"),
                                 format = "long")) %>%
    imap(~ .x %>% prior_posterior_plot() +
                    ggtitle(.y)) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_D_ht_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# posteriors are very well constrained

# 2.8.5 Predictions ####
O2_D_ht_predictions <- O2_D_ht_prior %>%
  map2(O2_D_ht_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("O2_0", "beta", "sigma"),
                               format = "short")) %>%
  map2(O2_D, ~ spread_continuous(.x, .y,
                                 predictor_name = "delta_t",
                                 length = 50)) %>%
  map(~ .x %>% mutate(mu = O2_0 * ( 1 - tanh( beta * delta_t / O2_0 ) ),
                      obs = rnorm( n(), mu, sigma )))

O2_D_ht_predictions_summary <- O2_D_ht_predictions %>%
  map(~ .x %>% group_by(distribution, delta_t) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_"))

rm(O2_D_ht_predictions) # remove raw predictions because they are too big

O2_D_ht_predictions_summary %>%
  map2(O2_D,
       ~ ggplot() +
            geom_point(data = .y, aes(delta_t, Value),
                       shape = 16, alpha = 0.05) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(delta_t, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            # geom_ribbon(data = .x %>% filter(distribution == "posterior"), # unhash to check
            #             aes(delta_t, ymin = obs_ymin, ymax = obs_ymax, # predicted observations
            #                 alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax), # unhash to check prior
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_D_ht_prediction.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# The hyperbolic tangent struggles to describe the approach of zero
# in most cases, usually by bending off too early.

# 2.9 Respiration: exponential decay ####
# 2.9.1 Prior simulation ####
tibble(n = 1:1e3,
       O2_0 = rgamma(n = 1e3, shape = O2_0^2 / 50^2, rate = O2_0 / 50^2),
       beta = rgamma(n = 1e3, shape = beta_r^2 / 1^2, rate = beta_r / 1^2)) %>%
  expand_grid(t = seq(0, 120)) %>%
  mutate(O2 = O2_0 * exp( -beta * t / O2_0 )) %>%
  ggplot(aes(t, O2, group = n)) +
    geom_hline(yintercept = 0) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# 2.9.2 Run model ####
O2_D_ed_stan <- "
data{
  int n;
  vector<lower=0>[n] Value;
  vector<lower=0>[n] delta_t;
}

parameters{
  real<lower=0> O2_0;
  real<lower=0> beta;
  real<lower=0> sigma;
}

model{
  // Priors
  O2_0 ~ gamma( 230^2 / 50^2 , 230 / 50^2 ); // reparameterised with mean and sd
  beta ~ gamma( 2^2 / 1^2 , 2 / 1^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = O2_0 * exp( -beta * delta_t[i] / O2_0 );
  }

  // Likelihood
  Value ~ normal( mu , sigma );
}
"

O2_D_ed_mod <- O2_D_ed_stan %>%
  write_stan_file() %>%
  cmdstan_model()

O2_D_ed_samples <- O2_D %>%
  map(~ O2_D_ed_mod$sample(
    data = .x %>%
      select(Value, delta_t) %>%
      compose_data(),
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4,
    iter_sampling = 1e4))

# 2.9.3 Model checks ####
# check Rhat, effective sample size and chains
O2_D_ed_samples %>%
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

ggsave(
  O2_D_ed_samples %>%
    map(~ .x$draws(format = "df") %>%
          mcmc_rank_overlay()) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_D_ed_rank.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look fine, but the same wiggliness is present

# 2.9.4 Prior-posterior comparison ####
# sample prior
O2_D_ed_prior <- O2_D %>%
  map(~ prior_samples(model = O2_D_ed_mod,
                      data = .x %>%
                        select(Value, delta_t) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
ggsave(
  O2_D_ed_prior %>%
    map2(O2_D_ed_samples,
         ~ prior_posterior_draws(prior_samples = .x,
                                 posterior_samples = .y,
                                 group = list(NA), # no groups so this has to be an empty list or tibble
                                 parameters = c("O2_0", "beta", "sigma"),
                                 format = "long")) %>%
    imap(~ .x %>% prior_posterior_plot() +
                    ggtitle(.y)) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_D_ed_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# posteriors are very well constrained

# 2.9.5 Predictions ####
O2_D_ed_predictions <- O2_D_ed_prior %>%
  map2(O2_D_ed_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("O2_0", "beta", "sigma"),
                               format = "short")) %>%
  map2(O2_D, ~ spread_continuous(.x, .y,
                                 predictor_name = "delta_t",
                                 length = 50)) %>%
  map(~ .x %>% mutate(mu = O2_0 * exp( -beta * delta_t / O2_0 ),
                      obs = rnorm( n(), mu, sigma )))

O2_D_ed_predictions_summary <- O2_D_ed_predictions %>%
  map(~ .x %>% group_by(distribution, delta_t) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_"))

rm(O2_D_ed_predictions) # remove raw predictions because they are too big

O2_D_ed_predictions_summary %>%
  map2(O2_D,
       ~ ggplot() +
            geom_point(data = .y, aes(delta_t, Value),
                       shape = 16, alpha = 0.05) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(delta_t, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            # geom_ribbon(data = .x %>% filter(distribution == "posterior"), # unhash to check
            #             aes(delta_t, ymin = obs_ymin, ymax = obs_ymax, # predicted observations
            #                 alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax), # unhash to check prior
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_D_ed_prediction.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# The exponential decay function struggles even more.

# 2.10 Respiration: rectangular hyperbola ####
# 2.10.1 Prior simulation ####
tibble(n = 1:1e3,
       O2_0 = rgamma(n = 1e3, shape = O2_0^2 / 50^2, rate = O2_0 / 50^2),
       beta = rgamma(n = 1e3, shape = beta_r^2 / 1^2, rate = beta_r / 1^2)) %>%
  expand_grid(t = seq(0, 120)) %>%
  mutate(O2 = O2_0^2 / ( O2_0 + beta * t )) %>%
  ggplot(aes(t, O2, group = n)) +
    geom_hline(yintercept = 0) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable, although predictions approach zero more slowly.

# 2.10.2 Run model ####
O2_D_rh_stan <- "
data{
  int n;
  vector<lower=0>[n] Value;
  vector<lower=0>[n] delta_t;
}

parameters{
  real<lower=0> O2_0;
  real<lower=0> beta;
  real<lower=0> sigma;
}

model{
  // Priors
  O2_0 ~ gamma( 230^2 / 50^2 , 230 / 50^2 ); // reparameterised with mean and sd
  beta ~ gamma( 2^2 / 1^2 , 2 / 1^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = O2_0^2 / ( O2_0 + beta * delta_t[i] );
  }

  // Likelihood
  Value ~ normal( mu , sigma );
}
"

O2_D_rh_mod <- O2_D_rh_stan %>%
  write_stan_file() %>%
  cmdstan_model()

O2_D_rh_samples <- O2_D %>%
  map(~ O2_D_rh_mod$sample(
    data = .x %>%
      select(Value, delta_t) %>%
      compose_data(),
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4,
    iter_sampling = 1e4))

# 2.10.3 Model checks ####
# check Rhat, effective sample size and chains
O2_D_rh_samples %>%
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

ggsave(
  O2_D_rh_samples %>%
    map(~ .x$draws(format = "df") %>%
          mcmc_rank_overlay()) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_D_rh_rank.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look fine, but same wiggliness

# 2.10.4 Prior-posterior comparison ####
# sample prior
O2_D_rh_prior <- O2_D %>%
  map(~ prior_samples(model = O2_D_rh_mod,
                      data = .x %>%
                        select(Value, delta_t) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
ggsave(
  O2_D_rh_prior %>%
    map2(O2_D_rh_samples,
         ~ prior_posterior_draws(prior_samples = .x,
                                 posterior_samples = .y,
                                 group = list(NA), # no groups so this has to be an empty list or tibble
                                 parameters = c("O2_0", "beta", "sigma"),
                                 format = "long")) %>%
    imap(~ .x %>% prior_posterior_plot() +
                    ggtitle(.y)) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_D_rh_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# posteriors are very well constrained

# 2.10.5 Predictions ####
O2_D_rh_predictions <- O2_D_rh_prior %>%
  map2(O2_D_rh_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("O2_0", "beta", "sigma"),
                               format = "short")) %>%
  map2(O2_D, ~ spread_continuous(.x, .y,
                                 predictor_name = "delta_t",
                                 length = 50)) %>%
  map(~ .x %>% mutate(mu = O2_0^2 / ( O2_0 + beta * delta_t ),
                      obs = rnorm( n(), mu, sigma )))

O2_D_rh_predictions_summary <- O2_D_rh_predictions %>%
  map(~ .x %>% group_by(distribution, delta_t) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_"))

rm(O2_D_rh_predictions) # remove raw predictions because they are too big

O2_D_rh_predictions_summary %>%
  map2(O2_D,
       ~ ggplot() +
            geom_point(data = .y, aes(delta_t, Value),
                       shape = 16, alpha = 0.05) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(delta_t, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            # geom_ribbon(data = .x %>% filter(distribution == "posterior"), # unhash to check
            #             aes(delta_t, ymin = obs_ymin, ymax = obs_ymax, # predicted observations
            #                 alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax), # unhash to check prior
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_D_rh_prediction.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# The rectangular hyperbola fits worst so far.

# 2.11 Respiration: piecewise linear ####
# 2.11.1 Prior simulation ####
tibble(n = 1:1e3,
       O2_0 = rgamma(n = 1e3, shape = O2_0^2 / 50^2, rate = O2_0 / 50^2),
       beta = rgamma(n = 1e3, shape = beta_r^2 / 1^2, rate = beta_r / 1^2)) %>%
  expand_grid(t = seq(0, 120)) %>%
  mutate(O2 = if_else(t <= O2_0 / beta, O2_0 - beta * t, 0)) %>%
  ggplot(aes(t, O2, group = n)) +
    geom_hline(yintercept = 0) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# 2.11.2 Run model ####
O2_D_pl_stan <- "
data{
  int n;
  vector<lower=0>[n] Value;
  vector<lower=0>[n] delta_t;
}

parameters{
  real<lower=0> O2_0;
  real<lower=0> beta;
  real<lower=0> sigma;
}

model{
  // Priors
  O2_0 ~ gamma( 230^2 / 50^2 , 230 / 50^2 ); // reparameterised with mean and sd
  beta ~ gamma( 2^2 / 1^2 , 2 / 1^2 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    if ( delta_t[i] <= O2_0 / beta ) {
      mu[i] = O2_0 - beta * delta_t[i];
    } else {
      mu[i] = 0;
    }
  }

  // Likelihood
  Value ~ normal( mu , sigma );
}
"

O2_D_pl_mod <- O2_D_pl_stan %>%
  write_stan_file() %>%
  cmdstan_model()

O2_D_pl_samples <- O2_D %>%
  map(~ O2_D_pl_mod$sample(
    data = .x %>%
      select(Value, delta_t) %>%
      compose_data(),
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4,
    iter_sampling = 1e4))

# 2.11.3 Model checks ####
# check Rhat, effective sample size and chains
O2_D_pl_samples %>%
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

ggsave(
  O2_D_pl_samples %>%
    map(~ .x$draws(format = "df") %>%
          mcmc_rank_overlay()) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_D_pl_rank.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look fine, same wiggliness

# 2.11.4 Prior-posterior comparison ####
# sample prior
O2_D_pl_prior <- O2_D %>%
  map(~ prior_samples(model = O2_D_pl_mod,
                      data = .x %>%
                        select(Value, delta_t) %>%
                        compose_data(),
                      chains = 8, samples = 1e4))

# plot prior-posterior comparison
ggsave(
  O2_D_pl_prior %>%
    map2(O2_D_pl_samples,
         ~ prior_posterior_draws(prior_samples = .x,
                                 posterior_samples = .y,
                                 group = list(NA), # no groups so this has to be an empty list or tibble
                                 parameters = c("O2_0", "beta", "sigma"),
                                 format = "long")) %>%
    imap(~ .x %>% prior_posterior_plot() +
                    ggtitle(.y)) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_D_pl_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# posteriors are very well constrained

# 2.11.5 Predictions ####
O2_D_pl_predictions <- O2_D_pl_prior %>%
  map2(O2_D_pl_samples,
       ~ prior_posterior_draws(prior_samples = .x,
                               posterior_samples = .y,
                               group = list(NA),
                               parameters = c("O2_0", "beta", "sigma"),
                               format = "short")) %>%
  map2(O2_D, ~ spread_continuous(.x, .y,
                                 predictor_name = "delta_t",
                                 length = 50)) %>%
  map(~ .x %>% mutate(mu = if_else(delta_t <= O2_0 / beta, 
                                   O2_0 - beta * delta_t, 
                                   0),
                      obs = rnorm( n(), mu, sigma )))

O2_D_pl_predictions_summary <- O2_D_pl_predictions %>%
  map(~ .x %>% group_by(distribution, delta_t) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_"))

rm(O2_D_pl_predictions) # remove raw predictions because they are too big

O2_D_pl_predictions_summary %>%
  map2(O2_D,
       ~ ggplot() +
            geom_point(data = .y, aes(delta_t, Value),
                       shape = 16, alpha = 0.05) +
            geom_line(data = .x %>% filter(distribution == "posterior"),
                      aes(delta_t, mu_y)) +
            geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax,
                            alpha = factor(mu_.width))) +
            # geom_ribbon(data = .x %>% filter(distribution == "posterior"), # unhash to check
            #             aes(delta_t, ymin = obs_ymin, ymax = obs_ymax, # predicted observations
            #                 alpha = factor(obs_.width))) +
            geom_ribbon(data = .x %>% filter(distribution == "prior", mu_.width == 0.9),
                        aes(delta_t, ymin = mu_ymin, ymax = mu_ymax), # unhash to check prior
                        colour = alpha("black", 0.3), fill = NA) +
            scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
            theme_minimal() +
            theme(panel.grid = element_blank())
       ) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_D_pl_prediction.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# The piecewise linear function seems to fit best in most cases.

# 2.12 Respiration: comparison of fit ####
O2_D_predictions_summary <-
  O2_D_ht_predictions_summary %>%
  map2(O2_D_ed_predictions_summary, ~ bind_rows(.x, .y)) %>%
  map2(O2_D_rh_predictions_summary, ~ bind_rows(.x, .y)) %>%
  map2(O2_D_pl_predictions_summary, ~ bind_rows(.x, .y) %>%
         mutate(Function = c(rep("Hyperbolic tangent", nrow(.y)),
                             rep("Exponential decay", nrow(.y)),
                             rep("Rectangular hyperbola", nrow(.y)),
                             rep("Piecewise linear", nrow(.y))) %>%
                  fct())
  )

ggsave(
  O2_D_predictions_summary %>%
    map2(O2_D,
         ~ ggplot() +
              geom_point(data = .y, aes(delta_t, Value),
                         shape = 16, alpha = 0.05) +
              geom_line(data = .x %>% filter(distribution == "posterior"),
                        aes(delta_t, mu_y, colour = Function)) +
              geom_ribbon(data = .x %>% filter(distribution == "posterior"),
                          aes(delta_t, ymin = mu_ymin, ymax = mu_ymax,
                              alpha = factor(mu_.width), fill = Function)) +
              scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
              theme_minimal() +
              theme(panel.grid = element_blank())
         ) %>%
    imap(~ .x + ggtitle(.y)) %>%
    wrap_plots() +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  filename = "O2_D_prediction.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)
# Piecewise linear is the best followed by hyperbolic tangent. Exponential decay and
# rectangular hyperbola are clearly worse.

# Pull out a few example plots for the supplement (M_4_L, J_5_L, A_6_L)
O2_D_predictions_summary_selection <- O2_D_predictions_summary %$%
  bind_rows(X240521_M_1_D, X240708_J_5_D, X240620_A_4_D) %>%
  mutate(Species = c("Metamastophora flabellata",
                     "Jania rosea",
                     "Amphiroa anceps") %>% 
           rep(each = O2_D_predictions_summary$X240521_M_1_D %>% nrow()) %>% 
           fct_relevel("Amphiroa anceps", "Jania rosea"))

O2_D_selection <- O2_D %$%
  bind_rows(X240521_M_1_D, X240708_J_5_D, X240620_A_4_D) %>%
  mutate(Species = case_when(
    ID %>% str_detect("M") ~ "Metamastophora flabellata",
    ID %>% str_detect("J") ~ "Jania rosea",
    ID %>% str_detect("A") ~ "Amphiroa anceps"
  ) %>% fct_relevel("Amphiroa anceps", "Jania rosea"))

Fig_S1b <- 
  O2_D_predictions_summary_selection %>%
    ggplot() +
      geom_point(data = O2_D_selection, 
                 aes(delta_t, Value),
                 shape = 16, alpha = 0.05) +
      geom_line(data = . %>% filter(distribution == "posterior"),
                aes(delta_t, mu_y, colour = Function),
                linewidth = 0.2) +
      geom_ribbon(data = . %>% filter(distribution == "posterior"), 
                  aes(delta_t, ymin = mu_ymin, ymax = mu_ymax,
                  alpha = factor(mu_.width), fill = Function)) +
      scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
      scale_colour_manual(values = c("#f5a54a", "#bc90c1", "#5bb5b5", "#2e4a5b")) +
      scale_fill_manual(values = c("#f5a54a", "#bc90c1", "#5bb5b5", "#2e4a5b")) +
      facet_grid(~ Species) +
      labs(x = "Incubation time (min)", y = expression("O"[2]*" (µM)")) +
      coord_cartesian(xlim = c(0, 120), ylim = c(0, 300),
                      clip = "off", expand = FALSE) +
      mytheme +
      theme(strip.text = element_text(face = "italic"))

Fig_S1b %>%
  ggsave(filename = "Fig_S1b.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# Retain only optimal estimates.
rm(list = setdiff(ls(), c("mytheme", "O2", "O2_B_D", "O2_B_L",
                          "O2_L_ht_samples", "O2_D_pl_samples", 
                          "prior_posterior_draws", "prior_posterior_plot",
                          "prior_samples", "spread_continuous")))

# 2.13 Light blanks: multilevel linear ####
# 2.13.1 Prior simulation ####
tibble(n = 1:1e3,
       alpha = rgamma(n = 1e3, shape = 230^2 / 50^2, rate = 230 / 50^2), # same as O2_0
       beta = rnorm(n = 1e3, mean = 0, sd = 1)) %>% # blank O2 can slightly increase or decrease
  expand_grid(t_c = seq(-60, 60)) %>% # use centred incubation time to decorrelate linear parameters
  mutate(O2 = alpha + beta * t_c) %>%
  ggplot(aes(t_c, O2, group = n)) +
    geom_hline(yintercept = 0) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# More than enough variation for filtered and UV-sterilised seawater blanks.

# 2.13.2 Run model ####
O2_B_L_stan <- "
data{
  int n;
  vector<lower=0>[n] Value;
  vector[n] delta_t_c; // centred time can take negative values
  array[n] int ID;
  int n_ID;
}

parameters{
  // Hyperparameters
  real beta_mu;
  real<lower=0> beta_sigma;
  
  // ID parameters
  vector<lower=0>[n_ID] alpha;
  vector[n_ID] beta;
  
  // Likelihood uncertainty parameter
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  beta_mu ~ normal( 0 , 1 );
  beta_sigma ~ exponential( 1 );

  // ID priors
  alpha ~ gamma( 230^2 / 50^2 , 230 / 50^2 ); // reparameterised with mean and sd
  beta ~ normal( beta_mu , beta_sigma );
  
  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );
  
  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = alpha[ID[i]] + beta[ID[i]] * delta_t_c[i];
  }

  // Likelihood
  Value ~ normal( mu , sigma );
}
"

O2_B_L_mod <- O2_B_L_stan %>%
  write_stan_file() %>%
  cmdstan_model()

O2_B_L_samples <- O2_B_L_mod$sample(
    data = O2_B_L %>%
      bind_rows() %>%
      select(Value, delta_t_c, ID) %>%
      compose_data(),
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4,
    iter_sampling = 1e4)

# 2.13.3 Model checks ####
# check Rhat, effective sample size and chains
O2_B_L_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

O2_B_L_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# chains look good

O2_B_L_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[1]", "beta[1]"))
O2_B_L_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[16]", "beta[16]"))
# no correlation between slope and intercept

# 2.13.4 Prior-posterior comparison ####
# sample prior
O2_B_L_prior <- prior_samples(
  model = O2_B_L_mod,
  data = O2_B_L %>%
    bind_rows() %>%
    select(Value, delta_t_c, ID) %>%
    compose_data(),
  chains = 8, samples = 1e4)
# The function prior_samples does not fully support multilevel models and produces
# somewhat noisy prior distributions, but they suffice for diagnostic purposes.

# plot prior-posterior comparison
O2_B_L_prior %>%
  prior_posterior_draws(posterior_samples = O2_B_L_samples,
                        group = O2_B_L %>%
                          bind_rows() %>%
                          select(ID),
                        parameters = c("beta_mu", "beta_sigma", 
                                       "alpha[ID]", "beta[ID]", "sigma"),
                        format = "long") %>%
  prior_posterior_plot(group_name = "ID", ridges = FALSE)
# posteriors are very well constrained

# 2.13.5 Predictions ####
O2_B_L_predictions <- O2_B_L_prior %>%
  prior_posterior_draws(posterior_samples = O2_B_L_samples,
                        group = O2_B_L %>%
                          bind_rows() %>%
                          select(ID),
                        parameters = c("alpha[ID]", "beta[ID]", "sigma"),
                        format = "short") %>%
  spread_continuous(O2_B_L %>%
                      bind_rows(),
                    predictor_name = "delta_t_c",
                    group_name = "ID",
                    length = 50) %>%
  mutate(mu = alpha + beta * delta_t_c,
         obs = rnorm( n(), mu, sigma ))

O2_B_L_predictions_summary <- O2_B_L_predictions %>%
  group_by(distribution, delta_t_c, ID) %>%
  reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
          obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
  unnest(c(mu, obs), names_sep = "_")

rm(O2_B_L_predictions) # remove raw predictions because they are too big

O2_B_L_predictions_summary %>%
  ggplot() +
    geom_point(data = O2_B_L %>%
                 bind_rows(), 
               aes(delta_t_c, Value),
               shape = 16, alpha = 0.05) +
    geom_line(data = . %>% filter(distribution == "posterior"),
              aes(delta_t_c, mu_y)) +
    geom_ribbon(data = . %>% filter(distribution == "posterior"),
                aes(delta_t_c, ymin = mu_ymin, ymax = mu_ymax,
                    alpha = factor(mu_.width))) +
    # geom_ribbon(data = . %>% filter(distribution == "posterior"), # unhash to check
    #             aes(delta_t_c, ymin = obs_ymin, ymax = obs_ymax, # predicted observations
    #                 alpha = factor(obs_.width))) +
    geom_ribbon(data = . %>% filter(distribution == "prior", mu_.width == 0.9),
                aes(delta_t_c, ymin = mu_ymin, ymax = mu_ymax), # unhash to check prior
                    colour = alpha("black", 0.3), fill = NA) +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    facet_wrap(~ ID) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Linear model fits well in all cases.

# 2.14 Dark blanks: multilevel linear ####
# 2.14.1 Prior simulation ####

# See above.

# 2.14.2 Run model ####
O2_B_D_stan <- "
data{
  int n;
  vector<lower=0>[n] Value;
  vector[n] delta_t_c; // centred time can take negative values
  array[n] int ID;
  int n_ID;
}

parameters{
  // Hyperparameters
  real beta_mu;
  real<lower=0> beta_sigma;
  
  // ID parameters
  vector<lower=0>[n_ID] alpha;
  vector[n_ID] beta;
  
  // Likelihood uncertainty parameter
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  beta_mu ~ normal( 0 , 1 );
  beta_sigma ~ exponential( 1 );

  // ID priors
  alpha ~ gamma( 230^2 / 50^2 , 230 / 50^2 ); // reparameterised with mean and sd
  beta ~ normal( beta_mu , beta_sigma );
  
  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );
  
  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = alpha[ID[i]] - beta[ID[i]] * delta_t_c[i];
  }

  // Likelihood
  Value ~ normal( mu , sigma );
}
"

O2_B_D_mod <- O2_B_D_stan %>%
  write_stan_file() %>%
  cmdstan_model()

O2_B_D_samples <- O2_B_D_mod$sample(
    data = O2_B_D %>%
      bind_rows() %>%
      select(Value, delta_t_c, ID) %>%
      compose_data(),
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4,
    iter_sampling = 1e4)

# 2.14.3 Model checks ####
# check Rhat, effective sample size and chains
O2_B_D_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

O2_B_D_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# chains look fine but somewhat wiggly, likely due to structure in data

O2_B_D_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[1]", "beta[1]"))
O2_B_D_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[16]", "beta[16]"))
# no correlation between slope and intercept

# 2.14.4 Prior-posterior comparison ####
# sample prior
O2_B_D_prior <- prior_samples(
  model = O2_B_D_mod,
  data = O2_B_D %>%
    bind_rows() %>%
    select(Value, delta_t_c, ID) %>%
    compose_data(),
  chains = 8, samples = 1e4)

# plot prior-posterior comparison
O2_B_D_prior %>%
  prior_posterior_draws(posterior_samples = O2_B_D_samples,
                        group = O2_B_D %>%
                          bind_rows() %>%
                          select(ID),
                        parameters = c("beta_mu", "beta_sigma", 
                                       "alpha[ID]", "beta[ID]", "sigma"),
                        format = "long") %>%
  prior_posterior_plot(group_name = "ID", ridges = FALSE)
# posteriors are very well constrained

# 2.14.5 Predictions ####
O2_B_D_predictions <- O2_B_D_prior %>%
  prior_posterior_draws(posterior_samples = O2_B_D_samples,
                        group = O2_B_D %>%
                          bind_rows() %>%
                          select(ID),
                        parameters = c("alpha[ID]", "beta[ID]", "sigma"),
                        format = "short") %>%
  spread_continuous(O2_B_D %>%
                      bind_rows(),
                    predictor_name = "delta_t_c",
                    group_name = "ID",
                    length = 50) %>%
  mutate(mu = alpha - beta * delta_t_c,
         obs = rnorm( n(), mu, sigma ))

O2_B_D_predictions_summary <- O2_B_D_predictions %>%
  group_by(distribution, delta_t_c, ID) %>%
  reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
          obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
  unnest(c(mu, obs), names_sep = "_")

rm(O2_B_D_predictions) # remove raw predictions because they are too big

O2_B_D_predictions_summary %>%
  ggplot() +
    geom_point(data = O2_B_D %>%
                 bind_rows(), 
               aes(delta_t_c, Value),
               shape = 16, alpha = 0.05) +
    geom_line(data = . %>% filter(distribution == "posterior"),
              aes(delta_t_c, mu_y)) +
    geom_ribbon(data = . %>% filter(distribution == "posterior"),
                aes(delta_t_c, ymin = mu_ymin, ymax = mu_ymax,
                    alpha = factor(mu_.width))) +
    # geom_ribbon(data = . %>% filter(distribution == "posterior"), # unhash to check
    #             aes(delta_t_c, ymin = obs_ymin, ymax = obs_ymax, # predicted observations
    #                 alpha = factor(obs_.width))) +
    geom_ribbon(data = . %>% filter(distribution == "prior", mu_.width == 0.9),
                aes(delta_t_c, ymin = mu_ymin, ymax = mu_ymax), # unhash to check prior
                    colour = alpha("black", 0.3), fill = NA) +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    facet_wrap(~ ID) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Linear model fits well in all cases.

# Retain raw estimates only.
rm(list = setdiff(ls(), c("mytheme", "O2", "O2_B_D", "O2_B_L", 
                          "O2_L_ht_samples", "O2_D_pl_samples",
                          "O2_B_L_samples", "O2_B_D_samples",
                          "prior_posterior_draws", "prior_posterior_plot",
                          "prior_samples", "spread_continuous")))

# 3. Confounders ####
# Potential confounders of O2 measurement include initial O2, temperature, pressure, 
# initial salinity and sample mass. Initial O2 is already estimated via the intercept.

# 3.1 Salinity ####
# There are two initial salinity records, one manually input into the O2 measurement 
# software, and another from the incubation metadata file:
incu_meta <- read.csv("Incubation.csv")

# Check if salinity records are identical
S <- O2 %>%
  bind_rows() %>%
  mutate(Round = if_else(
         str_length(ID) > 6,
         ID %>% 
           str_split(pattern = "_", n = 4) %>% 
           map_chr(~ str_c(.[2], .[4], sep = "_")),
         ID %>% 
           str_split(pattern = "_", n = 2) %>% 
           map_chr(2)
         ) %>% fct()) %>%
  group_by(Round) %>%
  summarise(Salinity = mean(Salinity)) %>%
  ungroup() %>%
  left_join(incu_meta %>% 
              filter(ID %>% str_detect("^0")) %>%
              mutate(Round = if_else(
                str_length(ID) > 6,
                ID %>% 
                  str_split(pattern = "_", n = 4) %>% 
                  map_chr(~ str_c(.[2], .[4], sep = "_")),
                ID %>% 
                  str_split(pattern = "_", n = 2) %>% 
                  map_chr(2) 
              ) %>% fct()) %>%
              select(Round, Salinity),
            by = "Round") %>%
  mutate(Identical = Salinity.x == Salinity.y) 

S %>% print(n = 34)
# Salinities are very similar but not always identical. The values from the 
# incubation metadata are more reliable.
S %<>%
  select(Round, Salinity.y) %>%
  rename(Salinity = Salinity.y)

# 3.2 Pressure ####
O2 %>%
  bind_rows() %>%
  group_by(ID) %>%
  summarise(T_mean = mean(Temp),
            T_sd = sd(Temp),
            P_mean = mean(Pressure),
            P_sd = sd(Pressure),
            n = length(Pressure)) %>%
  print(n = 128)

# For many rounds sensors were not grouped so the correct temperature and pressure are only
# displayed for one of the four channels, the others showing the default 20°C and 1013 hPa
# without variation. These defaults were also used during calibration, so temperature- and
# pressure-correction of O2 should still be fairly accurate. But for the confounders we are 
# interested in potential biological effects of incubation temperature and pressure. We can 
# extract meaningful measurements by filtering only for the master channel, which was B in
# most cases but A in round 14. B is missing for round 7 so any other channel will do.
O2 %>%
  bind_rows() %>%
  filter(ID %>% str_detect("B") & !ID %>% str_detect("B_14") | 
           ID %>% str_detect("A_14") | ID %>% str_detect("A_7")) %>%
  mutate(ID = ID %>% fct_drop()) %>%
  group_by(ID) %>%
  summarise(T_mean = mean(Temp),
            T_sd = sd(Temp),
            P_mean = mean(Pressure),
            P_sd = sd(Pressure),
            n = length(Pressure)) %>%
  print(n = 34)

# Unlike temperature, pressure doesn't have variation in three rounds (9_L, 13_D, 16_D). 
# It therefore cannot be modelled (sd = 0 means no variation to estimate) and I will 
# simply use the mean: 
P <- O2 %>%
  bind_rows() %>%
  filter(ID %>% str_detect("B") & !ID %>% str_detect("B_14") | 
           ID %>% str_detect("A_14") | ID %>% str_detect("A_7")) %>%
  mutate(Round = ID %>% 
           str_split(pattern = "_", n = 2) %>% 
           map_chr(2) %>% fct()
         ) %>%
  group_by(Round) %>%
  summarise(Pressure = mean(Pressure)) %>%
  ungroup()

P %>% print(n = 34)

# 3.3 Temperature ####
# Temperature always has some variation despite the incubations taking place in an 
# 18°C controlled-temperature room, so let's model it.

# 3.3.1 Prepare data ####
T_data <- O2 %>%
  bind_rows() %>%
  filter(ID %>% str_detect("B") & !ID %>% str_detect("B_14") | 
           ID %>% str_detect("A_14") | ID %>% str_detect("A_7")) %>%
  mutate(Round = ID %>% 
           str_split(pattern = "_", n = 2) %>% 
           map_chr(2) %>% fct()
         ) %>%
  select(Round, Temp)

str(T_data)

# 3.3.2 Prior simulation ####

ggplot() +
  geom_density(aes(rgamma(1e5, 18^2 / 2^2, 18 / 2^2))) + # 2°C is a large enough prior sd
  theme_minimal() +
  theme(panel.grid = element_blank()) 

# 3.3.3 Run model ####
T_stan <- "
data{
  int n;
  vector<lower=0>[n] Temp;
  array[n] int Round;
  int n_Round;
}

parameters{
  vector<lower=0>[n_Round] T_mu;
  real<lower=0> sigma;
}

model{
  T_mu ~ gamma( 18^2 / 2^2 , 18 / 2^2 );
  sigma ~ exponential( 1 );
  
  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = T_mu[Round[i]];
  }

  // Likelihood
  Temp ~ normal( mu , sigma );
}
"

T_mod <- T_stan %>%
  write_stan_file() %>%
  cmdstan_model()

T_samples <- T_mod$sample(
    data = T_data %>% compose_data(),
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4,
    iter_sampling = 1e4)

# 3.3.4 Model checks ####
T_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

T_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# chains look good

# 3.3.5 Prior-posterior comparison ####
T_prior <- prior_samples(
  model = T_mod,
  data = T_data %>% compose_data(),
  chains = 8, samples = 1e4)

# plot prior-posterior comparison
T_prior %>%
  prior_posterior_draws(posterior_samples = T_samples,
                        group = T_data %>% select(Round),
                        parameters = c("T_mu[Round]", "sigma"),
                        format = "long") %>%
  prior_posterior_plot(group_name = "Round", ridges = FALSE)
# posteriors are very well constrained

# 4. Conversion ####
# 4.1 Combine estimates ####
# 4.1.1 Combine sample estimates ####
O2_sample_estimates <- O2_L_ht_samples %>%
  imap(~ .x %>% 
         spread_draws(O2_0, beta) %>%
         ungroup() %>%
         mutate(Date = .y %>%
                  str_split(pattern = "_", n = 2) %>% 
                  map_chr(1) %>% str_sub(start = 2) %>%
                  ymd(),
                ID = .y %>% 
                  str_split(pattern = "_", n = 2) %>% 
                  map_chr(2) %>% fct(),
                Round = if_else(
                  str_length(.y) > 14,
                  .y %>% 
                    str_split(pattern = "_", n = 5) %>% 
                    map_chr(~ str_c(.[3], .[5], sep = "_")),
                  .y %>% 
                    str_split(pattern = "_", n = 3) %>% 
                    map_chr(3) 
                  ) %>% fct()
                )
       ) %>%
  map2(O2_D_pl_samples %>%
         imap(~ .x %>% 
                spread_draws(O2_0, beta) %>%
                ungroup() %>%
                mutate(Date = .y %>%
                         str_split(pattern = "_", n = 2) %>% 
                         map_chr(1) %>% str_sub(start = 2) %>%
                         ymd(),
                       ID = .y %>% 
                         str_split(pattern = "_", n = 2) %>% 
                         map_chr(2) %>% fct(),
                       Round = if_else(
                         str_length(.y) > 14,
                         .y %>% 
                           str_split(pattern = "_", n = 5) %>% 
                           map_chr(~ str_c(.[3], .[5], sep = "_")),
                         .y %>% 
                           str_split(pattern = "_", n = 3) %>% 
                           map_chr(3)
                         ) %>% fct()
                       )
              ), ~ .x %>% bind_rows(.y)
       ) %>% bind_rows()

# 4.1.2 Combine blank estimates ####
O2_blank_estimates <- O2_B_L_samples %>%
  recover_types(O2_B_L %>%
                  bind_rows() %>%
                  select(ID)) %>%
  spread_draws(beta[ID]) %>%
  ungroup() %>%
  bind_rows(
    O2_B_L_samples %>%
      recover_types(O2_B_L %>%
                      bind_rows() %>%
                      select(ID)) %>%
      spread_draws(beta_mu, beta_sigma) %>%
      mutate(beta = rnorm(n(), beta_mu, beta_sigma), # B_7_L had to be removed,
             ID = "B_7_L" %>% fct()) %>% # so is estimated using hyperparameters
      select(-c(beta_mu, beta_sigma))
    ) %>%
  bind_rows(
    O2_B_D_samples %>%
      recover_types(O2_B_D %>%
                      bind_rows() %>%
                      select(ID)) %>%
      spread_draws(beta[ID]) %>%
      ungroup() %>%
      bind_rows(
        O2_B_D_samples %>%
          recover_types(O2_B_D %>%
                          bind_rows() %>%
                          select(ID)) %>%
          spread_draws(beta_mu, beta_sigma) %>%
          mutate(beta = rnorm(n(), beta_mu, beta_sigma),
                 ID = "B_7_D" %>% fct()) %>% # see above
          select(-c(beta_mu, beta_sigma))
      )
  ) %>%
  mutate(Round = ID %>% 
           str_split(pattern = "_", n = 2) %>% 
           map_chr(2) %>% fct())

# 4.1.3 Add confounders ####
O2_sample_estimates %<>%
  left_join(
    T_samples %>%
      recover_types(T_data %>% select(Round)) %>%
      spread_draws(T_mu[Round]) %>%
      ungroup(),
    by = c("Round", ".chain", ".iteration", ".draw")
  ) %>%
  left_join(P, by = "Round") %>%
  left_join(S, by = "Round")

O2_blank_estimates %<>%
  left_join(
    T_samples %>%
      recover_types(T_data %>% select(Round)) %>%
      spread_draws(T_mu[Round]) %>%
      ungroup(),
    by = c("Round", ".chain", ".iteration", ".draw")
  ) %>%
  left_join(P, by = "Round") %>%
  left_join(S, by = "Round")

# 4.1.4 Check success of merge ####
str(O2_sample_estimates)
O2_sample_estimates %>%
  group_by(Date, Round, ID) %>%
  summarise(beta_mean = mean(beta),
            beta_sd = sd(beta),
            O2_0_mean = mean(O2_0),
            O2_0_sd = sd(O2_0),
            T_mu_mean = mean(T_mu),
            T_mu_sd = sd(T_mu),
            Pressure = mean(Pressure),
            Salinity = mean(Salinity),
            n = length(beta)) %>%
  print(n = 96)
# looks fine

str(O2_blank_estimates)
O2_blank_estimates %>%
  group_by(Round, ID) %>%
  summarise(beta_mean = mean(beta),
            beta_sd = sd(beta),
            T_mu_mean = mean(T_mu),
            T_mu_sd = sd(T_mu),
            Pressure = mean(Pressure),
            Salinity = mean(Salinity),
            n = length(beta)) %>%
  print(n = 34)
# looks fine (note greater sd for round 7)

# 4.2 Volume ####
# 4.2.1 Prepare data ####
V <- incu_meta %>%
  filter(!is.na(Volume)) %>%
  mutate(ID = ID %>% fct(),
         Species = if_else(is.na(Species), "Blank", Species) %>% fct(),
         Treatment = Treatment %>% fct()) %>%
  select(ID, Species, Treatment, Volume)
str(V)
# Volume is given as weigh of seawater in grams, later to be converted with temperature and salinity
# measured when water was weighed:
S_end <- incu_meta %>%
  filter(!ID %>% str_detect("^0")) %>%
  mutate(ID = ID %>% fct()) %>%
  select(ID, Salinity) %>%
  rename(Salinity_end = Salinity)

# Missing volume data for the first two rounds need to be substituted with modelled averages. Blanks 
# should logically have larger volumes, so should be modelled separately. Also, gas bubbles can displace
# water in light, so light and dark should be modelled separately. Species might also displace different
# volumes despite haphazard choice of differently sized individuals. This can all be done in one model.

# 4.2.2 Prior simulation ####
# By weighing ultrapure water, I know that the empty incubation jar volume is 175 mL.

ggplot() +
  geom_density(aes(rgamma(1e5, 175^2 / 10^2, 175 / 10^2))) + # 10 mL seems like a reasonable sd
  theme_minimal() +
  theme(panel.grid = element_blank())

# 4.2.3 Run model ####
V_stan <- "
data{
  int n;
  vector<lower=0>[n] Volume;
  array[n] int Species;
  int n_Species;
  array[n] int Treatment;
  int n_Treatment;
}

parameters{
  matrix<lower=0>[n_Species, n_Treatment] V_mu;
  real<lower=0> sigma;
}

model{
  to_vector(V_mu) ~ gamma( 175^2 / 10^2 , 175 / 10^2 );
  sigma ~ exponential( 1 );
  
  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = V_mu[Species[i], Treatment[i]];
  }

  // Likelihood
  Volume ~ normal( mu , sigma );
}
"

V_mod <- V_stan %>%
  write_stan_file() %>%
  cmdstan_model()

V_samples <- V_mod$sample(
    data = V %>%
      select(Volume, Species, Treatment) %>%
      compose_data(),
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4,
    iter_sampling = 1e4)

# 4.2.4 Model checks ####
V_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

V_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# chains look good

# 4.2.5 Prior-posterior comparison ####
# sample priors
V_prior <- prior_samples(
  model = V_mod,
  data = V %>%
      select(Volume, Species, Treatment) %>%
      compose_data(),
  chains = 8, samples = 1e4)

# plot prior-posterior comparison
V_prior %>%
  prior_posterior_draws(posterior_samples = V_samples,
                        group = V %>% select(Species, Treatment),
                        parameters = c("V_mu[Species, Treatment]", "sigma"),
                        format = "long") %>%
  ggplot() +
    geom_density(aes(.value, alpha = distribution),
                 colour = NA, fill = "black") +
    scale_alpha_manual(values = c(0.6, 0.2)) +
    ggh4x::facet_nested(Treatment ~ .variable + Species, 
                        scales = "free", nest_line = TRUE) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# posteriors are very well constrained

# 4.2.6 Add volume to estimates ####
# Where volume was actually measured, the exact measurement can be matched, while where 
# volume was not measured samples from the posterior corresponding to the respective
# Species and Treatment have to be matched. End salinity must also be matched to convert
# volume from g to L.

O2_sample_estimates %<>%
  mutate(Species = case_when( # first Species and Treatment need to be added to estimates
            ID %>% str_detect("A") ~ "Amphiroa anceps",
            ID %>% str_detect("J") ~ "Jania rosea",
            ID %>% str_detect("M") ~ "Metamastophora flabellata"
            ) %>% fct(),
         Treatment = case_when(
            ID %>% str_detect("L") ~ "Light",
            ID %>% str_detect("D") ~ "Dark"
            ) %>% fct()
         ) %>%
  left_join(V, by = c("ID", "Species", "Treatment")) %>% # then exact volumes can be matched
  left_join( # and modelled volumes to fill the gaps
    V_samples %>%
      recover_types(V %>% select(Species, Treatment)) %>%
      spread_draws(V_mu[Species, Treatment]) %>%
      ungroup(),
    by = c("Species", "Treatment", ".chain", ".iteration", ".draw")
  ) %>% 
  left_join(S_end, by = "ID") %>% # and finally salinities to convert volumes
  mutate(ID = ID %>% fct_drop(),
         Species = Species %>% fct_drop())

str(O2_sample_estimates)

O2_blank_estimates %<>%
  mutate(Species = "Blank" %>% fct(),
         Treatment = case_when(
            ID %>% str_detect("L") ~ "Light",
            ID %>% str_detect("D") ~ "Dark"
            ) %>% fct()
         ) %>%
  left_join(V, by = c("ID", "Species", "Treatment")) %>%
  left_join(
    V_samples %>%
      recover_types(V %>% select(Species, Treatment)) %>%
      spread_draws(V_mu[Species, Treatment]) %>%
      ungroup(),
    by = c("Species", "Treatment", ".chain", ".iteration", ".draw")
  ) %>% 
  left_join(S_end, by = "ID") %>% 
  mutate(ID = ID %>% fct_drop(),
         Species = Species %>% fct_drop())

str(O2_blank_estimates)

# 4.2.7 Convert volume to L ####
require(seacarb)
O2_sample_estimates %<>%
  mutate(Volume = coalesce(Volume, V_mu)) %>%
  select(-V_mu) %>% # merge Volume into one column
  # seacarb::rho calculates seawater density (g/L), so if we divide g by rho we get L
  mutate(g_L = rho(S = Salinity_end, T = T_mu) %>% as.numeric(),
         Volume_L = Volume / g_L)
  
str(O2_sample_estimates)

O2_blank_estimates %<>%
  mutate(Volume = coalesce(Volume, V_mu)) %>%
  select(-V_mu) %>%
  mutate(g_L = rho(S = Salinity_end, T = T_mu) %>% as.numeric(),
         Volume_L = Volume / g_L)

str(O2_blank_estimates)

# 4.2.8 Correct for volume ####
# Multiplying beta by Volume_L converts it from µM min^-1 to µmol min^-1 because
# µM min^-1 is equivalent to µmol L^-1 min^-1, which * L leaves µmol min^-1.
O2_sample_estimates %<>%
  mutate(beta_µmol_min = beta * Volume_L)

O2_blank_estimates %<>%
  mutate(beta_µmol_min = beta * Volume_L)

# 4.3 Blank correction ####
# 4.3.1 Combine estimates ####
O2_estimates <- 
  O2_sample_estimates %>%
  select(-c(beta, Volume, Salinity_end, g_L, Volume_L)) %>%
  left_join(
    O2_blank_estimates %>%
      select(-c(Species, beta, Volume, Salinity_end, g_L, Volume_L)) %>%
      rename(beta_µmol_min_blank = beta_µmol_min,
             ID_blank = ID),
    by = c(".chain", ".iteration", ".draw", "Round", 
           "T_mu", "Pressure", "Salinity", "Treatment")
  )

# 4.3.2 Check success of merge ####
O2_estimates %>%
  group_by(Round, ID, ID_blank, Species, Treatment) %>%
  summarise(beta_µmol_min_mean = mean(beta_µmol_min),
            beta_µmol_min_sd = sd(beta_µmol_min),
            beta_µmol_min_blank_mean = mean(beta_µmol_min_blank),
            beta_µmol_min_blank_sd = sd(beta_µmol_min_blank),
            n = length(beta_µmol_min)) %>%
  print(n = 96)

# 4.3.3 Subtract blank slope ####
O2_estimates %<>%
  mutate(beta_µmol_min_corrected = beta_µmol_min - beta_µmol_min_blank)

# 4.4 Mass ####
# 4.4.1 Add mass estimates ####
mass <- read.csv("Mass.csv") %>%
  select(ID, Species, DM) %>%
  mutate(ID = ID %>% fct(),
         Species = Species %>% fct()) %>%
  rename(Individual = ID,
         Mass = DM)

O2_estimates %<>%
  mutate(Individual = ID %>% str_sub(end = -3) %>% fct()) %>% # remove Treatment from ID
  left_join(mass, by = c("Individual", "Species")) # join mass by short ID

# 4.4.2 Correct for mass ####
# O2 evolution is currently given as µmol min^-1, but I want it as µmol g^-1 h^-1, so
# I need to divide by dry mass (g) and multiply by 60 min h^-1.
O2_estimates %<>%
  mutate(beta_µmol_g_h = beta_µmol_min_corrected / Mass * 60)

# 4.5 Photosynthesis and respiration ####
# 4.5.1 Check data ####
O2_estimates %>%
  group_by(Round, ID, Species, Treatment) %>%
  summarise(beta_µmol_g_h_mean = mean(beta_µmol_g_h),
            beta_µmol_g_h_sd = sd(beta_µmol_g_h),
            O2_0_mean = mean(O2_0),
            O2_0_sd = sd(O2_0),
            T_mu_mean = mean(T_mu),
            T_mu_sd = sd(T_mu),
            Pressure = unique(Pressure),
            Salinity = unique(Salinity),
            n = length(beta_µmol_g_h)) %>%
  print(n = 96)

# 4.5.2 Pivot data and calculate gross photosynthesis ####
O2_estimates %<>%
  mutate(ID = ID %>% str_sub(end = -3) %>% fct(), # remove Treatment from ID
         Round = Round %>% str_sub(end = -3) %>% fct(), # remove Treatment from Round
         Individual = Individual %>% str_sub(start = 3) %>% fct()) %>% # remove Species from Individual
  select(-c(beta_µmol_min, beta_µmol_min_blank, 
            beta_µmol_min_corrected, ID_blank)) %>%
  pivot_wider(names_from = Treatment, 
              values_from = c(beta_µmol_g_h, O2_0, T_mu,
                              Pressure, Salinity)) %>%
  rename(nP = beta_µmol_g_h_Light,
         R = beta_µmol_g_h_Dark,
         O2_L = O2_0_Light,
         O2_D = O2_0_Dark,
         T_L = T_mu_Light,
         T_D = T_mu_Dark,
         P_L = Pressure_Light,
         P_D = Pressure_Dark,
         S_L = Salinity_Light,
         S_D = Salinity_Dark) %>%
  mutate(gP = nP + R, # add net photosynthesis and respiration to get gross photosynthesis
         O2_LD = (O2_L + O2_D) / 2, # calculate the mean of all confounders that vary between
         T_LD = (T_L + T_D) / 2, # light and dark incubations
         P_LD = (P_L + P_D) / 2,
         S_LD = (S_L + S_D) / 2)

str(O2_estimates)

# 4.5.3 Clean up ####
rm(list = setdiff(ls(), "O2_estimates"))

# 4.5.4 Save estimates ####
O2_estimates %>% write_rds("O2_estimates.rds")