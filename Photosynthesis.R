# 1. Load data ####
require(tidyverse)
O2_estimates <- read_rds("O2_estimates.rds")
C_estimates <- read_rds("C_estimates.rds")

# Combine O2 and C estimates
O2_C_estimates <- C_estimates %>%
  select(starts_with("."), ID, Round, Species, Individual, G, D, nP, R, gP) %>%
  rename(nP_CO2 = nP, R_CO2 = R, gP_CO2 = gP) %>%
  full_join(O2_estimates %>%
              select(starts_with("."), ID, Round, Species, Individual, nP, R, gP) %>%
              rename(nP_O2 = nP, R_O2 = R, gP_O2 = gP),
            by = c(".chain", ".iteration", ".draw", "ID", "Round", "Species", "Individual")) 

# Load light data
PAR <- read.csv("Incubation.csv") %>%
  filter(Treatment == "Light" & !is.na(Species) & !is.na(PAR)) %>%
  mutate(Species = Species %>% fct(),
         Round = Round %>% as.character() %>% fct(),
         Treatment = Treatment %>% fct(),
         ID = ID %>% str_sub(end = -3) %>% fct(),
         Individual = ID %>% str_sub(start = 3) %>% fct(),
         PAR = PAR %>% as.numeric()) %>%
  select(Round, ID, Species, Individual, PAR)

# Add light data
require(magrittr)
O2_C_estimates %<>%
  left_join(PAR, by = c("Round", "ID", "Species", "Individual"))

# Summarise data
O2_C_estimates_summary <- O2_C_estimates %>%
  group_by(Round, ID, Species, Individual) %>%
  summarise(
    across(
      .cols = c(G, D, nP_O2, nP_CO2, R_O2, R_CO2, gP_O2, gP_CO2, PAR),
      .fns = list(mean = mean, sd = sd)
    )
  ) %>%
  ungroup() %>%
  keep(~ !all(.x == 0 | is.na(.x))) # remove columns that contain only zeros or NAs

# Check summary
O2_C_estimates_summary %>%
  select(Round, ID, Species, Individual, G_mean, D_mean,
         nP_O2_mean, nP_CO2_mean, R_O2_mean, R_CO2_mean) %>%
  print(n = 51)

# 2. Photosynthesis-irradiance ####
# Visualise
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

Fig_S4a <- O2_C_estimates %>%
  drop_na(nP_O2) %>%
  ggplot() +
    geom_violin(aes(PAR, nP_O2, fill = Species, colour = Species, group = ID),
                position = "identity", width = 50, alpha = 0.5) +
    scale_colour_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"),
                        guide = "none") +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"),
                      guide = "none") +
    facet_grid(~ Species) +
    labs(x = expression("PAR (µmol photons m"^-2*" s"^-1*")"),
         y = expression("P"["n"]*" (µmol O"[2]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(0, 800), ylim = c(0, 200), 
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(strip.text = element_text(face = "italic"))
    
Fig_S4a %>%
  ggsave(filename = "Fig_S4a.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

Fig_S4b <- O2_C_estimates %>%
  drop_na(nP_CO2) %>%
  ggplot() +
    geom_violin(aes(PAR, nP_CO2, fill = Species, colour = Species, group = ID),
                position = "identity", width = 50, alpha = 0.5) +
    scale_colour_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"),
                        guide = "none") +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"),
                      guide = "none") +
    facet_grid(~ Species) +
    labs(x = expression("PAR (µmol photons m"^-2*" s"^-1*")"),
         y = expression("P"["n"]*" (µmol CO"[2]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(0, 800), ylim = c(0, 100), 
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(strip.text = element_text(face = "italic"))
    
Fig_S4b %>%
  ggsave(filename = "Fig_S4b.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# 3. Photosynthetic quotient ####
# 3.1 Visualise data ####
require(ggdensity)
Fig_1 <- O2_C_estimates %>%
  drop_na(nP_O2, nP_CO2) %>%
  ggplot() +
    geom_hdr(aes(nP_CO2, nP_O2, fill = Species, group = ID),
             alpha = 0.5, n = 600, method = "mvnorm", probs = 0.999) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d")) +
    labs(x = expression("P"["n"]*" (µmol CO"[2]*" g"^-1*" h"^-1*")"),
         y = expression("P"["n"]*" (µmol O"[2]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(0, 80), ylim = c(0, 160),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

Fig_1 %>%
  ggsave(filename = "Fig_1.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# 3.2 Prepare data ####
# The model needs to be conditioned on a data summary and standard deviations
# are passed to it to account for measurement error.
PQ_data <- O2_C_estimates_summary %>%
  drop_na(nP_O2_mean, nP_CO2_mean) %>%
  select(Species, nP_CO2_mean, nP_CO2_sd, nP_O2_mean, nP_O2_sd)

# 3.3 Prior simulation ####
# In the regression of net O2 production on net CO2 fixation, the intercept
# is expected to be zero (no O2 production without CO2 fixation and vice versa) 
# and the slope is the photosynthetic quotient (PQ), which has to be positive 
# (truncated prior) and is usually assumed to be 1, i.e. one O2 produced for every 
# CO2 taken up.
require(truncnorm) # simulate from the truncated normal
tibble(n = 1:1e3, # simulate hierachical prior
       PQ_mu = rtruncnorm( 1e3 , mean = 1 , sd = 0.5 , a = 0 ), # a is lower bound
       PQ_sigma = rexp( 1e3 , 1 ),
       PQ = rtruncnorm( 1e3 , mean = PQ_mu , sd = PQ_sigma , a = 0 )) %>%
  expand_grid(nP_CO2 = PQ_data %$% seq(min(nP_CO2_mean), max(nP_CO2_mean))) %>%
  mutate(nP_O2 = PQ * nP_CO2) %>%
  ggplot(aes(nP_CO2, nP_O2, group = n)) +
    geom_hline(yintercept = PQ_data %$% c(0, max(nP_O2_mean))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(ylim = c(0, 180), expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# I also need a prior for estimates of the true predictor
PQ_data %$% max(nP_CO2_mean) %>% round / 2 # I'll take half the predictor range
ggplot() + # and choose a reasonable prior centred on it
  geom_density(aes(rnorm(n = 1e5, mean = 36, sd = 15))) +
  geom_vline(xintercept = PQ_data %$% c(min(nP_CO2_mean), max(nP_CO2_mean))) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# 3.4 Run model ####
PQ_stan <- "
data{
  int n;
  vector[n] nP_CO2_mean;
  vector<lower=0>[n] nP_CO2_sd;
  vector[n] nP_O2_mean;
  vector<lower=0>[n] nP_O2_sd;
  array[n] int Species;
  int n_Species;
}

parameters{
  // True estimates for measurement error
  vector[n] nP_CO2;
  vector[n] nP_O2;

  // Hyperparameters
  real<lower=0> PQ_mu;
  real<lower=0> PQ_sigma;
  
  // Species-specific parameters
  vector<lower=0>[n_Species] PQ;
  
  // Likelihood uncertainty
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  PQ_mu ~ normal( 1 , 0.5 ) T[0,];
  PQ_sigma ~ exponential( 1 );
  
  // Species-specific priors
  PQ ~ normal( PQ_mu , PQ_sigma ) T[0,];
  
  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );
  
  // Estimated predictor
  nP_CO2 ~ normal( 36 , 15 );
  nP_CO2_mean ~ normal( nP_CO2 , nP_CO2_sd );
  
  // Model
  vector[n] mu;
  for ( i in 1:n ) {
      mu[i] = PQ[Species[i]] * nP_CO2[i];
  }

  // Likelihood with measurement error
  nP_O2 ~ normal( mu , sigma );
  nP_O2_mean ~ normal( nP_O2 , nP_O2_sd );
}
"

require(cmdstanr)
PQ_mod <- PQ_stan %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
PQ_samples <- PQ_mod$sample(
  data = PQ_data %>% compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4,
  adapt_delta = 0.99,
  max_treedepth = 14) 
# increased adapt_delta and max_treedepth to make sampler go slower and reduce divergences

# 3.5 Model checks ####
PQ_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

require(bayesplot)
PQ_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "PQ_rank.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

# 3.6 Prior-posterior comparison ####
source("functions.R")
# sample priors
PQ_prior <- prior_samples(
  model = PQ_mod,
  data = PQ_data %>% compose_data(),
  chains = 8, samples = 1e4)
# Stan struggles to sample joint prior distributions of hierarchical models,
# but it's good enough for a quick prior-posterior check.

# plot prior-posterior comparison for main parameters
PQ_prior %>%
  prior_posterior_draws(posterior_samples = PQ_samples,
                        group = PQ_data %>% select(Species),
                        parameters = c("PQ[Species]", "PQ_mu", "PQ_sigma", "sigma"),
                        format = "long") %>%
  prior_posterior_plot(group_name = "Species", ridges = FALSE)
# posteriors are reasonably well constrained

# check prior-posterior of estimates of true predictor
ggsave(
  PQ_samples %>%
    spread_draws(nP_CO2[n]) %>%
    mutate(nP_CO2_prior = rnorm(max(.draw), 36 , 15)) %>%
    pivot_longer(cols = c(nP_CO2, nP_CO2_prior),
                 values_to = "nP_CO2", names_to = "distribution") %>%
    mutate(distribution = if_else(distribution %>% str_detect("prior"),
                                  "prior", "posterior") %>% fct()) %>%
    ggplot() +
      geom_density(aes(nP_CO2, alpha = distribution),
                   fill = "black", colour = NA) +
      scale_alpha_manual(values = c(0.6, 0.2)) +
      facet_wrap(~ n, scales = "free") +
      theme_minimal() +
      theme(panel.grid = element_blank()),
  filename = "PQ_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)

# 3.7 Predictions ####
# Simulate smooth prior distributions using R
PQ_prior <- tibble(
  PQ_mu = rtruncnorm( 8 * 1e4 , mean = 1 , sd = 0.5 , a = 0 ), 
  PQ_sigma = rexp( 8 * 1e4 , 1 ),
  PQ = rtruncnorm( 8 * 1e4 , mean = PQ_mu , sd = PQ_sigma , a = 0 ),
  sigma = rexp( 8 * 1e4 , 1 )
)

PQ_prior %>%
  ggplot() +
    geom_density(aes(PQ_mu), colour = "black") +
    geom_density(aes(PQ), colour = "grey") +
    scale_x_continuous(limits = c(0, 10), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# Simulate smooth posteriors with hyperparameters
PQ_hyperposterior <- PQ_samples %>%
  spread_draws(PQ_mu, PQ_sigma, sigma) %>%
  mutate(PQ = rtruncnorm( n() , mean = PQ_mu , sd = PQ_sigma , a = 0 )) # simulate for unobserved corallines

PQ_hyperposterior %>%
  ggplot() +
    geom_density(aes(PQ_mu), colour = "black") +
    geom_density(aes(PQ), colour = "grey") +
    scale_x_continuous(limits = c(0, 10), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# Join priors and posteriors
PQ_prior_posterior <- PQ_samples %>%
  recover_types(PQ_data %>% select(Species)) %>%
  spread_draws(PQ[Species], sigma) %>%
  ungroup() %>%
  bind_rows(
    PQ_hyperposterior %>% 
      mutate(Species = "Branching corallines" %>% fct()) %>%
      select(-c(PQ_mu, PQ_sigma)),
    PQ_prior %>%
      mutate(.chain = 1:8 %>% rep(each = 1e4),
             .iteration = 1:1e4 %>% rep(8),
             .draw = 1:(8*1e4),
             Species = "Prior" %>% fct()) %>%
      select(-c(PQ_mu, PQ_sigma))
  )

str(PQ_prior_posterior)

# Calculate probability of PQ > 1
PQ_P <- PQ_prior_posterior %>%
  group_by(Species) %>%
  summarise(mean = mean(PQ),
            median = median(PQ),
            mode = mode_qi(PQ)$y,
            sd = sd(PQ),
            P = mean(PQ > 1),
            n = length(PQ)) %>%
  mutate(P_percentage = signif(P * 100, digits = 2) %>%
           str_c("%"))

# Plot PQ
require(ggdist)
Fig_1b <- PQ_prior_posterior %>%
  mutate(Species = Species %>% 
           fct_rev() %>%
           fct_relevel("Prior", after = Inf)) %>%
  ggplot() +
    stat_slab(aes(y = Species, x = PQ, fill = Species),
              height = 5, alpha = 0.5) +
    geom_vline(xintercept = 1) +
    geom_text(data = PQ_P, aes(x = 3, y = Species, label = P_percentage),
              hjust = 1, vjust = -0.5, size = 3.5, family = "Futura") +
    scale_fill_manual(values = c("#363538", "#6b4d8d", "#f5a54a", "#ec7faa", "#b5b8ba"),
                      guide = guide_legend(reverse = TRUE)) +
    scale_x_continuous(breaks = seq(0, 3, 0.5),
                       labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1, 0.1, 1))) +
    labs(x = expression("PQ (µmol O"[2]*" µmol"^-1*" CO"[2]*")")) +
    coord_cartesian(xlim = c(0, 3), ylim = c(1, 6.5),
                    expand = FALSE, clip = "on") +
    mytheme +
    theme(legend.text = element_text(face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank())

# Spread across predictor
PQ_predictions <- PQ_prior_posterior %>%
  left_join(PQ_data %>%
              group_by(Species) %>%
              summarise(min = min(nP_CO2_mean),
                        max = max(nP_CO2_mean)),
            by = "Species") %>%
  mutate(min = if_else(is.na(min),
                       PQ_data %$% min(nP_CO2_mean),
                       min),
         max = if_else(is.na(max),
                       PQ_data %$% max(nP_CO2_mean),
                       max)) %>%
  rowwise() %>%
  mutate(nP_CO2_mean = list( seq(min, max, length.out = 100) )) %>%
  select(-c(min, max)) %>%
  unnest(nP_CO2_mean) %>%
  mutate(mu = PQ * nP_CO2_mean,
         obs = rnorm( n() , mu , sigma ))

PQ_predictions_summary <- PQ_predictions %>%
  group_by(Species, nP_CO2_mean) %>%
  reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
          obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
  unnest(c(mu, obs), names_sep = "_")

require(geomtextpath)
Fig_1a <- PQ_predictions_summary %>%
  ggplot() +
    geom_textabline(slope = 1, label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
    geom_hdr(data = O2_C_estimates %>% drop_na(nP_O2, nP_CO2),
             aes(nP_CO2, nP_O2, fill = Species, group = ID),
             n = 600, method = "mvnorm", probs = 0.999) +
    # geom_ribbon(data = . %>% filter(Species == "Branching corallines" & mu_.width == 0.9),
    #             aes(nP_CO2_mean, ymin = mu_ymin, ymax = mu_ymax), 
    #             colour = NA, fill = "#363538", alpha = 0.3) +
    geom_line(data = . %>% filter(!Species %in% c("Prior", "Branching corallines")),
              aes(nP_CO2_mean, mu_y, colour = Species)) +
    geom_ribbon(data = . %>% filter(!Species %in% c("Prior", "Branching corallines")),
                aes(nP_CO2_mean, ymin = mu_ymin, ymax = mu_ymax,
                    fill = Species, alpha = factor(mu_.width))) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"), guide = "none") +
    scale_colour_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.5, 0.4, 0.3), guide = "none") + # geom_hdr requires additional alpha level
    labs(x = expression("P"["n"]*" (µmol CO"[2]*" g"^-1*" h"^-1*")"),
         y = expression("P"["n"]*" (µmol O"[2]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(0, 80), ylim = c(0, 160),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

require(patchwork)
Fig_1 <- ( Fig_1a | Fig_1b ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

Fig_1 %>%
  ggsave(filename = "Fig_1mod.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# 4. Respiratory quotient ####
# 4.1 Visualise data ####
Fig_2 <- O2_C_estimates %>%
  drop_na(R_O2, R_CO2) %>%
  ggplot() +
    geom_hdr(aes(R_O2, R_CO2, fill = Species, group = ID),
             alpha = 0.5, n = 600, method = "mvnorm", probs = 0.999) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d")) +
    labs(x = expression("R (µmol O"[2]*" g"^-1*" h"^-1*")"),
         y = expression("R (µmol CO"[2]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(0, 15), ylim = c(0, 20),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

Fig_2 %>%
  ggsave(filename = "Fig_2.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# 4.2 Prepare data ####
RQ_data <- O2_C_estimates_summary %>%
  drop_na(R_O2_mean, R_CO2_mean) %>%
  select(Species, R_CO2_mean, R_CO2_sd, R_O2_mean, R_O2_sd)

# 4.3 Prior simulation ####
# In the regression of dark CO2 production on O2 consumption, the intercept
# is expected to be zero (no CO2 production without O2 consumption and vice versa) 
# and the slope is the respiratory quotient (RQ), which has to be positive 
# (truncated prior) and is usually assumed to be 1, i.e. one CO2 produced for every 
# O2 taken up.
tibble(n = 1:1e3, # simulate hierachical prior
       RQ_mu = rtruncnorm( 1e3 , mean = 1 , sd = 0.5 , a = 0 ), # a is lower bound
       RQ_sigma = rexp( 1e3 , 1 ),
       RQ = rtruncnorm( 1e3 , mean = RQ_mu , sd = RQ_sigma , a = 0 )) %>%
  expand_grid(R_O2 = RQ_data %$% seq(min(R_O2_mean), max(R_O2_mean))) %>%
  mutate(R_CO2 = RQ * R_O2) %>%
  ggplot(aes(R_O2, R_CO2, group = n)) +
    geom_hline(yintercept = RQ_data %$% c(0, max(R_CO2_mean))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(ylim = c(0, 25), expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# I also need a prior for estimates of the true predictor
RQ_data %$% max(R_O2_mean) %>% round / 2 # I'll take half the predictor range
ggplot() + # and choose a reasonable prior centred on it
  geom_density(aes(rnorm(n = 1e5, mean = 6.5, sd = 2.5))) +
  geom_vline(xintercept = RQ_data %$% c(min(R_O2_mean), max(R_O2_mean))) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# 4.4 Run model ####
RQ_stan <- "
data{
  int n;
  vector[n] R_O2_mean;
  vector<lower=0>[n] R_O2_sd;
  vector[n] R_CO2_mean;
  vector<lower=0>[n] R_CO2_sd;
  array[n] int Species;
  int n_Species;
}

parameters{
  // True estimates for measurement error
  vector[n] R_O2;
  vector[n] R_CO2;

  // Hyperparameters
  real<lower=0> RQ_mu;
  real<lower=0> RQ_sigma;
  
  // Species-specific parameters
  vector<lower=0>[n_Species] RQ;
  
  // Likelihood uncertainty
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  RQ_mu ~ normal( 1 , 0.5 ) T[0,];
  RQ_sigma ~ exponential( 1 );
  
  // Species-specific priors
  RQ ~ normal( RQ_mu , RQ_sigma ) T[0,];
  
  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );
  
  // Estimated predictor
  R_O2 ~ normal( 6.5 , 2.5 );
  R_O2_mean ~ normal( R_O2 , R_O2_sd );
  
  // Model
  vector[n] mu;
  for ( i in 1:n ) {
      mu[i] = RQ[Species[i]] * R_O2[i];
  }

  // Likelihood with measurement error
  R_CO2 ~ normal( mu , sigma );
  R_CO2_mean ~ normal( R_CO2 , R_CO2_sd );
}
"

RQ_mod <- RQ_stan %>%
  write_stan_file() %>%
  cmdstan_model()

RQ_samples <- RQ_mod$sample(
  data = RQ_data %>% compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4,
  adapt_delta = 0.999,
  max_treedepth = 15) 
# increased adapt_delta max_treedepth to make sampler go slower and reduce divergences

# 4.5 Model checks ####
RQ_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# <1% rhat above 1.001
# good effective sample size

RQ_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "RQ_rank.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

# 4.6 Prior-posterior comparison ####
# sample priors
RQ_prior <- prior_samples(
  model = RQ_mod,
  data = RQ_data %>% compose_data(),
  chains = 8, samples = 1e4)

# plot prior-posterior comparison for main parameters
RQ_prior %>%
  prior_posterior_draws(posterior_samples = RQ_samples,
                        group = RQ_data %>% select(Species),
                        parameters = c("RQ[Species]", "RQ_mu", "RQ_sigma", "sigma"),
                        format = "long") %>%
  prior_posterior_plot(group_name = "Species", ridges = FALSE)
# posteriors are reasonably well constrained

# check prior-posterior of estimates of true predictor
ggsave(
  RQ_samples %>%
    spread_draws(R_O2[n]) %>%
    mutate(R_O2_prior = rnorm(max(.draw), 6.5 , 2.5)) %>%
    pivot_longer(cols = c(R_O2, R_O2_prior),
                 values_to = "R_O2", names_to = "distribution") %>%
    mutate(distribution = if_else(distribution %>% str_detect("prior"),
                                  "prior", "posterior") %>% fct()) %>%
    ggplot() +
      geom_density(aes(R_O2, alpha = distribution),
                   fill = "black", colour = NA) +
      scale_alpha_manual(values = c(0.6, 0.2)) +
      facet_wrap(~ n, scales = "free") +
      theme_minimal() +
      theme(panel.grid = element_blank()),
  filename = "RQ_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)

# 4.7 Predictions ####
# Simulate smooth prior distributions using R
RQ_prior <- tibble(
  RQ_mu = rtruncnorm( 8 * 1e4 , mean = 1 , sd = 0.5 , a = 0 ), 
  RQ_sigma = rexp( 8 * 1e4 , 1 ),
  RQ = rtruncnorm( 8 * 1e4 , mean = RQ_mu , sd = RQ_sigma , a = 0 ),
  sigma = rexp( 8 * 1e4 , 1 )
)

RQ_prior %>%
  ggplot() +
    geom_density(aes(RQ_mu), colour = "black") +
    geom_density(aes(RQ), colour = "grey") +
    scale_x_continuous(limits = c(0, 10), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# Simulate smooth posteriors with hyperparameters
RQ_hyperposterior <- RQ_samples %>%
  spread_draws(RQ_mu, RQ_sigma, sigma) %>%
  mutate(RQ = rtruncnorm( n() , mean = RQ_mu , sd = RQ_sigma , a = 0 )) # simulate for unobserved corallines

RQ_hyperposterior %>%
  ggplot() +
    geom_density(aes(RQ_mu), colour = "black") +
    geom_density(aes(RQ), colour = "grey") +
    scale_x_continuous(limits = c(0, 10), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# Join priors and posteriors
RQ_prior_posterior <- RQ_samples %>%
  recover_types(RQ_data %>% select(Species)) %>%
  spread_draws(RQ[Species], sigma) %>%
  ungroup() %>%
  bind_rows(
    RQ_hyperposterior %>% 
      mutate(Species = "Branching corallines" %>% fct()) %>%
      select(-c(RQ_mu, RQ_sigma)),
    RQ_prior %>%
      mutate(.chain = 1:8 %>% rep(each = 1e4),
             .iteration = 1:1e4 %>% rep(8),
             .draw = 1:(8*1e4),
             Species = "Prior" %>% fct()) %>%
      select(-c(RQ_mu, RQ_sigma))
  )

str(RQ_prior_posterior)

# Calculate probability of RQ < 1
RQ_P <- RQ_prior_posterior %>%
  group_by(Species) %>%
  summarise(mean = mean(RQ),
            median = median(RQ),
            mode = mode_qi(RQ)$y,
            sd = sd(RQ),
            P = mean(RQ < 1),
            n = length(RQ)) %>%
  mutate(P_percentage = signif(P * 100, digits = 2) %>%
           str_c("%"))

Fig_2b <- RQ_prior_posterior %>%
  mutate(Species = Species %>% 
           fct_rev() %>%
           fct_relevel("Prior", after = Inf)) %>%
  ggplot() +
    stat_slab(aes(y = Species, x = RQ, fill = Species),
              height = 5, alpha = 0.5) +
    geom_vline(xintercept = 1) +
    geom_text(data = RQ_P, aes(x = 3, y = Species, label = P_percentage),
              hjust = 1, vjust = -0.5, size = 3.5, family = "Futura") +
    scale_fill_manual(values = c("#363538", "#6b4d8d", "#f5a54a", "#ec7faa", "#b5b8ba"),
                      guide = guide_legend(reverse = TRUE)) +
    scale_x_continuous(breaks = seq(0, 3, 0.5),
                       labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1, 0.1, 1))) +
    labs(x = expression("RQ (µmol CO"[2]*" µmol"^-1*" O"[2]*")")) +
    coord_cartesian(xlim = c(0, 3), ylim = c(1, 7),
                    expand = FALSE, clip = "on") +
    mytheme +
    theme(legend.text = element_text(face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank())

# Spread across predictor
RQ_predictions <- RQ_prior_posterior %>%
  left_join(RQ_data %>%
              group_by(Species) %>%
              summarise(min = min(R_O2_mean),
                        max = max(R_O2_mean)),
            by = "Species") %>%
  mutate(min = if_else(is.na(min),
                       RQ_data %$% min(R_O2_mean),
                       min),
         max = if_else(is.na(max),
                       RQ_data %$% max(R_O2_mean),
                       max)) %>%
  rowwise() %>%
  mutate(R_O2_mean = list( seq(min, max, length.out = 100) )) %>%
  select(-c(min, max)) %>%
  unnest(R_O2_mean) %>%
  mutate(mu = RQ * R_O2_mean,
         obs = rnorm( n() , mu , sigma ))

RQ_predictions_summary <- RQ_predictions %>%
  group_by(Species, R_O2_mean) %>%
  reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
          obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
  unnest(c(mu, obs), names_sep = "_")

Fig_2a <- RQ_predictions_summary %>%
  ggplot() +
    geom_textabline(slope = 1, label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
    geom_hdr(data = O2_C_estimates %>% drop_na(R_CO2, R_O2),
             aes(R_O2, R_CO2, fill = Species, group = ID),
             n = 600, method = "mvnorm", probs = 0.999) +
    # geom_ribbon(data = . %>% filter(Species == "Branching corallines" & mu_.width == 0.9),
    #             aes(R_O2_mean, ymin = mu_ymin, ymax = mu_ymax), 
    #             colour = NA, fill = "#363538", alpha = 0.3) +
    geom_line(data = . %>% filter(!Species %in% c("Prior", "Branching corallines")),
              aes(R_O2_mean, mu_y, colour = Species)) +
    geom_ribbon(data = . %>% filter(!Species %in% c("Prior", "Branching corallines")),
                aes(R_O2_mean, ymin = mu_ymin, ymax = mu_ymax,
                    fill = Species, alpha = factor(mu_.width))) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"), guide = "none") +
    scale_colour_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.5, 0.4, 0.3), guide = "none") + # geom_hdr requires additional alpha level
    labs(x = expression("R (µmol O"[2]*" g"^-1*" h"^-1*")"),
         y = expression("R (µmol CO"[2]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(0, 15), ylim = c(0, 20),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

Fig_2 <- ( Fig_2a | Fig_2b ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

Fig_2 %>%
  ggsave(filename = "Fig_2mod.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)