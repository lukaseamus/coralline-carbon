# 1. Load data ####
require(tidyverse)
C_estimates <- read_rds("C_estimates.rds")

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
C_estimates %<>%
  left_join(PAR, by = c("Round", "ID", "Species", "Individual"))

# Summarise data
C_estimates_summary <- C_estimates %>%
  group_by(Round, ID, Species, Individual) %>%
  summarise(
    across(
      .cols = c(G, D, nP, R, gP, PAR),
      .fns = list(mean = mean, sd = sd)
    )
  ) %>%
  ungroup() %>%
  keep(~ !all(.x == 0 | is.na(.x))) # remove columns that contain only zeros or NAs

# Check summary
C_estimates_summary %>%
  select(Round, ID, Species, Individual, G_mean, D_mean,
         nP_mean, R_mean, gP_mean) %>%
  print(n = 51)

# 2. Calcification-irradiance ####
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

Fig_S5 <- C_estimates %>%
  drop_na(G) %>%
  ggplot() +
    geom_violin(aes(PAR, G, fill = Species, colour = Species, group = ID),
                position = "identity", width = 50, alpha = 0.5) +
    scale_colour_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"),
                        guide = "none") +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"),
                      guide = "none") +
    facet_grid(~ Species) +
    labs(x = expression("PAR (µmol photons m"^-2*" s"^-1*")"),
         y = expression("G (µmol CaCO"[3]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(0, 800), ylim = c(0, 40), 
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(strip.text = element_text(face = "italic"))
    
Fig_S5 %>%
  ggsave(filename = "Fig_S5.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)


# 3. Calcification:photosynthesis ####
# 3.1 Visualise data ####
require(ggdensity)
Fig_3 <- C_estimates %>%
  drop_na(nP, G) %>%
  group_by(ID) %>% # randomly reshuffle within ID
  mutate(G = G %>% sample()) %>%
  ggplot() +
    geom_hdr(aes(nP, G, fill = Species, group = ID),
             alpha = 0.5, n = 600, method = "mvnorm", probs = 0.999) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d")) +
    labs(x = expression("P"["n"]*" (µmol CO"[2]*" g"^-1*" h"^-1*")"),
         y = expression("G (µmol CaCO"[3]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(0, 90), ylim = c(0, 40),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

Fig_3 %>%
  ggsave(filename = "Fig_3.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# 3.2 Prepare data ####
# The model needs to be conditioned on a data summary and standard deviations
# are passed to it to account for measurement error.
GP_data <- C_estimates_summary %>%
  drop_na(nP_mean, G_mean) %>%
  select(Species, nP_mean, nP_sd, G_mean, G_sd)

# 3.3 Prior simulation ####
# In the regression of calcification on net CO2 fixation, an intercept is imaginable
# (net calcification and dissolution are both possible with zero net photosynthesis)
# but is expected to be zero (no calcification without CO2 fixation) and the slope is 
# the ratio of calcification to photosynthesis (GP), which has to be positive 
# (truncated prior) and probably lies around 0.24 (doi: 10.1038/s41467-024-52697-5), 
# i.e. 0.24 CaCO3 produced for every CO2 taken up.

require(truncnorm) # simulate from the truncated normal
tibble(n = 1:1e3, # simulate hierachical prior
       GP_mu = rtruncnorm( 1e3 , mean = 0.24 , sd = 0.2 , a = 0 ), # a is lower bound
       GP_sigma = rexp( 1e3 , 1 ),
       GP = rtruncnorm( 1e3 , mean = GP_mu , sd = GP_sigma , a = 0 ),
       G0_mu = rnorm( 1e3 , 0 , 2 ),
       G0_sigma = rexp( 1e3 , 1 ),
       G0 = rnorm( 1e3 , mean = G0_mu , sd = G0_sigma )) %>%
  expand_grid(nP = GP_data %$% seq(min(nP_mean), max(nP_mean))) %>%
  mutate(G = G0 + GP * nP) %>%
  ggplot(aes(nP, G, group = n)) +
    geom_hline(yintercept = GP_data %$% c(min(G_mean), 0, max(G_mean))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(ylim = c(-5, 40), expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# I also need a prior for estimates of the true predictor
GP_data %$% max(nP_mean) %>% round / 2 # I'll take half the predictor range
ggplot() + # and choose a reasonable prior centred on it
  geom_density(aes(rnorm(n = 1e5, mean = 40.5, sd = 20))) +
  geom_vline(xintercept = GP_data %$% c(min(nP_mean), max(nP_mean))) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# 3.4 Run model ####
GP_stan <- "
data{
  int n;
  vector[n] nP_mean;
  vector<lower=0>[n] nP_sd;
  vector[n] G_mean;
  vector<lower=0>[n] G_sd;
  array[n] int Species;
  int n_Species;
}

parameters{
  // True estimates for measurement error
  vector[n] nP;
  vector[n] G;

  // Hyperparameters
  real<lower=0> GP_mu;
  real<lower=0> GP_sigma;
  real G0_mu;
  real<lower=0> G0_sigma;
  
  // Species-specific parameters
  vector<lower=0>[n_Species] GP;
  vector[n_Species] G0;
  
  // Likelihood uncertainty
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  GP_mu ~ normal( 0.24 , 0.2 ) T[0,];
  GP_sigma ~ exponential( 1 );
  G0_mu ~ normal( 0 , 2 );
  G0_sigma ~ exponential( 1 );
  
  // Species-specific priors
  GP ~ normal( GP_mu , GP_sigma ) T[0,];
  G0 ~ normal( G0_mu , G0_sigma );
  
  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );
  
  // Estimated predictor
  nP ~ normal( 40.5 , 20 );
  nP_mean ~ normal( nP , nP_sd );
  
  // Model
  vector[n] mu;
  for ( i in 1:n ) {
      mu[i] = G0[Species[i]] + GP[Species[i]] * nP[i];
  }

  // Likelihood with measurement error
  G ~ normal( mu , sigma );
  G_mean ~ normal( G , G_sd );
}
"

require(cmdstanr)
GP_mod <- GP_stan %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
GP_samples <- GP_mod$sample(
  data = GP_data %>% compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4,
  adapt_delta = 0.999,
  max_treedepth = 15) 
# increased adapt_delta and max_treedepth to make sampler go slower and reduce divergences

# 3.5 Model checks ####
GP_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# <1% rhat above 1.001
# good effective sample size

require(bayesplot)
GP_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "GP_rank.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

# 3.6 Prior-posterior comparison ####
source("functions.R")
# sample priors
GP_prior <- prior_samples(
  model = GP_mod,
  data = GP_data %>% compose_data(),
  chains = 8, samples = 1e4)
# Stan struggles to sample joint prior distributions of hierarchical models,
# but it's good enough for a quick prior-posterior check.

# plot prior-posterior comparison for main parameters
GP_prior %>%
  prior_posterior_draws(posterior_samples = GP_samples,
                        group = GP_data %>% select(Species),
                        parameters = c("GP[Species]", "GP_mu", "GP_sigma", 
                                       "G0[Species]", "G0_mu", "G0_sigma",
                                       "sigma"),
                        format = "long") %>%
  prior_posterior_plot(group_name = "Species", ridges = FALSE)
# posteriors are reasonably well constrained

# check prior-posterior of estimates of true predictor
ggsave(
  GP_samples %>%
    spread_draws(nP[n]) %>%
    mutate(nP_prior = rnorm(max(.draw), 40.5 , 20)) %>%
    pivot_longer(cols = c(nP, nP_prior),
                 values_to = "nP", names_to = "distribution") %>%
    mutate(distribution = if_else(distribution %>% str_detect("prior"),
                                  "prior", "posterior") %>% fct()) %>%
    ggplot() +
      geom_density(aes(nP, alpha = distribution),
                   fill = "black", colour = NA) +
      scale_alpha_manual(values = c(0.6, 0.2)) +
      facet_wrap(~ n, scales = "free") +
      theme_minimal() +
      theme(panel.grid = element_blank()),
  filename = "GP_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)

# 3.7 Predictions ####
# Simulate smooth prior distributions using R
GP_prior <- tibble(
  GP_mu = rtruncnorm( 8 * 1e4 , mean = 0.24 , sd = 0.2 , a = 0 ), 
  GP_sigma = rexp( 8 * 1e4 , 1 ),
  GP = rtruncnorm( 8 * 1e4 , mean = GP_mu , sd = GP_sigma , a = 0 ),
  G0_mu = rnorm( 8 * 1e4 , 0 , 2 ),
  G0_sigma = rexp( 8 * 1e4 , 1 ),
  G0 = rnorm( 8 * 1e4 , mean = G0_mu , sd = G0_sigma ),
  sigma = rexp( 8 * 1e4 , 1 )
)

GP_prior %>%
  ggplot() +
    geom_density(aes(GP_mu), colour = "black") +
    geom_density(aes(GP), colour = "grey") +
    scale_x_continuous(limits = c(0, 5), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

GP_prior %>%
  ggplot() +
    geom_density(aes(G0_mu), colour = "black") +
    geom_density(aes(G0), colour = "grey") +
    scale_x_continuous(limits = c(-10, 10), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# Simulate smooth posteriors with hyperparameters
GP_hyperposterior <- GP_samples %>%
  spread_draws(GP_mu, GP_sigma, G0_mu, G0_sigma, sigma) %>%
  mutate(GP = rtruncnorm( n() , mean = GP_mu , sd = GP_sigma , a = 0 ),
         G0 = rnorm( n(), mean = G0_mu , sd = G0_sigma )) # simulate for unobserved corallines

GP_hyperposterior %>%
  ggplot() +
    geom_density(aes(GP_mu), colour = "black") +
    geom_density(aes(GP), colour = "grey") +
    scale_x_continuous(limits = c(0, 1), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

GP_hyperposterior %>%
  ggplot() +
    geom_density(aes(G0_mu), colour = "black") +
    geom_density(aes(G0), colour = "grey") +
    scale_x_continuous(limits = c(-4, 6), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# Join priors and posteriors
GP_prior_posterior <- GP_samples %>%
  recover_types(GP_data %>% select(Species)) %>%
  spread_draws(GP[Species], G0[Species], sigma) %>%
  ungroup() %>%
  bind_rows(
    GP_hyperposterior %>% 
      mutate(Species = "Branching corallines" %>% fct()) %>%
      select(-c(GP_mu, GP_sigma, G0_mu, G0_sigma)),
    GP_prior %>%
      mutate(.chain = 1:8 %>% rep(each = 1e4),
             .iteration = 1:1e4 %>% rep(8),
             .draw = 1:(8*1e4),
             Species = "Prior" %>% fct()) %>%
      select(-c(GP_mu, GP_sigma, G0_mu, G0_sigma))
  )

str(GP_prior_posterior)

# Calculate probability of G0 > 0
G0_P <- GP_prior_posterior %>%
  group_by(Species) %>%
  summarise(mean = mean(G0),
            median = median(G0),
            mode = mode_qi(G0)$y,
            sd = sd(G0),
            P = mean(G0 > 0),
            n = length(G0)) %>%
  mutate(P_percentage = signif(P * 100, digits = 2) %>%
           str_c("%"))

# Plot GP
require(ggdist)
Fig_3b_1 <- GP_prior_posterior %>%
  mutate(Species = Species %>% 
           fct_rev() %>%
           fct_relevel("Prior", after = Inf)) %>%
  ggplot() +
    stat_slab(aes(y = Species, x = GP, fill = Species),
              height = 5, alpha = 0.5, n = 2e3) + # added samples to get smooth curvature of dist
    scale_fill_manual(values = c("#363538", "#6b4d8d", "#f5a54a", "#ec7faa", "#b5b8ba"),
                      guide = guide_legend(reverse = TRUE)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25),
                       labels = scales::label_number(accuracy = c(1, 0.01, 0.1, 0.01, 1))) +
    labs(x = expression("G:P (µmol CaCO"[3]*" µmol"^-1*" CO"[2]*")")) +
    coord_cartesian(xlim = c(0, 1), ylim = c(1, 6.5),
                    expand = FALSE, clip = "on") +
    mytheme +
    theme(legend.text = element_text(face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank())

# Plot G0
Fig_3b_2 <- GP_prior_posterior %>%
  mutate(Species = Species %>% 
           fct_rev() %>%
           fct_relevel("Prior", after = Inf)) %>%
  ggplot() +
    stat_slab(aes(y = Species, x = G0, fill = Species),
              height = 5, alpha = 0.5) +
    geom_vline(xintercept = 0) +
    geom_text(data = G0_P, aes(x = 10, y = Species, label = P_percentage),
              hjust = 1, vjust = -0.1, size = 3.5, family = "Futura") +
    scale_fill_manual(values = c("#363538", "#6b4d8d", "#f5a54a", "#ec7faa", "#b5b8ba"),
                      guide = guide_legend(reverse = TRUE)) +
    scale_x_continuous(labels = scales::label_number(style_negative = "minus")) +
    labs(x = expression("G"[0]*" (µmol CaCO"[3]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(-5, 10), ylim = c(1, 8.5),
                    expand = FALSE, clip = "on") +
    mytheme +
    theme(legend.text = element_text(face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank())

# Spread across predictor
GP_predictions <- GP_prior_posterior %>%
  left_join(GP_data %>%
              group_by(Species) %>%
              summarise(min = min(nP_mean),
                        max = max(nP_mean)),
            by = "Species") %>%
  mutate(min = if_else(is.na(min),
                       GP_data %$% min(nP_mean),
                       min),
         max = if_else(is.na(max),
                       GP_data %$% max(nP_mean),
                       max)) %>%
  rowwise() %>%
  mutate(nP_mean = list( seq(min, max, length.out = 100) )) %>%
  select(-c(min, max)) %>%
  unnest(nP_mean) %>%
  mutate(mu = G0 + GP * nP_mean,
         obs = rnorm( n() , mu , sigma ))

GP_predictions_summary <- GP_predictions %>%
  group_by(Species, nP_mean) %>%
  reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
          obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
  unnest(c(mu, obs), names_sep = "_")

require(geomtextpath)
Fig_3a <- GP_predictions_summary %>%
  ggplot() + # manually limit 1:1 line to 0-40 range
    geom_textline(data = tibble(x = c(0, 40), y = c(0, 40)), aes(x, y), 
                  label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
    geom_hdr(data = C_estimates %>% 
               drop_na(nP, G) %>% 
               group_by(ID) %>% # randomly reshuffle within ID
               mutate(G = G %>% sample()),
             aes(nP, G, fill = Species, group = ID),
             n = 600, method = "mvnorm", probs = 0.999) +
    # geom_ribbon(data = . %>% filter(Species == "Branching corallines" & mu_.width == 0.9),
    #             aes(nP_mean, ymin = mu_ymin, ymax = mu_ymax), 
    #             colour = NA, fill = "#363538", alpha = 0.3) +
    geom_line(data = . %>% filter(!Species %in% c("Prior", "Branching corallines")),
              aes(nP_mean, mu_y, colour = Species)) +
    geom_ribbon(data = . %>% filter(!Species %in% c("Prior", "Branching corallines")),
                aes(nP_mean, ymin = mu_ymin, ymax = mu_ymax,
                    fill = Species, alpha = factor(mu_.width))) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"), guide = "none") +
    scale_colour_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.5, 0.4, 0.3), guide = "none") + # geom_hdr requires additional alpha level
    scale_x_continuous(breaks = seq(0, 90, 30)) +
    labs(x = expression("P"["n"]*" (µmol CO"[2]*" g"^-1*" h"^-1*")"),
         y = expression("G (µmol CaCO"[3]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(0, 90), ylim = c(0, 40),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

require(patchwork)
Fig_3 <- ( Fig_3a | Fig_3b_1 / Fig_3b_2 ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

Fig_3 %>%
  ggsave(filename = "Fig_3mod.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# 4. Dissolution:respiration ####
# 4.1 Visualise data ####
Fig_4 <- C_estimates %>%
  drop_na(R, D) %>%
  group_by(ID) %>% # randomly reshuffle within ID
  mutate(D = D %>% sample()) %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_hdr(aes(R, D, fill = Species, group = ID),
             alpha = 0.5, n = 600, method = "mvnorm", probs = 0.999) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d")) +
    labs(x = expression("R (µmol CO"[2]*" g"^-1*" h"^-1*")"),
         y = expression("D (µmol CaCO"[3]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(-5, 25), ylim = c(-10, 15),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

Fig_4 %>%
  ggsave(filename = "Fig_4.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# 4.2 Prepare data ####
DR_data <- C_estimates_summary %>%
  drop_na(R_mean, D_mean) %>%
  select(Species, R_mean, R_sd, D_mean, D_sd)

# 4.3 Prior simulation ####
# In the regression of dissolution on dark CO2 production, an intercept is imaginable
# (net calcification and dissolution are both possible with zero CO2 production)
# but is expected to be zero (no dissolution without CO2 production) and the slope is 
# the ratio of dissolution to respiration (DR), which has to be positive 
# (truncated prior) and probably is similar to GP. From the previous model we know
# that branching corallines produce 0.32 CaCO3 for every CO2 taken up.

tibble(n = 1:1e3, # simulate hierachical prior
       DR_mu = rtruncnorm( 1e3 , mean = 0.32 , sd = 0.2 , a = 0 ), # a is lower bound
       DR_sigma = rexp( 1e3 , 1 ),
       DR = rtruncnorm( 1e3 , mean = DR_mu , sd = DR_sigma , a = 0 ),
       D0_mu = rnorm( 1e3 , 0 , 0.5 ), # dissolution is smaller than calcification, hence less variability
       D0_sigma = rexp( 1e3 , 1 ),
       D0 = rnorm( 1e3 , mean = D0_mu , sd = D0_sigma )) %>%
  expand_grid(R = DR_data %$% seq(min(R_mean), max(R_mean))) %>%
  mutate(D = D0 + DR * R) %>%
  ggplot(aes(R, D, group = n)) +
    geom_hline(yintercept = DR_data %$% c(min(D_mean), 0, max(D_mean))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(ylim = c(-10, 20), expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# I also need a prior for estimates of the true predictor
DR_data %$% ( max(R_mean) + min(R_mean) ) %>% round / 2 # I'll take half the predictor range
ggplot() + # and choose a reasonable prior centred on it
  geom_density(aes(rnorm(n = 1e5, mean = 9, sd = 6))) +
  geom_vline(xintercept = DR_data %$% c(min(R_mean), max(R_mean))) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# 4.4 Run model ####
DR_stan <- "
data{
  int n;
  vector[n] R_mean;
  vector<lower=0>[n] R_sd;
  vector[n] D_mean;
  vector<lower=0>[n] D_sd;
  array[n] int Species;
  int n_Species;
}

parameters{
  // True estimates for measurement error
  vector[n] R;
  vector[n] D;

  // Hyperparameters
  real<lower=0> DR_mu;
  real<lower=0> DR_sigma;
  real D0_mu;
  real<lower=0> D0_sigma;
  
  // Species-specific parameters
  vector<lower=0>[n_Species] DR;
  vector[n_Species] D0;
  
  // Likelihood uncertainty
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  DR_mu ~ normal( 0.32 , 0.2 ) T[0,];
  DR_sigma ~ exponential( 1 );
  D0_mu ~ normal( 0 , 1 );
  D0_sigma ~ exponential( 1 );
  
  // Species-specific priors
  DR ~ normal( DR_mu , DR_sigma ) T[0,];
  D0 ~ normal( D0_mu , D0_sigma );
  
  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );
  
  // Estimated predictor
  R ~ normal( 9 , 6 );
  R_mean ~ normal( R , R_sd );
  
  // Model
  vector[n] mu;
  for ( i in 1:n ) {
      mu[i] = D0[Species[i]] + DR[Species[i]] * R[i];
  }

  // Likelihood with measurement error
  D ~ normal( mu , sigma );
  D_mean ~ normal( D , D_sd );
}
"

DR_mod <- DR_stan %>%
  write_stan_file() %>%
  cmdstan_model()

DR_samples <- DR_mod$sample(
  data = DR_data %>% compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4,
  adapt_delta = 0.999,
  max_treedepth = 15) 
# increased adapt_delta and max_treedepth to make sampler go slower and reduce divergences

# 4.5 Model checks ####
DR_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# <5% rhat above 1.001
# good effective sample size

DR_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "DR_rank.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

# 4.6 Prior-posterior comparison ####
# sample priors
DR_prior <- prior_samples(
  model = DR_mod,
  data = DR_data %>% compose_data(),
  chains = 8, samples = 1e4)

# plot prior-posterior comparison for main parameters
DR_prior %>%
  prior_posterior_draws(posterior_samples = DR_samples,
                        group = DR_data %>% select(Species),
                        parameters = c("DR[Species]", "DR_mu", "DR_sigma", 
                                       "D0[Species]", "D0_mu", "D0_sigma",
                                       "sigma"),
                        format = "long") %>%
  prior_posterior_plot(group_name = "Species", ridges = FALSE)
# posteriors are reasonably well constrained

# check prior-posterior of estimates of true predictor
ggsave(
  DR_samples %>%
    spread_draws(R[n]) %>%
    mutate(R_prior = rnorm(max(.draw), 9 , 6)) %>%
    pivot_longer(cols = c(R, R_prior),
                 values_to = "R", names_to = "distribution") %>%
    mutate(distribution = if_else(distribution %>% str_detect("prior"),
                                  "prior", "posterior") %>% fct()) %>%
    ggplot() +
      geom_density(aes(R, alpha = distribution),
                   fill = "black", colour = NA) +
      scale_alpha_manual(values = c(0.6, 0.2)) +
      facet_wrap(~ n, scales = "free") +
      theme_minimal() +
      theme(panel.grid = element_blank()),
  filename = "DR_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)

# 4.7 Predictions ####
# Simulate smooth prior distributions using R
DR_prior <- tibble(
  DR_mu = rtruncnorm( 8 * 1e4 , mean = 0.32 , sd = 0.2 , a = 0 ), 
  DR_sigma = rexp( 8 * 1e4 , 1 ),
  DR = rtruncnorm( 8 * 1e4 , mean = DR_mu , sd = DR_sigma , a = 0 ),
  D0_mu = rnorm( 8 * 1e4 , 0 , 1 ),
  D0_sigma = rexp( 8 * 1e4 , 1 ),
  D0 = rnorm( 8 * 1e4 , mean = D0_mu , sd = D0_sigma ),
  sigma = rexp( 8 * 1e4 , 1 )
)

DR_prior %>%
  ggplot() +
    geom_density(aes(DR_mu), colour = "black") +
    geom_density(aes(DR), colour = "grey") +
    scale_x_continuous(limits = c(0, 3), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

DR_prior %>%
  ggplot() +
    geom_density(aes(D0_mu), colour = "black") +
    geom_density(aes(D0), colour = "grey") +
    scale_x_continuous(limits = c(-5, 5), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# Simulate smooth posteriors with hyperparameters
DR_hyperposterior <- DR_samples %>%
  spread_draws(DR_mu, DR_sigma, D0_mu, D0_sigma, sigma) %>%
  mutate(DR = rtruncnorm( n() , mean = DR_mu , sd = DR_sigma , a = 0 ),
         D0 = rnorm( n(), mean = D0_mu , sd = D0_sigma )) # simulate for unobserved corallines

DR_hyperposterior %>%
  ggplot() +
    geom_density(aes(DR_mu), colour = "black") +
    geom_density(aes(DR), colour = "grey") +
    scale_x_continuous(limits = c(0, 1), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

DR_hyperposterior %>%
  ggplot() +
    geom_density(aes(D0_mu), colour = "black") +
    geom_density(aes(D0), colour = "grey") +
    scale_x_continuous(limits = c(-3, 2), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# Join priors and posteriors
DR_prior_posterior <- DR_samples %>%
  recover_types(DR_data %>% select(Species)) %>%
  spread_draws(DR[Species], D0[Species], sigma) %>%
  ungroup() %>%
  bind_rows(
    DR_hyperposterior %>% 
      mutate(Species = "Branching corallines" %>% fct()) %>%
      select(-c(DR_mu, DR_sigma, D0_mu, D0_sigma)),
    DR_prior %>%
      mutate(.chain = 1:8 %>% rep(each = 1e4),
             .iteration = 1:1e4 %>% rep(8),
             .draw = 1:(8*1e4),
             Species = "Prior" %>% fct()) %>%
      select(-c(DR_mu, DR_sigma, D0_mu, D0_sigma))
  )

str(DR_prior_posterior)

# Calculate probability of D0 < 0
D0_P <- DR_prior_posterior %>%
  group_by(Species) %>%
  summarise(mean = mean(D0),
            median = median(D0),
            mode = mode_qi(D0)$y,
            sd = sd(D0),
            P = mean(D0 < 0),
            n = length(D0)) %>%
  mutate(P_percentage = signif(P * 100, digits = 2) %>%
           str_c("%"))

# Plot DR
Fig_4b_1 <- DR_prior_posterior %>%
  mutate(Species = Species %>% 
           fct_rev() %>%
           fct_relevel("Prior", after = Inf)) %>%
  ggplot() +
    stat_slab(aes(y = Species, x = DR, fill = Species),
              height = 5, alpha = 0.5, n = 2e3) + # added samples to get smooth curvature of dist
    scale_fill_manual(values = c("#363538", "#6b4d8d", "#f5a54a", "#ec7faa", "#b5b8ba"),
                      guide = guide_legend(reverse = TRUE)) +
    scale_x_continuous(breaks = seq(0, 1.5, 0.5),
                       labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1))) +
    labs(x = expression("D:R (µmol CaCO"[3]*" µmol"^-1*" CO"[2]*")")) +
    coord_cartesian(xlim = c(0, 1.5), ylim = c(1, 7.5),
                    expand = FALSE, clip = "on") +
    mytheme +
    theme(legend.text = element_text(face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank())

# Plot D0
Fig_4b_2 <- DR_prior_posterior %>%
  mutate(Species = Species %>% 
           fct_rev() %>%
           fct_relevel("Prior", after = Inf)) %>%
  ggplot() +
    stat_slab(aes(y = Species, x = D0, fill = Species),
              height = 5, alpha = 0.5, n = 1e3) +
    geom_vline(xintercept = 0) +
    geom_text(data = D0_P, aes(x = 4, y = Species, label = P_percentage),
              hjust = 1, vjust = -0.1, size = 3.5, family = "Futura") +
    scale_fill_manual(values = c("#363538", "#6b4d8d", "#f5a54a", "#ec7faa", "#b5b8ba"),
                      guide = guide_legend(reverse = TRUE)) +
    scale_x_continuous(labels = scales::label_number(style_negative = "minus")) +
    labs(x = expression("D"[0]*" (µmol CaCO"[3]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(-4, 4), ylim = c(1, 8.5),
                    expand = FALSE, clip = "on") +
    mytheme +
    theme(legend.text = element_text(face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank())

# Spread across predictor
DR_predictions <- DR_prior_posterior %>%
  left_join(DR_data %>%
              group_by(Species) %>%
              summarise(min = min(R_mean),
                        max = max(R_mean)),
            by = "Species") %>%
  mutate(min = if_else(is.na(min),
                       DR_data %$% min(R_mean),
                       min),
         max = if_else(is.na(max),
                       DR_data %$% max(R_mean),
                       max)) %>%
  rowwise() %>%
  mutate(R_mean = list( seq(min, max, length.out = 100) )) %>%
  select(-c(min, max)) %>%
  unnest(R_mean) %>%
  mutate(mu = D0 + DR * R_mean,
         obs = rnorm( n() , mu , sigma ))

DR_predictions_summary <- DR_predictions %>%
  group_by(Species, R_mean) %>%
  reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
          obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
  unnest(c(mu, obs), names_sep = "_")

Fig_4a <- DR_predictions_summary %>%
  ggplot() + # manually limit 1:1 line to -10-15 range
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_textline(data = tibble(x = c(-5, 15), y = c(-5, 15)), aes(x, y),
                  label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
    geom_hdr(data = C_estimates %>% 
               drop_na(R, D) %>% 
               group_by(ID) %>% # randomly reshuffle within ID
               mutate(D = D %>% sample()),
             aes(R, D, fill = Species, group = ID),
             n = 600, method = "mvnorm", probs = 0.999) +
    # geom_ribbon(data = . %>% filter(Species == "Branching corallines" & mu_.width == 0.9),
    #             aes(R_mean, ymin = mu_ymin, ymax = mu_ymax), 
    #             colour = NA, fill = "#363538", alpha = 0.3) +
    geom_line(data = . %>% filter(!Species %in% c("Prior", "Branching corallines")),
              aes(R_mean, mu_y, colour = Species)) +
    geom_ribbon(data = . %>% filter(!Species %in% c("Prior", "Branching corallines")),
                aes(R_mean, ymin = mu_ymin, ymax = mu_ymax,
                    fill = Species, alpha = factor(mu_.width))) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"), guide = "none") +
    scale_colour_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.5, 0.4, 0.3), guide = "none") + # geom_hdr requires additional alpha level
    labs(x = expression("R (µmol CO"[2]*" g"^-1*" h"^-1*")"),
         y = expression("D (µmol CaCO"[3]*" g"^-1*" h"^-1*")")) +
    coord_cartesian(xlim = c(-5, 20), ylim = c(-10, 15),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

Fig_4 <- ( Fig_4a | Fig_4b_1 / Fig_4b_2 ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

Fig_4 %>%
  ggsave(filename = "Fig_4mod.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)