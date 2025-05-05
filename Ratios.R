# 1. Load data ####
require(tidyverse)
C_estimates <- read_rds("C_estimates.rds") %>%
  mutate(G_P = G / nP) # PIC:POC metabolic ratio
    
# Load mass data
mass <- read.csv("Mass.csv") %>%
  mutate(ID = ID %>% fct(),
         Species = Species %>% fct(),
         Individual = Individual %>% fct(),
         # subtract empty tube mass from tube with organics
         POM = POM_Tube - Tube,
         # dry mass lost during acidification must be CaCO3
         CaCO3 = DMi - POM, # 2HCl + CaCO3 -> CO2 (gas) + CaCl2 (washed out)
         # g C = g CaCO3 / g mol^-1 CaCO3 * g mol^-1 C
         PIC = CaCO3 / 100.0869 * 12.0107,
         POC = POM * 0.3, # assuming 30% carbon content
         PIC_POC = PIC / POC) # PIC:POC mass ratio

# Join data
require(magrittr)
C_estimates %<>%
  left_join(mass %>% 
              select(ID, Species, Individual, PIC_POC),
            by = c("ID", "Species", "Individual"))
  
# Check merge
C_estimates %>%
  group_by(Round, ID, Species, Individual) %>%
  summarise(
    across(
      .cols = c(G_P, PIC_POC),
      .fns = list(mean = mean, sd = sd)
    )
  ) %>%
  ungroup() %>%
  keep(~ !all(.x == 0 | is.na(.x))) %>%
  print(n = 51)

# 2. Visualise data ####
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

require(geomtextpath)
Fig_6 <- # data cannot be specified here because geom_textvline acts up
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_textvline(xintercept = 1, label = "1:1", family = "Futura", 
                   size = 3.5, hjust = 1) +
    geom_violin(data = C_estimates %>% drop_na(PIC_POC, G_P), 
                aes(PIC_POC, G_P, fill = Species, colour = Species, group = ID),
                alpha = 0.5, position = "identity", width = 1.2) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d")) +
    scale_colour_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d")) +
    scale_y_continuous(breaks = seq(-0.5, 1, 0.5),
                       labels = scales::label_number(accuracy = c(0.1, 1, 0.1, 1),
                                                     style_negative = "minus")) +
    labs(x = "PIC:POC",
         y = expression("G:P"["n"]*" (µmol CaCO"[3]*" µmol"^-1*" CO"[2]*")")) +
    coord_cartesian(xlim = c(0, 40), ylim = c(-0.5, 1),
                    expand = FALSE) +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

Fig_6 %>%
  ggsave(filename = "Fig_6.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# 3. PIC:POC ####
# 3.1 Prepare data ####
# The first model will just look at PIC:POC. There was no good spread in mass because
# the acidification volume was limited to 20 mL, so we stuck with ~1 g dry mass. That
# means a regression would not yield a good estimate. Instead, I will run an intercept
# model on the ratio. The data need no further preparation.

# 3.2 Prior simulation ####
# Since PIC:POC is a ratio it has to be positive (gamma likelihood). Data from Haberman & Martone 
# (2023, doi: 10.3354/meps14341) for four corallines suggest a mean PIC:POC of 1.993545, or simply 2.

# simulate prior
tibble(PICPOC_mu_log = rnorm( 1e5 , mean = log(2) , sd = 0.8 ),
       PICPOC_mu = exp(PICPOC_mu_log)) %>%
  ggplot(aes(PICPOC_mu)) +
    geom_vline(xintercept = mass %$% PIC_POC, alpha = 0.05) +
    geom_density() +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# 3.3 Run model ####
PICPOC_stan <- "
data{
  int n;
  vector[n] PIC_POC;
  array[n] int Species;
  int n_Species;
}

parameters{
  // Hyperparameters
  real PICPOC_mu;
  real<lower=0> PICPOC_sigma;

  // Species-specific parameters
  vector[n_Species] PICPOC_z;

  // Likelihood uncertainty
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  PICPOC_mu ~ normal( log(2) , 0.8 );
  PICPOC_sigma ~ exponential( 1 );

  // Species-specific z-score priors
  PICPOC_z ~ normal( 0 , 1 );

  // Calculate PICPOC
  vector[n_Species] PICPOC;
  PICPOC = PICPOC_z * PICPOC_sigma + PICPOC_mu;

  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
      mu[i] = exp( PICPOC[Species[i]] );
  }

  // Likelihood reparameterised with mean and sd
  PIC_POC ~ gamma( square(mu) / square(sigma) , mu / square(sigma) );
}

generated quantities{
  // Save PICPOC
  vector[n_Species] PICPOC;
  PICPOC = PICPOC_z * PICPOC_sigma + PICPOC_mu;
}
"

require(cmdstanr)
PICPOC_mod <- PICPOC_stan %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
PICPOC_samples <- PICPOC_mod$sample(
  data = mass %>%
    select(Species, PIC_POC) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4,
  adapt_delta = 0.99)

# 3.4 Model checks ####
PICPOC_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

require(bayesplot)
PICPOC_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# chains look good

# 3.5 Prior-posterior comparison ####
source("functions.R")
# sample priors
PICPOC_prior <- prior_samples(
  model = PICPOC_mod,
  data = mass %>%
    select(Species, PIC_POC) %>%
    compose_data(),
  chains = 8, samples = 1e4)

# plot prior-posterior comparison for main parameters
PICPOC_prior %>%
  prior_posterior_draws(posterior_samples = PICPOC_samples,
                        group = mass %>% select(Species),
                        parameters = c("PICPOC[Species]", "PICPOC_mu",
                                       "PICPOC_sigma", "sigma"),
                        format = "long") %>%
  prior_posterior_plot(group_name = "Species", ridges = FALSE)
# posteriors are well constrained

# 4. Ratio comparison ####
# 4.1 Prepare data ####
# The model needs to be conditioned on a data summary and standard deviations
# are passed to it to account for measurement error. In the regression of metabolic 
# on mass ratio, a traditional intercept (x = 0) is meaningless since ratios can never 
# be zero. It is best to centre the data so that the intercept becomes an estimate of 
# the metabolic ratio at the mean mass ratio, we know to be roughly 0.3 from previous 
# models (see Calcification.R). 
R_data <- C_estimates %>%
  drop_na(PIC_POC, G_P) %>%
  group_by(ID, Species) %>%
  summarise(PIC_POC = unique(PIC_POC),
            G_P_mean = mean(G_P),
            G_P_sd = sd(G_P)) %>%
  ungroup() %>%
  mutate(PIC_POC_c = PIC_POC - mean(PIC_POC)) %>%
  select(-ID)

R_data

# 4.2 Prior simulation ####
# The slope is the relationship ratios, which is widely expected to be 1:1, but which
# can plausibly take negative values and is in reality completely unknown.
tibble(n = 1:1e3, # simulate hyperpriors
       alpha_mu = rnorm( 1e3 , mean = 0.3 , sd = 0.15 ),
       beta_mu = rnorm( 1e3 , mean = 0 , sd = 0.01 )) %>%
  expand_grid(PIC_POC_c = R_data %$% seq(min(PIC_POC_c), max(PIC_POC_c))) %>%
  mutate(G_P = alpha_mu + beta_mu * PIC_POC_c) %>%
  ggplot(aes(PIC_POC_c, G_P, group = n)) +
    geom_hline(yintercept = R_data %$% c(min(G_P_mean), 0, max(G_P_mean))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(ylim = c(-0.5, 1)) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# 4.3 Run model ####
R_stan <- "
data{
  int n;
  vector[n] PIC_POC_c;
  vector[n] G_P_mean;
  vector<lower=0>[n] G_P_sd;
  array[n] int Species;
  int n_Species;
}

parameters{
  // True estimates for measurement error
  vector[n] G_P;

  // Hyperparameters
  real GPmu_mu;
  real<lower=0> GPmu_sigma;
  real beta_mu;
  real<lower=0> beta_sigma;
  
  // Species-specific parameters
  vector[n_Species] GPmu_z;
  vector[n_Species] beta_z;
  
  // Likelihood uncertainty
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  GPmu_mu ~ normal( 0.3 , 0.15 );
  GPmu_sigma ~ exponential( 10 );
  beta_mu ~ normal( 0 , 0.01 );
  beta_sigma ~ exponential( 30 );
  
  // Species-specific z-score priors
  GPmu_z ~ normal( 0 , 1 );
  beta_z ~ normal( 0 , 1 );
  
  // Calculate GPmu and beta
  vector[n_Species] GPmu;
  GPmu = GPmu_z * GPmu_sigma + GPmu_mu;
  vector[n_Species] beta;
  beta = beta_z * beta_sigma + beta_mu;
  
  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );
  
  // Model
  vector[n] mu;
  for ( i in 1:n ) {
      mu[i] = GPmu[Species[i]] + beta[Species[i]] * PIC_POC_c[i];
  }

  // Likelihood with measurement error
  G_P ~ normal( mu , sigma );
  G_P_mean ~ normal( G_P , G_P_sd );
}

generated quantities{
  // Save GPmu and beta
  vector[n_Species] GPmu;
  GPmu = GPmu_z * GPmu_sigma + GPmu_mu;
  vector[n_Species] beta;
  beta = beta_z * beta_sigma + beta_mu;
}
"

R_mod <- R_stan %>%
  write_stan_file() %>%
  cmdstan_model()

R_samples <- R_mod$sample(
  data = R_data %>% compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4,
  adapt_delta = 0.999) 
# increased adapt_delta to make sampler go slower and reduce divergences

# 4.4 Model checks ####
R_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

R_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "R_rank.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

# 4.5 Prior-posterior comparison ####
# sample priors
R_prior <- prior_samples(
  model = R_mod,
  data = R_data %>% compose_data(),
  chains = 8, samples = 1e4)
# Stan has no problem to sample joint prior distributions of hierarchical models
# when non-centred parameterisation is used.

# plot prior-posterior comparison for main parameters
R_prior %>%
  prior_posterior_draws(posterior_samples = R_samples,
                        group = R_data %>% select(Species),
                        parameters = c("GPmu[Species]", "GPmu_mu", "GPmu_sigma", 
                                       "beta[Species]", "beta_mu", "beta_sigma", 
                                       "sigma"),
                        format = "long") %>%
  prior_posterior_plot(group_name = "Species", ridges = FALSE)
# posteriors are reasonably well constrained

# 3.7 Predictions ####
R_prior_posterior <- R_samples %>%
  spread_draws(GPmu[Species], beta[Species], sigma) %>%
  ungroup() %>%
  bind_rows(
    R_prior %>%
      prior_posterior_draws(posterior_samples = R_samples,
                            group = tibble(NULL),
                            parameters = c("GPmu_mu", "beta_mu", 
                                           "GPmu_sigma", "beta_sigma", 
                                           "sigma"),
                            format = "short") %>%
      mutate(GPmu = rnorm( n() , GPmu_mu , GPmu_sigma ),
             beta = rnorm( n() , beta_mu , beta_sigma )) %>%
      mutate(Species = if_else(distribution == "prior",
                               "Prior", "Branching corallines") %>%
                          fct()) %>%
      select(-c(GPmu_mu, beta_mu, 
                GPmu_sigma, beta_sigma,
                distribution))
  )

str(R_prior_posterior)

# Calculate probability of beta < 0
beta_P <- R_prior_posterior %>%
  group_by(Species) %>%
  summarise(mean = mean(beta),
            median = median(beta),
            mode = mode_qi(beta)$y,
            sd = sd(beta),
            P = mean(beta < 0),
            n = length(beta)) %>%
  mutate(P_percentage = signif(P * 100, digits = 2) %>%
           str_c("%"))

# Plot beta
require(ggdist)
Fig_6b_1 <- R_prior_posterior %>%
  mutate(Species = Species %>% 
           fct_rev() %>%
           fct_relevel("Prior", after = Inf)) %>%
  ggplot() +
    stat_slab(aes(y = Species, x = beta, fill = Species),
              height = 5, alpha = 0.5) +
    geom_vline(xintercept = 0) +
    geom_text(data = beta_P, aes(x = 0.04, y = Species, label = P_percentage),
              hjust = 1, vjust = -0.1, size = 3.5, family = "Futura") +
    scale_fill_manual(values = c("#363538", "#6b4d8d", "#f5a54a", "#ec7faa", "#b5b8ba"),
                      guide = guide_legend(reverse = TRUE)) +
    scale_x_continuous(breaks = seq(-0.04, 0.04, 0.02),
                       labels = scales::label_number(accuracy = c(0.01, 0.01, 1, 0.01, 0.01),
                                                     style_negative = "minus")) +
    labs(x = expression(italic(beta)*" (G:P"["n"]*" PIC:POC"^-1*")")) +
    coord_cartesian(xlim = c(-0.04, 0.04), ylim = c(1, 8),
                    expand = FALSE, clip = "on") +
    mytheme +
    theme(legend.text = element_text(face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank())

# Plot GPmu
Fig_6b_2 <- R_prior_posterior %>%
  mutate(Species = Species %>% 
           fct_rev() %>%
           fct_relevel("Prior", after = Inf)) %>%
  ggplot() +
    stat_slab(aes(y = Species, x = GPmu, fill = Species),
              height = 5, alpha = 0.5) +
    scale_fill_manual(values = c("#363538", "#6b4d8d", "#f5a54a", "#ec7faa", "#b5b8ba"),
                      guide = guide_legend(reverse = TRUE)) +
    scale_x_continuous(breaks = seq(0, 0.6, 0.2),
                       labels = scales::label_number(accuracy = c(1, 0.1, 0.1, 0.1))) +
    labs(x = expression(bar("G:P"["n"])*" (µmol CaCO"[3]*" µmol"^-1*" CO"[2]*")")) +
    coord_cartesian(xlim = c(0, 0.6), ylim = c(1, 8.5),
                    expand = FALSE, clip = "on") +
    mytheme +
    theme(legend.text = element_text(face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank())

# Spread across predictor
R_predictions <- R_prior_posterior %>%
  left_join(R_data %>%
              group_by(Species) %>%
              summarise(min = min(PIC_POC_c),
                        max = max(PIC_POC_c)),
            by = "Species") %>%
  mutate(min = if_else(is.na(min),
                       R_data %$% min(PIC_POC_c),
                       min),
         max = if_else(is.na(max),
                       R_data %$% max(PIC_POC_c),
                       max)) %>%
  rowwise() %>%
  mutate(PIC_POC_c = list( seq(min, max, length.out = 100) )) %>%
  select(-c(min, max)) %>%
  unnest(PIC_POC_c) %>%
  mutate(mu = GPmu + beta * PIC_POC_c,
         obs = rnorm( n() , mu , sigma ))

R_predictions_summary <- R_predictions %>%
  group_by(Species, PIC_POC_c) %>%
  reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
          obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
  unnest(c(mu, obs), names_sep = "_") %>%
  mutate(PIC_POC = PIC_POC_c + R_data %$% mean(PIC_POC)) # undo centring

# Load posteriors from Calcification.R
GP_prior_posterior <- read_rds("GP_prior_posterior.rds")

Fig_6a <- # data cannot be specified here because geom_textvline acts up
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_violin(data = C_estimates %>% drop_na(PIC_POC, G_P), 
                aes(PIC_POC, G_P, fill = Species, colour = Species, group = ID),
                alpha = 0.5, position = "identity", width = 1.2) +
    geom_textvline(xintercept = 1, label = "1:1", family = "Futura", 
                   size = 3.5, hjust = 1) +
    # geom_ribbon(data = . %>% filter(Species == "Branching corallines" & mu_.width == 0.9),
    #             aes(PIC_POC, ymin = mu_ymin, ymax = mu_ymax), 
    #             colour = NA, fill = "#363538", alpha = 0.3) +
    geom_line(data = R_predictions_summary %>% 
                filter(!Species %in% c("Prior", "Branching corallines")),
              aes(PIC_POC, mu_y, colour = Species)) +
    geom_ribbon(data = R_predictions_summary %>% 
                  filter(!Species %in% c("Prior", "Branching corallines")),
                aes(PIC_POC, ymin = mu_ymin, ymax = mu_ymax,
                    fill = Species, alpha = factor(mu_.width))) +
    # stat_slab(data = mass, aes(x = PIC_POC, y = -0.5, fill = Species), 
    #           colour = NA, alpha = 0.5, scale = 0.15, trim = FALSE, expand = TRUE) +
    stat_slab(data = GP_prior_posterior %>%
                filter(!Species %in% c("Prior", "Branching corallines")),
              aes(y = GP, x = 0, fill = Species), 
              colour = NA, alpha = 0.5, scale = -4) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"), guide = "none") +
    scale_colour_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_y_continuous(breaks = seq(-0.4, 0.8, 0.4),
                       labels = scales::label_number(accuracy = c(0.1, 1, 0.1, 0.1),
                                                     style_negative = "minus")) +
    labs(x = "PIC:POC",
         y = expression("G:P"["n"]*" (µmol CaCO"[3]*" µmol"^-1*" CO"[2]*")")) +
    coord_cartesian(xlim = c(0, 40), ylim = c(-0.4, 0.8),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

require(patchwork)
Fig_6 <- ( Fig_6a | Fig_6b_1 / Fig_6b_2 ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

Fig_6 %>%
  ggsave(filename = "Fig_6mod.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)