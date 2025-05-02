# 1. Load data ####
require(tidyverse)
O2_estimates <- read_rds("O2_estimates.rds")
C_estimates <- read_rds("C_estimates.rds")

# Combine O2 and C estimates
O2_C_estimates <- C_estimates %>%
  select(starts_with("."), ID, Round, Species, Individual, G, D, nP, R, gP, 
         TA_L, TA_D, TA_LD, DIC_L, DIC_D, DIC_LD) %>%
  rename(nP_CO2 = nP, R_CO2 = R, gP_CO2 = gP) %>%
  full_join(O2_estimates %>%
              select(starts_with("."), ID, Round, Species, Individual, nP, R, gP, 
                     T_L, T_D, T_LD, P_L, P_D, P_LD) %>%
              rename(nP_O2 = nP, R_O2 = R, gP_O2 = gP),
            by = c(".chain", ".iteration", ".draw", "ID", "Round", "Species", "Individual")) 

# Add salinity data
S <- read.csv("Incubation.csv") %>%
  select(ID, Round, Species, Treatment, Salinity) %>%
  mutate(ID = ID %>% str_sub(end = -3) %>% fct(),
         Round = Round %>% as.character() %>% fct(),
         Species = Species %>% fct(),
         Treatment = Treatment %>% fct()) %>%
  pivot_wider(names_from = Treatment, values_from = Salinity) %>%
  rename(S_L = Light, S_D = Dark) %>%
  mutate(S_LD = ( S_L + S_D ) / 2)

require(magrittr)
S %<>%
  drop_na(Species) %>%
  left_join(S %>%
              filter(ID %>% str_detect("^0")) %>%
              rename(S_0_L = S_L, S_0_D = S_D, S_0_LD = S_LD) %>%
              select(-c(ID, Species)),
            by = "Round")

O2_C_estimates %<>%
  left_join(S, by = c("ID", "Round", "Species"))

str(O2_C_estimates)

# Check merge
O2_C_estimates %>%
  group_by(Round, ID, Species, Individual) %>%
  summarise(
    across(
      .cols = c(G, D, nP_CO2, R_CO2, gP_CO2, nP_O2, R_O2, gP_O2, S_L, S_D, S_LD, 
                S_0_L, S_0_D, S_0_LD, TA_L, TA_D, TA_LD, DIC_L, DIC_D, DIC_LD, 
                T_L, T_D, T_LD, P_L, P_D, P_LD),
      .fns = list(mean = mean, sd = sd)
    )
  ) %>%
  ungroup() %>%
  keep(~ !all(.x == 0 | is.na(.x))) %>%
  print(n = 51)

# 2. Calculate CO2 budget ####
# Calculate O2 budget and seawater density
O2_C_estimates %<>%
  mutate(CDR_naive = ( nP_O2 - R_O2 ) * 12, # assuming PQ = 1 and G = 0
         kg_L = rho(S = S_0_LD, T = T_LD) %>% as.numeric() * 1e-3) # convert g L^-1 to kg L^-1

# Calculate psi
require(seacarb)
# seacarb::psi() cannot handle NAs for TA, DIC, S and T, so it is important to drop
# NAs beforehand. There are no NAs in S and TA and DIC have the same missing rows so
O2_C_estimates %<>% drop_na(TA_LD, T_LD)

# seacarb::psi() also takes very long for computation, longer than seacarb::carb(),
# so it is necessary to split the tibble into chunks and run them separately
O2_C_estimates_split <- O2_C_estimates %>%
  mutate(chunk = 1:parallel::detectCores() %>% 
           rep(each = n() / parallel::detectCores()) %>% 
           as.character() %>% 
           fct()) %>%
  group_split(chunk)

require(furrr)
plan(multisession)

O2_C_estimates_split %<>%
  future_map(~ .x %>%
               mutate(psi = psi(flag = 15,
                                var1 = TA_LD * 1e-6 / kg_L, # convert µM to mol kg^-1
                                var2 = DIC_LD * 1e-6 / kg_L,
                                S = S_0_LD,
                                T = T_LD,
                                Patm = P_LD / 1013.25) %>% # convert hPa to atm
                              as.numeric())
             )

O2_C_estimates <- O2_C_estimates_split %>%
  map(~ .x %>% select(-chunk)) %>%
  bind_rows()

# Calculate CO2 budget
O2_C_estimates %<>%
  mutate(CDR_true = ( nP_CO2 - R_CO2 ) * 12 - ( G - D ) * 12 * psi)
    

# 3. Predict CO2 budget ####
# 3.1 Visualise data ####
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

require(ggdensity)
Fig_5 <- O2_C_estimates %>%
  ggplot() +
    geom_hdr(aes(CDR_naive/1e3, CDR_true/1e3, fill = Species, group = ID),
             alpha = 0.5, n = 600, method = "mvnorm", probs = 0.999) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d")) +
    scale_x_continuous(breaks = seq(0, 2, 0.5),
                       labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1))) +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       labels = scales::label_number(accuracy = c(1, 0.01, 0.1, 0.01, 1))) +
    labs(x = expression("CDR"["naive"]*" (mmol CO"[2]*" g"^-1*" d"^-1*")"),
         y = expression("CDR"["true"]*" (mmol CO"[2]*" g"^-1*" d"^-1*")")) +
    coord_cartesian(xlim = c(0, 2), ylim = c(0, 1),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

Fig_5 %>%
  ggsave(filename = "Fig_5.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# 3.2 Prepare data ####
# The model needs to be conditioned on a data summary and standard deviations
# are passed to it to account for measurement error.
CDR_data <- O2_C_estimates %>%
  group_by(ID, Species) %>%
  summarise(
    across(
      .cols = c(CDR_naive, CDR_true),
      .fns = list(mean = mean, sd = sd)
    )
  ) %>%
  ungroup() %>%
  select(-ID)

# 3.3 Prior simulation ####
# In the regression of true on naive CO2 fixation, an intercept is imaginable since
# dissolution may ameliorate the CO2 efflux caused by respiration, for example
# causing true CDR to be positive when naive CDR is zero.

# The slope is the ratio of true CDR to naive CDR, which represents the proportion of 
# naive CDR that is actually removed and as a proportion has to be positive and smaller than 1
# (1:1 relationship is impossible given that PQ is factored in and G is subtracted).
# The beta distribution is ideal for proportions but is hard to use as a hierarchical prior.
# The normal distribution truncated to the 0-1 space will do the job.


require(truncnorm) # simulate from the truncated normal
tibble(n = 1:1e3, # simulate hierachical prior
       p_mu = rtruncnorm( 1e3 , mean = 0.5 , sd = 0.2 , a = 0 , b = 1 ), # a is lower bound, b is upper bound
       p_sigma = rexp( 1e3 , 1 ),
       p = rtruncnorm( 1e3 , mean = p_mu , sd = p_sigma , a = 0 , b = 1 )) %>%
  expand_grid(CDR_naive = CDR_data %$% seq(min(CDR_naive_mean), max(CDR_naive_mean))) %>%
  mutate(CDR_true = p * CDR_naive) %>%
  ggplot(aes(CDR_naive, CDR_true, group = n)) +
    geom_abline(slope = 1) +
    geom_hline(yintercept = CDR_data %$% c(min(CDR_true_mean), 0, max(CDR_true_mean))) +
    geom_line(alpha = 0.05) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# I also need a prior for estimates of the true predictor
CDR_data %$% max(CDR_naive_mean) %>% round / 2 # I'll take half the predictor range
ggplot() + # and choose a reasonable prior centred on it
  geom_density(aes(rnorm(n = 1e5, mean = 860, sd = 400))) +
  geom_vline(xintercept = CDR_data %$% c(min(CDR_naive_mean), max(CDR_naive_mean))) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# 3.4 Run model ####
CDR_stan <- "
data{
  int n;
  vector[n] CDR_naive_mean;
  vector<lower=0>[n] CDR_naive_sd;
  vector[n] CDR_true_mean;
  vector<lower=0>[n] CDR_true_sd;
  array[n] int Species;
  int n_Species;
}

parameters{
  // True estimates for measurement error
  vector[n] CDR_naive;
  vector[n] CDR_true;

  // Hyperparameters
  real<lower=0, upper=1> p_mu;
  real<lower=0> p_sigma;
  
  // Species-specific parameters
  vector<lower=0, upper=1>[n_Species] p;
  
  // Likelihood uncertainty
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  p_mu ~ normal( 0.5 , 0.2 ) T[0,1];
  p_sigma ~ exponential( 1 );
  
  // Species-specific priors
  p ~ normal( p_mu , p_sigma ) T[0,1];
  
  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );
  
  // Estimated predictor
  CDR_naive ~ normal( 860 , 400 );
  CDR_naive_mean ~ normal( CDR_naive , CDR_naive_sd );
  
  // Model
  vector[n] mu;
  for ( i in 1:n ) {
      mu[i] = p[Species[i]] * CDR_naive[i];
  }

  // Likelihood with measurement error
  CDR_true ~ normal( mu , sigma );
  CDR_true_mean ~ normal( CDR_true , CDR_true_sd );
}
"

require(cmdstanr)
CDR_mod <- CDR_stan %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
CDR_samples <- CDR_mod$sample(
  data = CDR_data %>% compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4,
  adapt_delta = 0.999,
  max_treedepth = 15) 
# increased adapt_delta and max_treedepth to make sampler go slower and reduce divergences

# 3.5 Model checks ####
CDR_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

require(bayesplot)
CDR_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "CDR_rank.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# chains look good

# 3.6 Prior-posterior comparison ####
source("functions.R")
# sample priors
CDR_prior <- prior_samples(
  model = CDR_mod,
  data = CDR_data %>% compose_data(),
  chains = 8, samples = 1e4)
# Stan struggles to sample joint prior distributions of hierarchical models,
# but it's good enough for a quick prior-posterior check.

# plot prior-posterior comparison for main parameters
CDR_prior %>%
  prior_posterior_draws(posterior_samples = CDR_samples,
                        group = CDR_data %>% select(Species),
                        parameters = c("p[Species]", "p_mu", "p_sigma", 
                                       "sigma"),
                        format = "long") %>%
  prior_posterior_plot(group_name = "Species", ridges = FALSE)
# posteriors are reasonably well constrained

# check prior-posterior of estimates of true predictor
ggsave(
  CDR_samples %>%
    spread_draws(CDR_naive[n]) %>%
    mutate(CDR_naive_prior = rnorm(max(.draw), 860 , 400)) %>%
    pivot_longer(cols = c(CDR_naive, CDR_naive_prior),
                 values_to = "CDR_naive", names_to = "distribution") %>%
    mutate(distribution = if_else(distribution %>% str_detect("prior"),
                                  "prior", "posterior") %>% fct()) %>%
    ggplot() +
      geom_density(aes(CDR_naive, alpha = distribution),
                   fill = "black", colour = NA) +
      scale_alpha_manual(values = c(0.6, 0.2)) +
      facet_wrap(~ n, scales = "free") +
      theme_minimal() +
      theme(panel.grid = element_blank()),
  filename = "CDR_prior_posterior.pdf", path = "Plots",
  width = 80, height = 40, unit = "cm", device = cairo_pdf)

# 3.7 Predictions ####
# Simulate smooth prior distributions using R
CDR_prior <- tibble(
  p_mu = rtruncnorm( 8 * 1e4 , mean = 0.5 , sd = 0.2 , a = 0 , b = 1 ), 
  p_sigma = rexp( 8 * 1e4 , 1 ),
  p = rtruncnorm( 8 * 1e4 , mean = p_mu , sd = p_sigma , a = 0 , b = 1 ),
  sigma = rexp( 8 * 1e4 , 1 )
)

CDR_prior %>%
  ggplot() +
    geom_density(aes(p_mu), colour = "black") +
    geom_density(aes(p), colour = "grey") +
    scale_x_continuous(limits = c(0, 1), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# Simulate smooth posteriors with hyperparameters
CDR_hyperposterior <- CDR_samples %>%
  spread_draws(p_mu, p_sigma, sigma) %>% # simulate for unobserved corallines
  mutate(p = rtruncnorm( n() , mean = p_mu , sd = p_sigma , a = 0 , b = 1 )) 

CDR_hyperposterior %>%
  ggplot() +
    geom_density(aes(p_mu), colour = "black") +
    geom_density(aes(p), colour = "grey") +
    scale_x_continuous(limits = c(0, 1), oob = scales::oob_keep) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# Join priors and posteriors
CDR_prior_posterior <- CDR_samples %>%
  recover_types(CDR_data %>% select(Species)) %>%
  spread_draws(p[Species], sigma) %>%
  ungroup() %>%
  bind_rows(
    CDR_hyperposterior %>% 
      mutate(Species = "Branching corallines" %>% fct()) %>%
      select(-c(p_mu, p_sigma)),
    CDR_prior %>%
      mutate(.chain = 1:8 %>% rep(each = 1e4),
             .iteration = 1:1e4 %>% rep(8),
             .draw = 1:(8*1e4),
             Species = "Prior" %>% fct()) %>%
      select(-c(p_mu, p_sigma))
  )

str(CDR_prior_posterior)

# Calculate percentage loss relative to naive CDR
p_loss_summary <- CDR_prior_posterior %>%
  group_by(Species) %>%
  mean_qi(mean = 1 - p, .width = 0.9) %>% # 1 - p represents the proportional loss
  mutate(percentage = str_c(
    signif(mean * 100, digits = 2), "[",
    signif(.lower * 100, digits = 2), ",",
    signif(.upper * 100, digits = 2), "]", "%"
    ))

# Plot 1 - p (complementary proportion which here is CDR loss)
require(ggdist)
Fig_5b <- CDR_prior_posterior %>%
  mutate(Species = Species %>% 
           fct_rev() %>%
           fct_relevel("Prior", after = Inf)) %>%
  ggplot() +
    stat_slab(aes(y = Species, x = ( 1 - p ) * 100, fill = Species),
              height = 5, alpha = 0.5) +
    geom_text(data = p_loss_summary, aes(x = 100, y = Species, label = percentage),
              hjust = 1, vjust = -1, size = 3.5, family = "Futura") +
    scale_fill_manual(values = c("#363538", "#6b4d8d", "#f5a54a", "#ec7faa", "#b5b8ba"),
                      guide = guide_legend(reverse = TRUE)) +
    labs(x = "CDR loss (%)") +
    coord_cartesian(xlim = c(0, 100), ylim = c(1, 7.2),
                    expand = FALSE, clip = "on") +
    mytheme +
    theme(legend.text = element_text(face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank())


# Spread across predictor
CDR_predictions <- CDR_prior_posterior %>%
  left_join(CDR_data %>%
              group_by(Species) %>%
              summarise(min = min(CDR_naive_mean),
                        max = max(CDR_naive_mean)),
            by = "Species") %>%
  mutate(min = if_else(is.na(min),
                       CDR_data %$% min(CDR_naive_mean),
                       min),
         max = if_else(is.na(max),
                       CDR_data %$% max(CDR_naive_mean),
                       max)) %>%
  rowwise() %>%
  mutate(CDR_naive_mean = list( seq(min, max, length.out = 100) )) %>%
  select(-c(min, max)) %>%
  unnest(CDR_naive_mean) %>%
  mutate(mu = p * CDR_naive_mean,
         obs = rnorm( n() , mu , sigma ))

CDR_predictions_summary <- CDR_predictions %>%
  group_by(Species, CDR_naive_mean) %>%
  reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
          obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
  unnest(c(mu, obs), names_sep = "_")

require(geomtextpath)
Fig_5a <- CDR_predictions_summary %>%
  ggplot() + # manually limit 1:1 line to 0-40 range
    geom_textline(data = tibble(x = c(0, 1), y = c(0, 1)), aes(x, y),
                  label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
    geom_hdr(data = O2_C_estimates, aes(CDR_naive/1e3, CDR_true/1e3, fill = Species, group = ID),
             alpha = 0.5, n = 600, method = "mvnorm", probs = 0.999) +
    # geom_ribbon(data = . %>% filter(Species == "Branching corallines" & mu_.width == 0.9),
    #             aes(CDR_naive_mean, ymin = mu_ymin, ymax = mu_ymax), 
    #             colour = NA, fill = "#363538", alpha = 0.3) +
    geom_line(data = . %>% filter(!Species %in% c("Prior", "Branching corallines")),
              aes(CDR_naive_mean/1e3, mu_y/1e3, colour = Species)) +
    geom_ribbon(data = . %>% filter(!Species %in% c("Prior", "Branching corallines")),
                aes(CDR_naive_mean/1e3, ymin = mu_ymin/1e3, ymax = mu_ymax/1e3,
                    fill = Species, alpha = factor(mu_.width))) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"), guide = "none") +
    scale_colour_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.5, 0.4, 0.3), guide = "none") + # geom_hdr requires additional alpha level
    scale_x_continuous(breaks = seq(0, 2, 0.5),
                       labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1))) +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       labels = scales::label_number(accuracy = c(1, 0.01, 0.1, 0.01, 1))) +
    labs(x = expression("CDR"["naive"]*" (mmol CO"[2]*" g"^-1*" d"^-1*")"),
         y = expression("CDR"["true"]*" (mmol CO"[2]*" g"^-1*" d"^-1*")")) +
    coord_cartesian(xlim = c(0, 2), ylim = c(0, 1),
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

require(patchwork)
Fig_5 <- ( Fig_5a | Fig_5b ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

Fig_5 %>%
  ggsave(filename = "Fig_5mod.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)