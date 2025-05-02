# 1. Load data ####
require(tidyverse)
C_estimates <- read_rds("C_estimates.rds") %>%
  mutate(G_P = G / nP) # PIC:POC metabolic ratio
    
# Load mass data
mass <- read.csv("Mass.csv") %>%
  mutate(ID = ID %>% fct(),
         Species = Species %>% fct(),
         Individual = Individual %>% fct(),
         POM = POM_Tube - Tube,
         # dry mass lost during acidification must be CaCO3
         CaCO3 = DMi - POM,
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
    geom_violin(data = C_estimates, aes(PIC_POC, G_P, fill = Species, colour = Species, group = ID),
                alpha = 0.5, position = "identity", width = 1.5) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d")) +
    scale_colour_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d")) +
    scale_y_continuous(breaks = seq(-0.5, 1, 0.5),
                       labels = scales::label_number(accuracy = c(0.1, 1, 0.1, 1),
                                                     style_negative = "minus")) +
    labs(x = "PIC:POC",
         y = expression("G:P"["n"])) +
    coord_cartesian(xlim = c(0, 40), ylim = c(-0.5, 1),
                    expand = FALSE) +
    mytheme +
    theme(legend.text = element_text(face = "italic"))

Fig_6 %>%
  ggsave(filename = "Fig_6.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

###################

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

