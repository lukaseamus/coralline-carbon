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
require(seacarb)
O2_C_estimates %<>%
  mutate(CDR_naive = ( nP_O2 - R_O2 ) * 12, # assuming PQ = 1 and G = 0
         kg_L = rho(S = S_0_LD, T = T_LD) %>% as.numeric() * 1e-3) # convert g L^-1 to kg L^-1

# Calculate psi
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
  chains = 8, samples = 1e4,
  adapt_delta = 0.999,
  max_treedepth = 15
  )
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
CDR_samples %>%
  spread_draws(CDR_naive[n]) %>%
  mutate(CDR_naive_prior = rnorm(max(.draw), 860 , 400)) %>%
  pivot_longer(cols = c(CDR_naive, CDR_naive_prior),
               values_to = "CDR_naive", names_to = "distribution") %>%
  mutate(distribution = if_else(distribution %>% str_detect("prior"),
                                "prior", "posterior") %>% fct()) %>%
  { ggplot(.) +
      geom_density(aes(CDR_naive, alpha = distribution),
                   fill = "black", colour = NA) +
      scale_alpha_manual(values = c(0.6, 0.2)) +
      facet_wrap(~ n, scales = "free") +
      theme_minimal() +
      theme(panel.grid = element_blank()) } %>% 
  ggsave(filename = "CDR_prior_posterior.pdf", path = "Plots",
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
  ggplot() +
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

# 3.8 Animation ####
# 3.8.1 Build geometries ####
# Load internal ggdensity functions from 
# https://github.com/jamesotto852/ggdensity/tree/main/R.
fix_probs <- function(probs) {
  stopifnot("Probabilities must be between 0 and 1, exclusive" = all(probs > 0) & all(probs < 1))
  
  sort(probs, decreasing = TRUE)
}

tibble0 <- function(...) {
  tibble::tibble(..., .name_repair = "minimal")
}

unique0 <- function(x, ...) {
  if (is.null(x)) x else vctrs::vec_unique(x, ...)
}

isoband_z_matrix <- function(data) {
  x_pos <- as.integer(factor(data$x, levels = sort(unique0(data$x))))
  y_pos <- as.integer(factor(data$y, levels = sort(unique0(data$y))))
  nrow <- max(y_pos)
  ncol <- max(x_pos)
  raster <- matrix(NA_real_, nrow = nrow, ncol = ncol)
  raster[cbind(y_pos, x_pos)] <- data$z
  raster
}

xyz_to_isobands <- function(data, breaks) {
  isoband::isobands(x = sort(unique0(data$x)), y = sort(unique0(data$y)),
                    z = isoband_z_matrix(data), levels_low = breaks[-length(breaks)],
                    levels_high = breaks[-1])
}

iso_to_polygon <- function(iso, group = 1) {
  lengths <- vapply(iso, function(x) length(x$x), integer(1))
  if (all(lengths == 0)) {
    warning("Zero contours were generated")
    return(tibble0())
  }
  levels <- names(iso)
  xs <- unlist(lapply(iso, "[[", "x"), use.names = FALSE)
  ys <- unlist(lapply(iso, "[[", "y"), use.names = FALSE)
  ids <- unlist(lapply(iso, "[[", "id"), use.names = FALSE)
  item_id <- rep(seq_along(iso), lengths)
  groups <- paste(group, sprintf("%03d", item_id), sep = "-")
  groups <- factor(groups)
  tibble0(level = rep(levels, lengths), x = xs, y = ys,
          piece = as.integer(groups), group = groups, subgroup = ids,
          .size = length(xs))
}

xyz_to_isolines <- function(data, breaks) {
  isoband::isolines(x = sort(unique0(data$x)), y = sort(unique0(data$y)),
                    z = isoband_z_matrix(data), levels = breaks)
}

iso_to_path <- function(iso, group = 1) {
  lengths <- vapply(iso, function(x) length(x$x), integer(1))
  if (all(lengths == 0)) {
    warning("Zero contours were generated")
    return(tibble0())
  }
  levels <- names(iso)
  xs <- unlist(lapply(iso, "[[", "x"), use.names = FALSE)
  ys <- unlist(lapply(iso, "[[", "y"), use.names = FALSE)
  ids <- unlist(lapply(iso, "[[", "id"), use.names = FALSE)
  item_id <- rep(seq_along(iso), lengths)
  groups <- paste(group, sprintf("%03d", item_id), sprintf("%03d",
                                                           ids), sep = "-")
  groups <- factor(groups)
  tibble0(level = rep(levels, lengths), x = xs, y = ys,
          piece = as.integer(groups), group = groups, .size = length(xs))
}

res_to_df <- function(res, probs, group, output) {
  
  probs <- fix_probs(probs)
  
  # Need z for xyz_to_isobands/lines()
  res$df_est$z <- res$df_est$fhat
  
  if (output == "bands") {
    
    isobands <- xyz_to_isobands(res$df_est, res$breaks)
    names(isobands) <- scales::percent_format(accuracy = 1)(probs)
    df <- iso_to_polygon(isobands, group)
    df$probs <- ordered(df$level, levels = names(isobands))
    df$level <- NULL
    
  } else if (output == "lines") {
    
    isolines <- xyz_to_isolines(res$df_est, res$breaks)
    names(isolines) <- scales::percent_format(accuracy = 1)(probs)
    df <- iso_to_path(isolines, group)
    df$probs <- ordered(df$level, levels = names(isolines))
    df$level <- NULL
    
  } else if (output == "points") {
    
    df <- res$data
    df$hdr_membership <- scales::percent_format(accuracy = 1)(df$hdr_membership)
    df$probs <- ordered(df$hdr_membership, levels = scales::percent_format(accuracy = 1)(c(1, probs)))
    df$hdr_membership <- NULL
    
  }
  
  df
  
}

# O2_C_estimates <- tibble(
#   CDR_naive = c(rnorm(1e3, 1, 0.5), rnorm(1e3, 2, 0.8), rnorm(1e3, 0.5, 0.4),
#                 rnorm(1e3, 0.4, 0.1), rnorm(1e3, 0.3, 0.2), rnorm(1e3, 3, 0.8)),
#   CDR_true = c(rnorm(1e3, 0.3, 0.2), rnorm(1e3, 2.5, 1), rnorm(1e3, 0.7, 0.25),
#                rnorm(1e3, 1.5, 0.45), rnorm(1e3, 4, 2), rnorm(1e3, 0.1, 0.1)),
#   ID = c("A1", "A2", "A3", "A4", "A5", "A6") %>% rep(each = 1e3),
#   Species = c("bug1", "bug2", "bug3") %>% rep(each = 2e3)
# ) %T>%
#   print()

ellipses <- O2_C_estimates %>%
  select(ID, Species, CDR_naive, CDR_true) %>%
  rename(x = CDR_naive, y = CDR_true) %>%
  mutate(x = x/1e3, y = y/1e3) %>% # convert µmol to mmol
  group_by(ID, Species) %>%
  nest() %>%
  mutate(
    hdr = map2(
      data, ID,
      ~ .x %>%
        get_hdr(n = 600, method = "mvnorm", probs = 0.999,
                rangex = .x %$% c(min(x) - 0.1, max(x) + 0.1), # add buffer
                rangey = .x %$% c(min(y) - 0.1, max(y) + 0.1),
                hdr_membership = FALSE) %>%
        res_to_df(probs = 0.999, group = .y, output = "bands")
    )
  ) %>%
  select(-data) %>%
  unnest(hdr) %>%
  ungroup() %>%
  mutate(
    fill = case_when(
      Species == "Amphiroa anceps" ~ "#ec7faa",
      Species == "Jania rosea" ~ "#f5a54a",
      Species == "Metamastophora flabellata" ~ "#6b4d8d"
    ),
    alpha = 0.5,
    id = group %>% 
      str_c(subgroup) %>% 
      str_split_i(pattern = "_",
                  i = 2)
  ) %>%
  select(-c(piece, .size, probs)) %T>%
  print()

lines <- CDR_predictions_summary %>%
  filter(!Species %in% c("Prior", "Branching corallines")) %>%
  droplevels() %>%
  select(Species, CDR_naive_mean, mu_y) %>%
  rename(x = CDR_naive_mean, y = mu_y) %>%
  mutate(
    x = x/1e3, y = y/1e3,
    colour = case_when(
      Species == "Amphiroa anceps" ~ "#ec7faa",
      Species == "Jania rosea" ~ "#f5a54a",
      Species == "Metamastophora flabellata" ~ "#6b4d8d"
    )
  ) %T>%
  print()

ribbons <- CDR_predictions_summary %>%
  filter(!Species %in% c("Prior", "Branching corallines")) %>%
  droplevels() %>%
  group_by(Species, mu_.width) %>%
  reframe(x = c(CDR_naive_mean, rev(CDR_naive_mean)),
          y = c(mu_ymax, rev(mu_ymin))) %>%
  mutate(
    x = x/1e3, y = y/1e3,
    alpha = case_when(
      mu_.width == 0.9 ~ 0.5, 
      mu_.width == 0.8 ~ 0.4, 
      mu_.width == 0.5 ~ 0.3
    ),
    fill = case_when(
      Species == "Amphiroa anceps" ~ "#ec7faa",
      Species == "Jania rosea" ~ "#f5a54a",
      Species == "Metamastophora flabellata" ~ "#6b4d8d"
    )
  ) %>%
  select(-mu_.width) %T>%
  print()

densities <- CDR_prior_posterior %>%
  filter(!Species %in% c("Prior", "Branching corallines")) %>%
  droplevels() %>%
  mutate(p_c = ( 1 - p ) * 100) %>% # calculate complementary proportion as percentage
  group_by(Species) %>%
  reframe(x = c(0, density(p_c, n = 2^10, from = 0, to = 100, bw = 100 * 0.01)$x, 100),
          y = c(0, density(p_c, n = 2^10, from = 0, to = 100, bw = 100 * 0.01)$y, 0)) %>%
  group_by(Species) %>%
  mutate(
    y = y * 10 / ( sum(y) * 100 / 2^10 ), # Riemann sum (range/n = x[3] - x[2])
    fill = case_when(
      Species == "Amphiroa anceps" ~ "#ec7faa",
      Species == "Jania rosea" ~ "#f5a54a",
      Species == "Metamastophora flabellata" ~ "#6b4d8d"
    )
  ) %>%
  ungroup() %T>%
  print()

# 3.8.2 Static plots ####
require(ggnewscale)
ggplot() +
  geom_textline(data = tibble(x = c(0, 0.75), y = c(0, 0.75)), aes(x, y),
                label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
  geom_polygon(data = ellipses,
               aes(x = x, y = y, fill = fill, alpha = alpha, 
                   group = group, subgroup = subgroup)) +
  geom_line(data = lines,
            aes(x = x, y = y, colour = colour)) +
  geom_polygon(data = ribbons,
               aes(x = x, y = y, fill = fill, alpha = alpha,
                   # Interaction grouping is only needed in the static version
                   group = interaction(alpha, Species))) +
  scale_alpha_identity() +
  scale_colour_identity() +
  scale_fill_identity() +
  new_scale_fill() +
  geom_polygon(data = ellipses %>% mutate(x = Inf), # Add new fill geometry for legend
               aes(x = x, y = y, fill = Species)) +
  scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d")) +
  scale_x_continuous(breaks = seq(0, 2, 0.5),
                     labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1))) +
  scale_y_continuous(breaks = seq(0, 0.75, 0.25),
                     labels = scales::label_number(accuracy = c(1, 0.01, 0.1, 0.01))) +
  labs(x = expression("CDR"["naive"]*" (mmol CO"[2]*" g"^-1*" d"^-1*")"),
       y = expression("CDR"["true"]*" (mmol CO"[2]*" g"^-1*" d"^-1*")")) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 0.75),
                  expand = FALSE, clip = "off") +
  mytheme +
  theme(legend.text = element_text(face = "italic"))

ggplot() +
  geom_polygon(data = densities,
               aes(x = x, y = y, fill = fill)) +
  scale_fill_identity() +
  xlab("CDR loss (%)") +
  coord_cartesian(expand = FALSE, clip = "off") +
  mytheme +
  theme(legend.text = element_text(face = "italic"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank())
# All looks to be in order.

# 3.8.3 Tween geometries ####
# As per usual there is an uneven sample size, so not all points can be paired.
# Define enter/exit functions
enter <- function(data){
  data$alpha <- 0
  data
}

exit <- function(data){
  data$alpha <- 0
  data
}

require(tweenr)
require(transformr)
ellipses_ani <- bind_rows(
  tween_polygon(ellipses %>% filter(Species == "Amphiroa anceps"),
                ellipses %>% filter(Species == "Jania rosea"),
                ease = "cubic-in-out", nframes = 100,
                id = id, enter = enter, exit = exit) %>%
    keep_state(nframes = 50),
  tween_polygon(ellipses %>% filter(Species == "Jania rosea"),
                ellipses %>% filter(Species == "Metamastophora flabellata"),
                ease = "cubic-in-out", nframes = 100,
                id = id, enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_polygon(ellipses %>% filter(Species == "Metamastophora flabellata"),
                ellipses %>% filter(Species == "Amphiroa anceps"),
                ease = "cubic-in-out", nframes = 100,
                id = id, enter = enter, exit = exit) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()

lines_ani <- bind_rows(
  tween_path(lines %>% filter(Species == "Amphiroa anceps"),
             lines %>% filter(Species == "Jania rosea"),
             ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_path(lines %>% filter(Species == "Jania rosea"),
             lines %>% filter(Species == "Metamastophora flabellata"),
             ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_path(lines %>% filter(Species == "Metamastophora flabellata"),
             lines %>% filter(Species == "Amphiroa anceps"),
             ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()

ribbons_ani <- bind_rows(
  tween_path(ribbons %>% filter(Species == "Amphiroa anceps"),
             ribbons %>% filter(Species == "Jania rosea"),
             ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_path(ribbons %>% filter(Species == "Jania rosea"),
             ribbons %>% filter(Species == "Metamastophora flabellata"),
             ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_path(ribbons %>% filter(Species == "Metamastophora flabellata"),
             ribbons %>% filter(Species == "Amphiroa anceps"),
             ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()

densities_ani <- bind_rows(
  tween_polygon(densities %>% filter(Species == "Amphiroa anceps"),
                densities %>% filter(Species == "Jania rosea"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50),
  tween_polygon(densities %>% filter(Species == "Jania rosea"),
                densities %>% filter(Species == "Metamastophora flabellata"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 150),
  tween_polygon(densities %>% filter(Species == "Metamastophora flabellata"),
                densities %>% filter(Species == "Amphiroa anceps"),
                ease = "cubic-in-out", nframes = 100) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 300)
) %T>%
  print()

# 6.7 Dynamic plots ####
require(gganimate)
require(ggtext)
( ggplot() +
    geom_textline(data = tibble(x = c(0, 0.75), y = c(0, 0.75)), aes(x, y),
                  label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
    geom_polygon(data = ellipses_ani,
                 aes(x = x, y = y, fill = fill, alpha = alpha, 
                     group = group, subgroup = subgroup)) +
    geom_line(data = lines_ani,
              aes(x = x, y = y, colour = colour)) +
    geom_polygon(data = ribbons_ani,
                 aes(x = x, y = y, fill = fill,
                     alpha = alpha, group = alpha)) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_fill_identity() +
    new_scale_fill() +
    geom_polygon(data = ellipses %>% mutate(x = Inf), # Add new fill geometry for legend
                 aes(x = x, y = y, fill = Species)) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d")) +
    scale_x_continuous(breaks = seq(0, 2, 0.5),
                       labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1))) +
    scale_y_continuous(breaks = seq(0, 0.75, 0.25),
                       labels = scales::label_number(accuracy = c(1, 0.01, 0.1, 0.01))) +
    labs(x = "CDR<sub><span style='font-size:8.4pt'>naive</span></sub> (mmol CO<sub><span style=
         'font-size:8.4pt'>2</span></sub> g<sup><span style='font-size:8.4pt'>−1</span></sup> 
         d<sup><span style='font-size:8.4pt'>−1</span></sup>)",
         y = "CDR<sub><span style='font-size:8.4pt'>true</span></sub> (mmol CO<sub><span style=
         'font-size:8.4pt'>2</span></sub> g<sup><span style='font-size:8.4pt'>−1</span></sup> 
         d<sup><span style='font-size:8.4pt'>−1</span></sup>)") +
    coord_cartesian(xlim = c(0, 2), ylim = c(0, 0.75),
                    expand = FALSE, clip = "off") +
    transition_manual(.frame) +
    mytheme +
    theme(axis.title = element_markdown()) ) %>%
  animate(nframes = 450, duration = 15,
          width = 21 * 1/2, height = 10,
          units = "cm", res = 300, renderer = gifski_renderer()) %>%
  anim_save(filename = "Fig_5a.gif", path = "Figures")
# The legend doesn't fit in the frame. Plot legend separately.
# Also make plot a bit shorter.

( ggplot() +
    geom_polygon(data = ellipses %>% mutate(x = Inf, y = Inf),
                 aes(x = x, y = y, fill = Species)) +
    scale_fill_manual(values = c("#ec7faa", "#f5a54a", "#6b4d8d")) +
    mytheme +
    theme(legend.text = element_text(face = "italic"),
          axis.title = element_blank()) ) %>%
  ggsave(filename = "Fig_5_legend.png", path = "Figures",
         units = "cm", width = 21, height = 1, dpi = 300)

( ggplot() +
    geom_textline(data = tibble(x = c(0, 0.75), y = c(0, 0.75)), aes(x, y),
                  label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
    geom_polygon(data = ellipses_ani,
                 aes(x = x, y = y, fill = fill, alpha = alpha, 
                     group = group, subgroup = subgroup)) +
    geom_line(data = lines_ani,
              aes(x = x, y = y, colour = colour)) +
    geom_polygon(data = ribbons_ani,
                 aes(x = x, y = y, fill = fill,
                     alpha = alpha, group = alpha)) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_x_continuous(breaks = seq(0, 2, 0.5),
                       labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1))) +
    scale_y_continuous(breaks = seq(0, 0.75, 0.25),
                       labels = scales::label_number(accuracy = c(1, 0.01, 0.1, 0.01))) +
    labs(x = "CDR<sub><span style='font-size:8.4pt'>naive</span></sub> (mmol CO<sub><span style=
         'font-size:8.4pt'>2</span></sub> g<sup><span style='font-size:8.4pt'>−1</span></sup> 
         d<sup><span style='font-size:8.4pt'>−1</span></sup>)",
         y = "CDR<sub><span style='font-size:8.4pt'>true</span></sub> (mmol CO<sub><span style=
         'font-size:8.4pt'>2</span></sub> g<sup><span style='font-size:8.4pt'>−1</span></sup> 
         d<sup><span style='font-size:8.4pt'>−1</span></sup>)") +
    coord_cartesian(xlim = c(0, 2), ylim = c(0, 0.75),
                    expand = FALSE, clip = "off") +
    transition_manual(.frame) +
    mytheme +
    theme(axis.title = element_markdown()) ) %>%
  animate(nframes = 450, duration = 15,
          width = 21 * 1/2, height = 9,
          units = "cm", res = 300, renderer = gifski_renderer()) %>%
  anim_save(filename = "Fig_5a.gif", path = "Figures")

( ggplot() +
    geom_polygon(data = densities_ani,
                 aes(x = x, y = y, fill = fill)) +
    scale_fill_identity() +
    coord_cartesian(expand = FALSE, clip = "off") +
    xlab("CDR loss (%)<sub><span style='color:white;font-size:8.4pt'>2</span></sub>
         <sup><span style='color:white;font-size:8.4pt'>−1</span></sup>") +
    transition_manual(.frame) +
    mytheme +
    theme(axis.title = element_markdown(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) ) %>%
  animate(nframes = 450, duration = 15,
          width = 21 * 1/2, height = 9,
          units = "cm", res = 300, renderer = gifski_renderer()) %>%
  anim_save(filename = "Fig_5b.gif", path = "Figures")
