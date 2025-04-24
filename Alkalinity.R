# 1. Load data ####
# Load output from Metrohm Titrando PC control V1
require(tidyverse)
require(here)
files <- here("Alkalinity") %>% list.files(pattern = "\\.txt$", full.names = TRUE)

TA <- files %>%
  map(~ read.delim(.x, skip = 27, header = FALSE) %>% # read text file, skipping some rows
        mutate(pH_initial = V2[V1 == "Initial pH"] %>% as.numeric(), # extract relevant numbers
               HCl_M = V2[V1 == "CONC"] %>% as.numeric(), # in the lower section of the text file
               Name = V3[V1 == "Total_Alkalinity"][1] %>%
                 str_split_i(pattern = "-", i = 1),
               Date_time = V4[V1 == "Total_Alkalinity"][1] %>%
                 ymd_hms(tz = "Asia/Singapore"),
               pH_slope = V2[V1 == "pH electrode"] %>% as.numeric(),
               pH_cal = V6[V1 == "pH electrode"] %>% 
                 ymd_hms(tz = "Asia/Singapore")) %>%
        separate_wider_delim(Name, delim = "_", names = c("Date", "ID"), 
                             too_many = "merge", too_few = "align_end") %>%
        slice(1:(which(V2 == "")[1] - 1)) %>% # remove the lower section of the text file
        mutate(across(c(V1, V2, V3, V4, V5, V6), as.numeric),
               Date = Date %>% ymd(), ID = ID %>% fct()) %>%
        rename(n = V1, t = V2, pH = V3, mL = V4, µL_min = V5, Temp = V6)) %>%
  set_names(files %>% basename() %>% # extract list names from file names
              str_remove("PC_LIMS_Report-") %>%
              str_remove("\\.txt$") %>%
              str_split_i(pattern = "-", i = 1) %>%
              make.names) %>%  
  imap(~ .x %>% mutate(Date_file = if_else( # extract Date and ID again from list names
                         .y %>% str_detect("_"),
                         .y %>%
                           str_split_i(pattern = "_", i = 1) %>% 
                           str_sub(start = 2) %>%
                           ymd(),
                         NA_Date_),
                        ID_file = if_else(
                          .y %>% str_detect("_"),
                          str_split(.y, pattern = "_", n = 2)[[1]][2],
                          .y) %>% fct()
                       )
       )

str(TA$X240524_0_1_D) # examples
str(TA$X241204_A_13_D)

# Check redundancy of Date, Date_file and ID_file
TA %>%
  map(~ .x %>% 
        mutate(Date_check_1 = identical(Date, Date_file),
               Date_check_2 = identical(Date_file, Date_time %>% as_date()),
               ID_check = identical(ID, ID_file))
      ) %>%
  bind_rows() %>%
  group_by(ID, ID_file, Date, Date_file, Date_time) %>%
  summarise(Date_check_1 = unique(Date_check_1),
            Date_check_2 = unique(Date_check_2),
            ID_check = unique(ID_check)) %>%
  print(n = 182)
# There are some cases that don't check out. These cases are either typos, 
# e.g. O for 0, wrong species identifier, date written dmy instead of ymd,
# and round and individual numbers are swapped, or cases where date was not
# recorded in the file name. ID_file is correct, so ID can be removed. 
# Date_file is correct, so Date can be removed. Since Date_time as date is 
# identical to Date_file, Date_file can also be removed, leaving ID_file
# and Date_time.

# Remove Date, Date_file and ID
require(magrittr)
TA %<>%
  map(~ .x %>% 
        select(-c(Date, Date_file, ID)) %>%
        rename(ID = ID_file))

# Compare time to last calibration
TA %>%
  map(~ .x %>% 
        mutate(Time = pH_cal %--% Date_time %>% 
                 as.duration() / ddays())
      ) %>%
  bind_rows() %>%
  group_by(Date_time, ID) %>%
  summarise(Time = mean(Time)) %>%
  arrange(Time) %>%
  print(n = 182)
# Never more than 2 days between calibration and measurement; 
# mostly done on the same day.

# Calculate pH slope summary statistics
TA %>%
  bind_rows() %>%
  summarise(pH_slope_min = min(pH_slope),
            pH_slope_max = max(pH_slope),
            pH_slope_mean = mean(pH_slope),
            pH_slope_sd = sd(pH_slope))


# 2. Raw models ####
# 2.1 Visualise data ####
require(patchwork)
TA %>%
  imap(~ .x %>%
        ggplot() +
          geom_hline(yintercept = 3) + # endpoint pH
          geom_point(aes(t, pH), shape = 16, alpha = 0.2) +
          theme_minimal() +
          theme(panel.grid = element_blank()) +
          ggtitle(.y)
       ) %>%
  wrap_plots() %>%
  ggsave(filename = "TA_data.pdf", path = "Plots", 
         width = 100, height = 50, unit = "cm", device = cairo_pdf)
# Most curves look fine but A_6_L, J_6_L and M_6_L were overacidified, likely because
# their TAs were very low. A_9_L never had any acid added so its pH remains unchanged.
# TA cannot directly be estimated for any of these samples. However, given an average
# pH vs. acid volume slope, one could get estimates for the overacidified samples.
# The initial pH is also still valid for the sample that had no acid added to it.
# Save these data to a separate list:
TA_fails <- TA %>%
  map(~ .x %>% 
        filter(ID %in% c("A_6_L", "J_6_L", "M_6_L", "A_9_L"))
      ) %>%
  keep(~ nrow(.x) > 0)

# Clean data
TA %<>%
  map(~ .x %>% 
        filter(!ID %in% c("A_6_L", "J_6_L", "M_6_L", "A_9_L"))
  ) %>%
  keep(~ nrow(.x) > 0)
  
# Visualise clean data
TA %>%
  imap(~ .x %>%
        ggplot() +
          geom_hline(yintercept = 3) + # endpoint pH
          geom_point(aes(t, pH), shape = 16, alpha = 0.2) +
          theme_minimal() +
          theme(panel.grid = element_blank()) +
          ggtitle(.y)
       ) %>%
  wrap_plots() %>%
  ggsave(filename = "TA_data_cleaned.pdf", path = "Plots", 
         width = 100, height = 50, unit = "cm", device = cairo_pdf)
# Looks fine

# 2.2 Transform data ####
# pH responds non-linearly with increasing acid volume:
TA %>%
  imap(~ .x %>%
        ggplot() +
          geom_hline(yintercept = 3) + # endpoint pH
          geom_point(aes(mL, pH), shape = 16, alpha = 0.2) +
          theme_minimal() +
          theme(panel.grid = element_blank()) +
          ggtitle(.y)
       ) %>%
  wrap_plots() %>%
  ggsave(filename = "TA_data_mL.pdf", path = "Plots", 
         width = 100, height = 50, unit = "cm", device = cairo_pdf)
# pH therefore needs to be transformed to [H+] using the Gran function. 
# Added acid volume correspondingly needs to be converted to µM or µeq
# (HCl is monoprotic so N and M are equivalent). Expressing both variables
# in units of µM in terms of sample volume makes their relationship  
# interpretable.

TA %<>%
  map(~ .x %>%
        # free [H+] in µM = total volume in µL * [H+] in M / sample volume in L
        mutate(H_free = ( 10 + mL ) * 1e3 * 10^-pH / 0.01, # all samples had a volume of 10 mL
               # added [H+] in µM = added acid in µmol / sample volume in L
               H_added = mL * HCl_M * 1e3 / 0.01)
  )

TA %>%
  imap(~ .x %>%
        ggplot() +
          geom_hline(yintercept = 0) + # zero free protons
          geom_point(aes(H_added, H_free), shape = 16, alpha = 0.2) +
          theme_minimal() +
          theme(panel.grid = element_blank()) +
          ggtitle(.y)
       ) %>%
  wrap_plots() %>%
  ggsave(filename = "TA_data_transformed.pdf", path = "Plots", 
         width = 100, height = 50, unit = "cm", device = cairo_pdf)

# As is typical for the Gran plot, due to buffering effect of seawater, 
# only free protons increase linearly with total protons added. To fit
# a linear model, only the portion of data that are still nonlinear
# need to be removed. Rather than doing this arbitrarily on a case-by-case
# basis, we can set a cutoff pH. We know the equivalence point to lie
# around pH 4 since all CO3 2- and HCO3- are fully converted to CO2.
# But since there is often not a smooth transition from buffered to
# free states around pH 4, which is precisely why we are using the Gran
# function to estimate the equivalence point, I will set a lower cutoff 
# of pH 3.5.

TA %>%
  imap(~ .x %>%
        filter(pH <= 3.5) %>% # filter pH to be equal to or less than 3.5
        ggplot() +
          geom_hline(yintercept = 0) + # zero free protons
          geom_point(aes(H_added, H_free), shape = 16, alpha = 0.2) +
          theme_minimal() +
          theme(panel.grid = element_blank()) +
          ggtitle(.y)
       ) %>%
  wrap_plots() %>%
  ggsave(filename = "TA_data_linear.pdf", path = "Plots", 
         width = 100, height = 50, unit = "cm", device = cairo_pdf)
# Looks linear in all cases

# 2.3 Prepare data ####
# Since I want to estimate an average slope of free protons over added protons 
# in solution with which to salvage three of the four TA_fails, I need a multi-
# level model. All data need to be combined into one tibble to achieve this. 
# My key assumption is that the increase in free protons with protons added
# beyond the buffering capacity is comparable across samples with different
# buffering capacities.

TA_data <- TA %>% 
  map(~ .x %>% filter(pH <= 3.5) ) %>%
  bind_rows() %>%
  filter(ID != "standard") %>%
  mutate(ID = ID %>% fct_drop())

# Standards need to be estimated in a separate multilevel model since I am
# interested in deriving an average standard TA. Standards need unique IDs 
# for this step.

TA_standard_data <- TA %>% 
  map(~ .x %>% filter(pH <= 3.5) ) %>%
  bind_rows() %>%
  filter(ID == "standard") %>%
  mutate(ID = ID %>% 
           str_c(Date_time, sep = "_") %>% 
           fct())

# 2.4 Sample model ####
# 2.4.1 Prior simulation ####
# The typical linear regression is parameterised as alpha + beta * x, but in
# this case alpha is not meaningful. It is bound to be negative because the
# x-intercept is always positive and negative free protons are not useful.
# The x-intercept is our estimate of total alkalinity (TA) in µM because
# it gives [H+] in µM at the equivalence point (transition from buffered
# to non-buffered state) and is therefore best estimated directly. f(0), i.e.
# all protons are buffered, thus yields 0 = alpha + beta * TA which allows
# us to re-express alpha as -beta * TA. The reparameterised function therefore
# is -beta * TA + beta * x = beta * ( x - TA ). Now the intercept term TA is
# directly interpretable: for average seawater we expect 2300 µM. beta is
# theoretically expected to lie around 1 since in the non-buffered state
# all added protons are expected to be free, i.e. 1:1 relationship between
# added and free protons.

tibble(n = 1:1e3,
       TA = rgamma(n = 1e3, shape = 2300^2 / 500^2, rate = 2300 / 500^2),
       beta = rgamma(n = 1e3, shape = 1^2 / 0.5^2, rate = 1 / 0.5^2)) %>%
  expand_grid(H_added = seq(0, 6000)) %>%
  mutate(H_free = beta * ( H_added - TA )) %>%
  ggplot(aes(H_added, H_free, group = n)) +
    geom_hline(yintercept = 0) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable.

# 2.4.2 Run model ####
TA_stan <- "
data{
  int n;
  vector<lower=0>[n] H_free;
  vector<lower=0>[n] H_added;
  array[n] int ID;
  int n_ID;
}

parameters{
  // Hyperparameters
  real<lower=0> beta_mu;
  real<lower=0> beta_sigma;
  
  // Titration-specific parameters
  vector<lower=0>[n_ID] TA;
  vector<lower=0>[n_ID] beta;
  
  // Likelihood uncertainty
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  beta_mu ~ gamma( 1^2 / 0.5^2 , 1 / 0.5^2 ); // reparameterised with mean and sd
  beta_sigma ~ exponential( 1 );
  
  // Titration-specific priors
  TA ~ gamma( 2300^2 / 500^2 , 2300 / 500^2 ); 
  beta ~ gamma( beta_mu^2 / beta_sigma^2 , beta_mu / beta_sigma^2 );
  
  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
      mu[i] = beta[ID[i]] * ( H_added[i] - TA[ID[i]] );
  }

  // Likelihood
  H_free ~ normal( mu , sigma );
}
"

require(cmdstanr)
TA_mod <- TA_stan %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
TA_samples <- TA_mod$sample(
  data = TA_data %>%
    select(H_free, H_added, ID) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4)

# 2.4.3 Model checks ####
TA_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

require(bayesplot)
TA_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "TA_rank.pdf", path = "Plots",
         width = 100, height = 50, unit = "cm", device = cairo_pdf)
# chains look good

TA_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("beta[1]", "TA[1]"))
TA_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("beta[100]", "TA[100]"))
# parameters are somewhat correlated but this is not an issue for the sampler

# 2.4.4 Prior-posterior comparison ####
source("functions.R")
# sample priors
TA_prior <- prior_samples(
  model = TA_mod,
  data = TA_data %>%
    select(H_free, H_added, ID) %>%
    compose_data(),
  chains = 8, samples = 1e4)

# plot prior-posterior comparison
TA_prior %>%
  prior_posterior_draws(posterior_samples = TA_samples,
                        group = TA_data %>% select(ID),
                        parameters = c("beta[ID]", "TA[ID]",
                                       "beta_mu", "beta_sigma", 
                                       "sigma"),
                        format = "long") %>%
  prior_posterior_plot(group_name = "ID", ridges = FALSE) %>%
  ggsave(filename = "TA_prior_posterior.pdf", path = "Plots",
         width = 80, height = 80, unit = "cm", device = cairo_pdf)
# posteriors are very well constrained

# 2.4.5 Predictions ####
# Split grouped tibble into list of tibbles for smoother computation
TA_data_list <- TA_data %>%
  group_by(ID) %>%
  group_split() %>%
  set_names(
    map(., ~ .x %$% 
          as.character(ID) %>% 
          unique()
    )
  )

TA_prior_posterior_list <- TA_prior %>%
  prior_posterior_draws(posterior_samples = TA_samples,
                        group = TA_data %>% select(ID),
                        parameters = c("beta[ID]", "TA[ID]", "sigma"),
                        format = "short") %>%
  group_by(ID) %>%
  group_split() %>%
  set_names(
    map(., ~ .x %$% 
          as.character(ID) %>% 
          unique()
        )
  )

TA_predictions_list <- TA_prior_posterior_list %>%
  map2(TA_data_list,
       ~ .x %>% expand_grid(H_added = .y %$%
                              seq(.x %>% # lower bound is minimum x-intercept
                                    filter(distribution == "posterior") %$%
                                    min(TA), 
                                  max(H_added), # upper bound is maximum [H+]
                                  length.out = 15))
       ) %>%
  map(~ .x %>%
        mutate(mu = beta * ( H_added - TA ),
               obs = rnorm( n(), mu, sigma ))
      )

TA_predictions_summary <- TA_predictions_list %>%
  map(~ .x %>%
        group_by(ID, distribution, H_added) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_")
      ) %>%
  bind_rows()

# Remove raw predictions because they are too big
rm(TA_data_list, TA_prior_posterior_list, TA_predictions_list) 

ggsave(
  TA_predictions_summary %>%
    ggplot() +
      geom_hline(yintercept = 0) +
      geom_point(data = TA_data, 
                 aes(H_added, H_free),
                 shape = 16, alpha = 0.1) +
      geom_line(data = . %>% filter(distribution == "posterior"),
                aes(H_added, mu_y)) +
      geom_ribbon(data = . %>% filter(distribution == "posterior"),
                  aes(H_added, ymin = mu_ymin, ymax = mu_ymax,
                      alpha = factor(mu_.width))) +
      # geom_ribbon(data = . %>% filter(distribution == "posterior"), # unhash to check
      #             aes(H_added, ymin = obs_ymin, ymax = obs_ymax, # predicted observations
      #                 alpha = factor(obs_.width))) +
      # geom_ribbon(data = . %>% filter(distribution == "prior", mu_.width == 0.9),
      #             aes(H_added, ymin = mu_ymin, ymax = mu_ymax), # unhash to check prior
      #                 colour = alpha("black", 0.3), fill = NA) +
      scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
      facet_wrap(~ ID, scales = "free") +
      theme_minimal() +
      theme(panel.grid = element_blank()),
  filename = "TA_prediction.pdf", path = "Plots",
  width = 100, height = 50, unit = "cm", device = cairo_pdf)
# Linear model fits well in all cases.

# 2.4.6 Visualisations for supplement ####
# Pull out a few examples for the supplement (A_16_D, J_4_D, M_12_D)
TA_predictions_summary_selection <- TA_predictions_summary %>%
  filter(ID %in% c("A_16_D", "J_4_D", "M_12_D")) %>%
  mutate(Species = case_when(
                      ID %>% str_detect("A") ~ "Amphiroa anceps",
                      ID %>% str_detect("J") ~ "Jania rosea",
                      ID %>% str_detect("M") ~ "Metamastophora flabellata"
                      ) %>% fct_relevel("Amphiroa anceps")
         )

TA_data_selection <- TA_data %>% # these are the data the model was conditioned on
  filter(ID %in% c("A_16_D", "J_4_D", "M_12_D")) %>%
  mutate(Species = case_when(
                      ID %>% str_detect("A") ~ "Amphiroa anceps",
                      ID %>% str_detect("J") ~ "Jania rosea",
                      ID %>% str_detect("M") ~ "Metamastophora flabellata"
                      ) %>% fct_relevel("Amphiroa anceps")
         )

TA_data_selection_raw <- TA %$% # these are the unfiltered data (full pH range)
  bind_rows(X250227_A_16_D, X240621_J_4_D, X241122_M_12_D) %>%
  mutate(Species = case_when(
                      ID %>% str_detect("A") ~ "Amphiroa anceps",
                      ID %>% str_detect("J") ~ "Jania rosea",
                      ID %>% str_detect("M") ~ "Metamastophora flabellata"
                      ) %>% fct()
         )

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

Fig_S2a <- 
  TA_data_selection_raw %>%
    ggplot() +
      geom_point(aes(mL / (mL + 10) * 100, pH), # re-express added acid as % v/v 
                 shape = 16, alpha = 0.1) +
      facet_grid(~ Species) +
      labs(x = "0.012 M HCl (% v/v)", 
           y = expression("pH"["F"])) +
      scale_x_continuous(breaks = seq(16, 28, 3)) +
      scale_y_continuous(breaks = seq(3, 5, 0.5),
                         labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1))) +
      coord_cartesian(xlim = c(16, 28), ylim = c(3, 5),
                      clip = "off", expand = FALSE) +
      mytheme +
      theme(strip.text = element_text(face = "italic"))

Fig_S2a %>%
  ggsave(filename = "Fig_S2a.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

Fig_S2b <- 
  TA_data_selection_raw %>%
    ggplot() +
      geom_point(aes(H_added, H_free),
                 shape = 16, alpha = 0.1) +
      facet_grid(~ Species) +
      labs(x = expression("[H"^"+"*"]"["A"]*" (µM)"), 
           y = expression("[H"^"+"*"]"["F"]*" (µM)")) +
      coord_cartesian(xlim = c(2e3, 4.5e3), ylim = c(0, 1.5e3),
                      clip = "off", expand = FALSE) +
      mytheme +
      theme(strip.text = element_text(face = "italic"),
            panel.spacing = unit(1, "cm"))

Fig_S2b %>%
  ggsave(filename = "Fig_S2b.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

Fig_S2c <- 
  TA_data_selection %>%
    ggplot() +
      geom_point(aes(H_added, H_free),
                 shape = 16, alpha = 0.1) +
      facet_grid(~ Species) +
      labs(x = expression("[H"^"+"*"]"["A"]*" (µM)"), 
           y = expression("[H"^"+"*"]"["F"]*" (µM)")) +
      coord_cartesian(xlim = c(2e3, 4.5e3), ylim = c(0, 1.5e3),
                      clip = "off", expand = FALSE) +
      mytheme +
      theme(strip.text = element_text(face = "italic"),
            panel.spacing = unit(1, "cm"))

Fig_S2c %>%
  ggsave(filename = "Fig_S2c.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

Fig_S2d <- 
  TA_predictions_summary_selection %>%
    ggplot() +
      geom_point(data = TA_data_selection, 
                 aes(H_added, H_free),
                 shape = 16, alpha = 0.1) +
      geom_line(data = . %>% filter(distribution == "posterior"),
                aes(H_added, mu_y), linewidth = 0.2) +
      geom_ribbon(data = . %>% filter(distribution == "posterior"), 
                  aes(H_added, ymin = mu_ymin, ymax = mu_ymax,
                      alpha = factor(mu_.width))) +
      scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
      facet_grid(~ Species) +
      labs(x = expression("[H"^"+"*"]"["A"]*" (µM)"), 
           y = expression("[H"^"+"*"]"["F"]*" (µM)")) +
      coord_cartesian(xlim = c(2e3, 4.5e3), ylim = c(0, 1.5e3),
                      clip = "off", expand = FALSE) +
      mytheme +
      theme(strip.text = element_text(face = "italic"),
            panel.spacing = unit(1, "cm"))

Fig_S2d %>%
  ggsave(filename = "Fig_S2d.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# 2.5 Standard model ####
# Priors as above
# 2.5.1 Run model ####
TA_standard_stan <- "
data{
  int n;
  vector<lower=0>[n] H_free;
  vector<lower=0>[n] H_added;
  array[n] int ID;
  int n_ID;
}

parameters{
  // Hyperparameters
  real<lower=0> TA_mu;
  real<lower=0> TA_sigma;
  
  // Titration-specific parameters
  vector<lower=0>[n_ID] TA;
  vector<lower=0>[n_ID] beta;
  
  // Likelihood uncertainty
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  TA_mu ~ gamma( 2300^2 / 500^2 , 2300 / 500^2 ); // reparameterised with mean and sd
  TA_sigma ~ exponential( 1 );
  
  // Titration-specific priors
  TA ~ gamma( TA_mu^2 / TA_sigma^2 , TA_mu / TA_sigma^2 ); 
  beta ~ gamma( 1^2 / 0.5^2 , 1 / 0.5^2 );
  
  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
      mu[i] = beta[ID[i]] * ( H_added[i] - TA[ID[i]] );
  }

  // Likelihood
  H_free ~ normal( mu , sigma );
}
"

TA_standard_mod <- TA_standard_stan %>%
  write_stan_file() %>%
  cmdstan_model()

TA_standard_samples <- TA_standard_mod$sample(
  data = TA_standard_data %>%
    select(H_free, H_added, ID) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4)

# 2.5.2 Model checks ####
TA_standard_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat), # proportion > 1.001
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# no rhat above 1.001
# good effective sample size

TA_standard_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(filename = "TA_standard_rank.pdf", path = "Plots",
         width = 40, height = 20, unit = "cm", device = cairo_pdf)
# chains look good

TA_standard_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("beta[1]", "TA[1]"))
TA_standard_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("beta[10]", "TA[10]"))
# parameters are somewhat correlated but this is not an issue for the sampler

# 2.5.3 Prior-posterior comparison ####
# sample priors
TA_standard_prior <- prior_samples(
  model = TA_standard_mod,
  data = TA_standard_data %>%
    select(H_free, H_added, ID) %>%
    compose_data(),
  chains = 8, samples = 1e4)

# plot prior-posterior comparison
TA_standard_prior %>%
  prior_posterior_draws(posterior_samples = TA_standard_samples,
                        group = TA_standard_data %>% select(ID),
                        parameters = c("beta[ID]", "TA[ID]",
                                       "TA_mu", "TA_sigma", 
                                       "sigma"),
                        format = "long") %>%
  prior_posterior_plot(group_name = "ID", ridges = FALSE) %>%
  ggsave(filename = "TA_standard_prior_posterior.pdf", path = "Plots",
         width = 80, height = 40, unit = "cm", device = cairo_pdf)
# posteriors are very well constrained

# 2.5.4 Predictions ####
# Split grouped tibble into list of tibbles for smoother computation
TA_standard_data_list <- TA_standard_data %>%
  group_by(ID) %>%
  group_split() %>%
  set_names(
    map(., ~ .x %$% 
          as.character(ID) %>% 
          unique()
    )
  )

TA_standard_prior_posterior_list <- TA_standard_prior %>%
  prior_posterior_draws(posterior_samples = TA_standard_samples,
                        group = TA_standard_data %>% select(ID),
                        parameters = c("beta[ID]", "TA[ID]", "sigma"),
                        format = "short") %>%
  group_by(ID) %>%
  group_split() %>%
  set_names(
    map(., ~ .x %$% 
          as.character(ID) %>% 
          unique()
        )
  )

TA_standard_predictions_list <- TA_standard_prior_posterior_list %>%
  map2(TA_standard_data_list,
       ~ .x %>% expand_grid(H_added = .y %$%
                              seq(.x %>% # lower bound is minimum x-intercept
                                    filter(distribution == "posterior") %$%
                                    min(TA), 
                                  max(H_added), # upper bound is maximum [H+]
                                  length.out = 15))
       ) %>%
  map(~ .x %>%
        mutate(mu = beta * ( H_added - TA ),
               obs = rnorm( n(), mu, sigma ))
      )

TA_standard_predictions_summary <- TA_standard_predictions_list %>%
  map(~ .x %>%
        group_by(ID, distribution, H_added) %>%
        reframe(mu = mu %>% mean_qi(.width = c(.5, .8, .9)),
                obs = obs %>% mean_qi(.width = c(.5, .8, .9))) %>%
        unnest(c(mu, obs), names_sep = "_")
      ) %>%
  bind_rows()

# Remove raw predictions because they are too big
rm(TA_standard_data_list, TA_standard_prior_posterior_list, TA_standard_predictions_list) 

ggsave(
  TA_standard_predictions_summary %>%
    ggplot() +
      geom_hline(yintercept = 0) +
      geom_point(data = TA_standard_data, 
                 aes(H_added, H_free),
                 shape = 16, alpha = 0.1) +
      geom_line(data = . %>% filter(distribution == "posterior"),
                aes(H_added, mu_y)) +
      geom_ribbon(data = . %>% filter(distribution == "posterior"),
                  aes(H_added, ymin = mu_ymin, ymax = mu_ymax,
                      alpha = factor(mu_.width))) +
      # geom_ribbon(data = . %>% filter(distribution == "posterior"), # unhash to check
      #             aes(H_added, ymin = obs_ymin, ymax = obs_ymax, # predicted observations
      #                 alpha = factor(obs_.width))) +
      # geom_ribbon(data = . %>% filter(distribution == "prior", mu_.width == 0.9),
      #             aes(H_added, ymin = mu_ymin, ymax = mu_ymax), # unhash to check prior
      #                 colour = alpha("black", 0.3), fill = NA) +
      scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
      facet_wrap(~ ID, scales = "free") +
      theme_minimal() +
      theme(panel.grid = element_blank()),
  filename = "TA_standard_prediction.pdf", path = "Plots",
  width = 40, height = 20, unit = "cm", device = cairo_pdf)
# Linear model fits well in all cases.

# 3. Salvage failed TA ####
# Summarise salvageable titration data
TA_salvaged <- TA_fails %>%
  bind_rows() %>%
  filter(ID != "A_9_L") %>% # filter out measurement that had no acid added
  group_by(ID) %>%
  summarise(pH = mean(pH), # summarise overacidified pH measurements
            mL = mean(mL),
            Temp = unique(Temp),
            pH_initial = unique(pH_initial),
            HCl_M = unique(HCl_M)) %>%
  mutate(H_free = ( 10 + mL ) * 1e3 * 10^-pH / 0.01, # apply Gran function as before
         H_added = mL * HCl_M * 1e3 / 0.01)

TA_salvaged

# Estimate TA using hyperparameters
# The linear equation used before is H_free = beta * ( H_added - TA ), which solved 
# for TA becomes TA = H_added - H_free / beta. In overacidified samples H_added and 
# H_free are measured for a point beyond the endpoint of pH 3, and a global beta can
# be estimated using the hyperparameters beta_mu and beta_sigma, so we have all the
# ingredients to estimate TA with the adequate uncertainty.

TA_salvaged %<>%
  cross_join(
    TA_samples %>%
      spread_draws(beta_mu, beta_sigma) %>% # estimate grand mean for beta
      mutate(beta = rgamma( n() , beta_mu^2 / beta_sigma^2 , beta_mu / beta_sigma^2 ))
    ) %>%
  mutate(TA = H_added - H_free / beta) %>% # calculate TA as describe above
  select(ID, Temp, pH_initial, starts_with("."), TA)

TA_salvaged

# Combine all pH and TA data
pH_TA <- TA_data %>%
  group_by(ID) %>%
  summarise(Temp = unique(Temp),
            pH_initial = unique(pH_initial)) %>%
  left_join(
    TA_samples %>%
      recover_types(TA_data %>% select(ID)) %>%
      spread_draws(TA[ID]),
    by = "ID"
  ) %>%
  bind_rows(TA_salvaged)

# Check that merges worked correctly
pH_TA %>%
  group_by(ID) %>%
  summarise(Temp = unique(Temp),
            pH_initial = unique(pH_initial),
            TA = mean(TA)) %>%
  print(n = 169)

# 4. Conversion ####
# 4.1 Correct TA with standard ####
# The standard used is Scripps Institution of Oceanography CRM #205 - 0344
# with the following properties:
# Salinity = 33.443‰
# TA = 2202.05 ± 0.98 μmol kg^–1 (mean ± s.d.)
# DIC = 2011.85 ± 0.99 μmol kg^–1 (mean ± s.d.)
# Known TA is given in µmol kg^-1 rather than µM, so needs to be converted
# to µM by multiplying by seawater density (kg L^-1), 
# i.e. µmol kg^-1 * kg L^-1 = µmol L^-1 = µM.

require(seacarb)
TA_standard <- TA_standard_samples %>%
  spread_draws(TA_mu, TA_sigma) %>%
  mutate(TA_measured = rgamma( n() , TA_mu^2 / TA_sigma^2 , TA_mu / TA_sigma^2 ),
         TA_true = rnorm( n() , 2202.05 , 0.98 ) * rho(S = 33.443, T = 25)[1] * 1e-3,
         TA_ratio = TA_true / TA_measured,
         TA_diff = TA_measured - TA_true)

TA_standard %>%
  summarise(across(
              .cols = c(TA_measured, TA_true, TA_ratio, TA_diff),
              .fns = list(mean = mean, sd = sd)
              )) %>%
  pivot_longer(cols = everything(), names_to = c("variable", "statistic"), 
               values_to = "value", names_pattern = "^(.*)_(mean|sd)$") %>%
  pivot_wider(values_from = value, names_from = statistic)

# Subtract difference between measured and expected from sample TA estimates
pH_TA %<>%
  left_join(
    TA_standard %>%
      select(starts_with("."), TA_diff),
    by = c(".chain", ".iteration", ".draw")
  ) %>%
  mutate(TA_corrected = TA - TA_diff) %>%
  select(-c(TA, TA_diff))

# 4.2 DIC from pH and TA ####
# Get salinity data
incu_meta <- read.csv("Incubation.csv") %>%
  mutate(ID = ID %>% fct(),
         Round = if_else(
                    str_length(ID) > 6,
                    ID %>% 
                      str_split(pattern = "_", n = 4) %>% 
                      map_chr(2),
                    ID %>% 
                      str_split(pattern = "_", n = 3) %>% 
                      map_chr(2)
                    ) %>% fct(),
         Species = case_when(
                      is.na(Species) & ID %>% str_detect("^0") ~ "Initial",
                      is.na(Species) & ID %>% str_detect("B") ~ "Blank",
                      TRUE ~ Species
                      ) %>% fct_relevel("Initial", "Blank", "Amphiroa anceps"),
         Treatment = Treatment %>% fct())
  
pH_TA %<>%
  left_join(
    incu_meta %>% 
      select(ID, Round, Species, Treatment, Salinity),
    by = "ID"
    )

# Calculate DIC
DIC_TA <- pH_TA %>% # seacarb::carb() takes TA in mol kg^-1,
  mutate(DIC_calculated = carb(flag = 8, 
                               var1 = pH_initial, # so the input needs to be converted from µM to mol kg^-1
                               var2 = TA_corrected * 1e-6 / ( rho(S = Salinity, T = Temp)[1] * 1e-3 ),
                               S = Salinity, 
                               T = Temp, # and the output needs to be back-converted from mol kg^-1 to µM
                               pHscale = "F")$DIC * 1e6 * rho(S = Salinity, T = Temp)[1] * 1e-3)

rm(pH_TA)

# Viusalise TA and DIC using Deffeyes diagram
# Generate values for pCO2 contours with µM TA and DIC input
pCO2_contour <- DIC_TA %$% # simulate TA and DIC based on empirical TA and DIC ranges in µM
  expand_grid(TA_sim = seq(floor( min(TA_corrected) / 100 ) * 100, # round down to nearest 100
                           ceiling( max(TA_corrected) / 100 ) * 100, # round up to nearest 100
                           length.out = 1e3),
              DIC_sim = seq(floor( min(DIC_calculated) / 100 ) * 100,
                            ceiling( max(DIC_calculated) / 100 ) * 100,
                            length.out = 1e3)) %>%
  mutate(pCO2_sim = DIC_TA %$% 
                      carb(flag = 15, # convert µM * 1e-6 = M and M / kg L^-1 = mol L^-1 * L kg^-1 = mol kg^-1
                           var1 = TA_sim * 1e-6 / ( rho(S = mean(Salinity), T = unique(Temp))[1] * 1e-3 ), 
                           var2 = DIC_sim * 1e-6 / ( rho(S = mean(Salinity), T = unique(Temp))[1] * 1e-3 ), 
                           S = mean(Salinity),
                           T = unique(Temp))$pCO2) # pH does not factor in the calculation, so no pH scale defined
require(geomtextpath)
require(ggdensity)
Fig_S3 <- DIC_TA %>%
  ggplot() +
    geom_textcontour(data = pCO2_contour, aes(x = DIC_sim, y = TA_sim, z = log10(pCO2_sim)),
                     breaks = c(seq(-2.5, 4, 0.5), 4.2, seq(4.4, 5, 0.1)), colour = "#c9d2d7",
                     size = 3.5, family = "Futura") + # 3.5 is equivalent to 10 pt
    geom_hdr(data = . %>%
               group_by(ID) %>% # randomly reshuffle within ID
               mutate(DIC_calculated = DIC_calculated %>% sample(),
                      Modified = if_else(
                                  Species %in% c("Blank", "Initial"),
                                  "Control", "Coralline"
                                  )
                      ),
             aes(DIC_calculated, TA_corrected, fill = Modified, group = ID),
             alpha = 0.5, n = 500, method = "mvnorm", probs = 0.999) +
    geom_hline(yintercept = Inf, lineend = "square") +
    geom_vline(xintercept = Inf) +
    annotate("segment", x = 2700, y = 1700, xend = 2700 + 500, yend = 1700 + 500 * 2,
             linewidth = 0.5, lineend = "square", linejoin = "mitre",
             arrow = arrow(angle = 20, length = unit(0.3, "cm"), ends = "both", type = "closed")) +
    annotate("segment", x = 2700, y = 2200, xend = 2700 + 500, yend = 2200,
             linewidth = 0.5, lineend = "square", linejoin = "mitre",
             arrow = arrow(angle = 20, length = unit(0.3, "cm"), ends = "both", type = "closed")) +
    annotate("label", x = 330, y = 3300, label = "italic(p)*'CO'[2]*' (log'[10]*' µatm)'", 
             parse = TRUE, colour = "#c9d2d7", family = "Futura", size = 3.5, hjust = 0,
             label.size = 0) +
    scale_fill_manual(values = c("#5bb5b5", "#bc90c1")) +
    scale_x_continuous(breaks = pCO2_contour %$% seq(min(DIC_sim), max(DIC_sim), 300)) +
    scale_y_continuous(breaks = pCO2_contour %$% seq(min(TA_sim), max(TA_sim), 300)) +
    labs(x = expression("C"["T"]*" (µM)"), y = expression("A"["T"]*" (µM)")) +
    coord_cartesian(xlim = pCO2_contour %$% c(min(DIC_sim), max(DIC_sim)), 
                    ylim = pCO2_contour %$% c(min(TA_sim), max(TA_sim)),
                    expand = FALSE, clip = "off") +
    mytheme

Fig_S3 %>%
  ggsave(filename = "Fig_S3.pdf", path = "Figures",
         width = 22, height = 10, unit = "cm", device = cairo_pdf)

# 4.3 ΔTA and ΔDIC ####
DIC_TA %<>% 
  rename(TA = TA_corrected, DIC = DIC_calculated) %>%
  select(-c(pH_initial, Temp)) %>%
  filter(Species != "Initial") %>%
  mutate(Species = Species %>% fct_drop(),
         ID = ID %>% fct_drop()) %>%
  left_join( # join initial and blank/coralline final samples horizontally
    DIC_TA %>%
      filter(Species == "Initial") %>%
      rename(TA_initial = TA_corrected, DIC_initial = DIC_calculated) %>%
      select(starts_with("."), Round, Treatment, TA_initial, DIC_initial),
    by = c(".chain", ".iteration", ".draw", "Round", "Treatment")
  ) %>%
  mutate(delta_TA = (TA - TA_initial) / 2, # calculate ΔTA and ΔDIC h^-1
         delta_DIC = (DIC - DIC_initial) / 2) %>%
  select(-c(TA, TA_initial, DIC, DIC_initial))

DIC_TA

# 4.4 Volume ####
# See details in Oxygen.R under 4.2 Volume.
# 4.4.1 Prepare data ####
V <- incu_meta %>%
  filter(!is.na(Volume)) %>%
  mutate(Species = Species %>% fct_drop()) %>%
  select(ID, Species, Treatment, Volume)

# 4.4.2 Prior simulation ####
ggplot() +
  geom_density(aes(rgamma(1e5, 175^2 / 10^2, 175 / 10^2))) + # 10 mL seems like a reasonable sd
  theme_minimal() +
  theme(panel.grid = element_blank())

# 4.4.3 Run model ####
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
      select(Species, Treatment, Volume) %>%
      compose_data(),
    chains = 8,
    parallel_chains = parallel::detectCores(),
    iter_warmup = 1e4,
    iter_sampling = 1e4)

# 4.4.4 Model checks ####
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
# chains look fine

# 4.4.5 Prior-posterior comparison ####
# sample priors
V_prior <- prior_samples(
  model = V_mod,
  data = V %>% 
    select(Species, Treatment, Volume) %>%
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

# 4.4.6 Add volume to estimates ####
DIC_TA %<>%
  left_join( # match exact volumes
    V, by = c("ID", "Species", "Treatment")
  ) %>%
  left_join( # match estimated volumes
    V_samples %>%
      recover_types(V %>% select(Species, Treatment)) %>%
      spread_draws(V_mu[Species, Treatment]) %>%
      ungroup(),
    by = c("Species", "Treatment", ".chain", ".iteration", ".draw")
  ) %>% # merge Volume into one column
  mutate(Volume = coalesce(Volume, V_mu)) %>%
  select(-V_mu)

# 4.4.7 Convert volume to L ####
# The function seacarb::rho() needs temperature and salinity. The temperature measurement
# needs to be taken wen the water is weighed so I cannot use the titration temperature. The
# best source of this temperature is that measured by the oxygen meter and can be got from
# O2_estimates.
O2_estimates <- read.csv("O2_estimates.csv")

DIC_TA %<>%
  left_join(
    O2_estimates %>%
      select(starts_with("."), Round, T_nP, T_R) %>%
      pivot_longer(cols = c(T_nP, T_R),
                   values_to = "Temp", 
                   names_to = "Treatment") %>%
      mutate(Treatment = if_else(Treatment == "T_nP",
                                 "Light", "Dark") %>% 
                           fct()) %>%
      distinct(.chain, .iteration, .draw, Round, Treatment, .keep_all = TRUE),
    by = c(".chain", ".iteration", ".draw", "Round", "Treatment")
  ) %>%
  mutate(Volume_L = Volume / rho(Salinity, Temp)[1])

# 4.4.8 Correct for volume ####
# Multiplying ΔTA and ΔDIC by Volume_L converts them from µM h^-1 to µmol h^-1 because
# µM h^-1 is equivalent to µmol L^-1 h^-1, which * L leaves µmol h^-1.
DIC_TA %<>%
  mutate(delta_TA_µmol = delta_TA * Volume_L,
         delta_DIC_µmol = delta_DIC * Volume_L) %>%
  select(-c(delta_TA, delta_DIC, Volume, Volume_L))

# 4.4 Blank correction ####
DIC_TA %<>% 
  filter(Species != "Blank") %>%
  mutate(Species = Species %>% fct_drop(),
         ID = ID %>% fct_drop()) %>%
  left_join( # join blank and coralline samples horizontally
    DIC_TA %>%
      filter(Species == "Blank") %>%
      rename(delta_TA_µmol_blank = delta_TA_µmol, delta_DIC_µmol_blank = delta_DIC_µmol) %>%
      select(starts_with("."), Round, Treatment, delta_TA_µmol_blank, delta_DIC_µmol_blank),
    by = c(".chain", ".iteration", ".draw", "Round", "Treatment")
  ) %>%
  mutate(delta_TA_µmol_corrected = delta_TA_µmol - delta_TA_µmol_blank, # subtract blank
         delta_DIC_µmol_corrected = delta_DIC_µmol - delta_DIC_µmol_blank) %>%
  select(-c(delta_TA_µmol, delta_TA_µmol_blank, delta_DIC_µmol, delta_DIC_µmol_blank))

# 4.5 Mass ####






# 4.6 Calcification and CO2 fixation ####
G_P <- DIC_TA %>%
  select(starts_with("."), ID, Round, Species, Treatment, 
         delta_TA_corrected, delta_DIC_corrected) %>% # all incubations ran for 2 h, hence division by 2
  mutate(G_µM_h = -delta_TA_corrected / 2 / 2, # G is inversely related to TA and frees up two protons per CaCO3
         CO2_µM_h = -delta_DIC_corrected / 2 - G_µM_h) # some C is fixed by G instead of P and needs to be subtracted








