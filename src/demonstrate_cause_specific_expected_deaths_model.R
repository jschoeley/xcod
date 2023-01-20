#' Estimate Weekly Expected Deaths by Cause over Strata
#' Given Time-series of Total Deaths and Cause-Proportions
#'
#' First data is simulated and split into training-test chunks.
#' Then a GAM is fit to the training time series of total death
#' counts and total deaths are predicted for the test period.
#' The cause specific proportions are transformed by the additive
#' log-ratio transform and then modeled using a GAM, predicted
#' for the test period, and back transformed to proportions.
#' Given the expected total deaths and expected cause-specific
#' proportions for the test period, one can derive expected
#' deaths by cause over the test period.

# Init ------------------------------------------------------------

library(shiny)
library(tidyverse)

# Setup data grid -------------------------------------------------

sim <- list()
sim$dims <- list(
  # t: months since origin
  t   = 1:(12*3),
  # w: months since start of year
  w   = 1:12,
  # sex
  sex = relevel(as.factor(c('F', 'M')), ref = 'F'),
  # only three age groups for demonstration,
  # real analysis should have more age strata
  age = relevel(as.factor(c('0-40', '40-80', '80+')), ref = '0-40')
)
# create an empty frame of time and strata without data
sim$X_df <-
  expand.grid(
    # months since study origin
    t = sim$dims$t,
    # discrete sex
    sex = sim$dims$sex,
    # discrete age
    age = sim$dims$age
  )
# months into year
sim$X_df$w <- sim$X_df$t%%12+1
# we use the following model to specify the expected number of deaths
# for a given cause by time, month of year, sex, and age group
# lambda_t = exp( b0 + b1*t + b2*sin(2*pi*w/12) + b3*sex + b4*sex*t +
#                b5*age1 + b6*age2 + b7*age3 )
sim$X <- model.matrix(~t + sin(2*pi*w/12) + sex + sex*t + age, data = sim$X_df)

# Visualize simulation parameter choices --------------------------

ui <- fluidPage(
  
  sidebarLayout(
    
    sidebarPanel(
      sliderInput('b0', 'Intercept', min = 0, max = 12, value = 0),
      sliderInput('b1', 'Time slope', min = -0.01, max = 0.01, value = 0, step = 0.001),
      sliderInput('b2', 'Seasonality', min = -0.1, max = 0.1, value = 0, step = 0.0001),
      sliderInput('b3', 'Rate ratio Males to Females', min = 1, max = 1.2, value = 1, step = 0.01),
      sliderInput('b4', 'Rate ratio 40-80 to 0:40', min = 1, max = 10, value = 1, step = 0.01),
      sliderInput('b5', 'Rate ratio 80+ to 0:40', min = 1, max = 10, value = 1, step = 0.01),
      sliderInput('b6', 'exp(Sex*Time) Interaction', min = 0.9, max = 1.1, value = 1, step = 0.001)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

server <- function(input, output) {
  
  output$distPlot <- renderPlot({
    betas <- c(input$b0, input$b1, input$b2, log(input$b3), log(input$b4), log(input$b5), log(input$b6))
    lambda <- exp(sim$X%*%betas)
    y <- rpois(length(lambda), lambda)
    expected_deaths_X <- data.frame(sim$X_df, lambda = lambda, y = y)
    ggplot(expected_deaths_X) +
      aes(x = t, y = lambda) +
      geom_line() +
      geom_point(aes(y = y)) +
      facet_grid(age~sex, scales = 'free_y')
  })
}

# start shiny app
#shinyApp(ui, server)

# Simulate expected deaths ----------------------------------------

# simulate expected and observed deaths for three causes of deaths
sim$causeA <- list()
sim$causeA$para <-
  c(b0 = 4, b1 = 0.01, b2 = 0.1, b3 = log(0.95),
    b4 = log(4.87), b5 = log(8.44), b6 = log(0.98))
sim$causeA$pred <- data.frame(
  sim$X_df, lambda = exp(sim$X%*%sim$causeA$para)
)
sim$causeA$pred$deaths <-
  rpois(nrow(sim$X_df), sim$causeA$pred$lambda)

sim$causeB <- list()
sim$causeB$para <-
  c(b0 = 6, b1 = -0.003, b2 = 0.001, b3 = log(1.17),
    b4 = log(4.14), b5 = log(7.66), b6 = log(0.99))
sim$causeB$pred <- data.frame(
  sim$X_df, lambda = exp(sim$X%*%sim$causeB$para)
)
sim$causeB$pred$deaths <-
  rpois(nrow(sim$X_df), sim$causeB$pred$lambda)

sim$causeC <- list()
sim$causeC$para <-
  c(b0 = 6, b1 = 0.003, b2 = -0.2, b3 = log(1.01),
    b4 = log(3.11), b5 = log(3.55), b6 = log(1.001))
sim$causeC$pred <- data.frame(
  sim$X_df, lambda = exp(sim$X%*%sim$causeC$para)
)
sim$causeC$pred$deaths <-
  rpois(nrow(sim$X_df), sim$causeC$pred$lambda)

sim$causeD <- list()
sim$causeD$para <-
  c(b0 = 6, b1 = 0.003, b2 = 0.2, b3 = log(1.01),
    b4 = log(3.11), b5 = log(3.55), b6 = log(1.001))
sim$causeD$pred <- data.frame(
  sim$X_df, lambda = exp(sim$X%*%sim$causeD$para)
)
sim$causeD$pred$deaths <-
  rpois(nrow(sim$X_df), sim$causeD$pred$lambda)


# Sum causes to total number of deaths ----------------------------

sim$total <-
  cbind(
    sim$X_df,
    deathsA = sim$causeA$pred$deaths,
    deathsB = sim$causeB$pred$deaths,
    deathsC = sim$causeC$pred$deaths,
    deathsD = sim$causeD$pred$deaths
  )
sim$total$deathsTotal <-
  sim$total$deathsA + sim$total$deathsB + sim$total$deathsC +
  sim$total$deathsD

# Derive weekly death proportions by cause ------------------------

sim$prop <-
  cbind(
    sim$X_df,
    # total deaths in this week, sex, and age
    deathsTotal = sim$total$deathsTotal,
    # proportion of deaths due to cause A on all deaths
    # in given week, sex, and age
    pA = sim$total$deathsA / sim$total$deathsTotal,
    pB = sim$total$deathsB / sim$total$deathsTotal,
    pC = sim$total$deathsC / sim$total$deathsTotal,
    pD = sim$total$deathsD / sim$total$deathsTotal
  )

# The dataset as we would use it in the analysis ------------------

the_analysis_data <-
  sim$prop
the_analysis_data$cv_flag <- ifelse(
  the_analysis_data$t <= 24, 'training', 'test'
)
the_analysis_data$stratum <-
  interaction(the_analysis_data$sex, the_analysis_data$age)

# Predict expected deaths by cause --------------------------------

#' Predict Expected Deaths By Cause
#'
#' @param df
#'   A data frame.
#' @param cols_prop
#'   Character vector of column names giving the weekly death
#'   proportions by cause.
#' @param col_total
#'   Quoted column name for total deaths.
#' @param col_stratum 
#'   Quoted column name for stratum.
#' @param col_origin_time
#'   Quoted column name for numeric time since origin
#'   (e.g. months since Jan 2015).
#' @param col_seasonal_time
#'   Quoted column name for numeric seasonal time
#'   (e.g. months into year).
#' @param col_cvflag
#'   Quoted column name for cross-validation flag column. Column
#'   must be character with "training" for data used to fit the models
#'   and test for time points to make predictions for.
#' @param nsim
#'   Number of simulation draws. Default = 1000.
#' @param quantiles
#'   Numeric vector of quantiles to report for predicted distribution.
#' @param basis
#'   The basis for the coda transformation of the data. Can be "ilr"
#'   (default), "alr", or "cdp". See ?coda.base::coordinates.
XCOD <- function (
  df,
  cols_prop, col_total, col_stratum, col_origin_time,
  col_seasonal_time, col_cvflag,
  nsim = 1000, quantiles = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975),
  basis = 'ilr'
) {

  ## preparation
    
  require(mgcv)
  require(coda.base)
  
  N = nrow(df)
  p = length(cols_prop)-1
  
  df_training <- df[df$cv_flag == 'training',]
  df_prediction <- df
  
  ## model expected totals
  
  deathsTotal_form <- as.formula(paste0(
    col_total, "~",
    # log-linear time trend separate by stratum
    col_origin_time, "*", col_stratum, "+",
    # cyclical spline for seasonality separate
    # by stratum
    "s(",
    col_seasonal_time,
    ", bs = 'cp', by = ", col_stratum, ")"
  ))
  deathsTotal_fit <- gam(
    deathsTotal_form,
    data = df_training,
    family = poisson(link = 'log')
  )
  
  ## simulate from expected totals posterior predictive dist
  
  deathsTotal_Xprd <-
    predict(deathsTotal_fit, newdata = df_prediction, type = 'lpmatrix')
  deathsTotal_beta <- coef(deathsTotal_fit)
  deathsTotal_beta_sim <- MASS::mvrnorm(
    nsim, deathsTotal_beta,
    vcov(deathsTotal_fit, freq = FALSE, unconditional = TRUE)
  )
  deathsTotal_lambda_sim <- c(apply(
    deathsTotal_beta_sim, 1,
    FUN = function (b) exp(deathsTotal_Xprd%*%b)
  ))
  L <- matrix(deathsTotal_lambda_sim, nrow = N, ncol = nsim)
  L <- cbind(exp(deathsTotal_Xprd%*%deathsTotal_beta), L)

  ## model proportions
  
  P <- as.matrix(df_training[,cols_prop])
  Plr <- coordinates(P, basis = basis)
  Plr_hat <- array(NA, dim = list(i = N, j = nsim+1, p = p))
  
  props_form <- update.formula(deathsTotal_form, plr ~ .)
  for (k in 1:p) {
    the_data <- cbind(df_training, plr = Plr[,k])
    exclude_zero_props <- is.infinite(the_data$plr)
    the_data <- the_data[!exclude_zero_props,]
    prop_gam_fit <- gam(props_form,
                        family = gaussian(link = 'identity'),
                        data = the_data
    )
    prop_gam_Xprd <-
      predict(prop_gam_fit, newdata = df_prediction, type = 'lpmatrix')
    prop_gam_beta <- coef(prop_gam_fit)
    prop_gam_beta_sim <- MASS::mvrnorm(
      nsim, prop_gam_beta,
      vcov(prop_gam_fit, freq = FALSE, unconditional = TRUE)
    )
    prop_gam_plr_sim <- c(apply(
      prop_gam_beta_sim, 1,
      FUN = function (b) prop_gam_Xprd%*%b
    ))
    Plr_hat[,-1,k] <-
      matrix(prop_gam_plr_sim, nrow = N, ncol = nsim)
    Plr_hat[,1,k] <- prop_gam_Xprd%*%prop_gam_beta
  }
  
  P_hat <- array(NA, dim = list(i = N, j = nsim+1, p = p+1))
  for (j in 1:nsim) {
    P_hat_j <- composition(Plr_hat[,j,], basis = basis)
    for (k in 1:(p+1))
      P_hat[,j,k] <- P_hat_j[,k]
  }
  
  ## derive statistics of interest
  
  # Expected deaths by cause (mean + simulations)
  Dk_hat <- array(NA, dim = list(i = N, j = nsim+1, p = p+1))
  for (k in 1:(p+1)) {
    Dk_hat[,-1,k] <-
      apply(P_hat[,-1,k]*L[,-1], 2, function (lambda_k) rpois(N, lambda_k))
    Dk_hat[,1,k] <- P_hat[,1,k]*L[,1]
  }
  # quantiles
  Dk_hat_quantiles <-
    aperm(
      apply(Dk_hat[,-1,], c(1,3), quantile, probs = quantiles,
            na.rm = TRUE, type = 1),
      c(2,1,3)
    )
  
  # Observed deaths by cause
  Dk_obs <- round(df[,cols_prop]*df[,col_total], digits = 0)

  # Excess deaths by cause (mean + simulations)
  Dk_xcs <- array(NA, dim = list(i = N, j = nsim+1, p = p+1))
  for (k in 1:(p+1)) {
    Dk_xcs[,,k] <- Dk_obs[,k]-Dk_hat[,,k]
  }
  # quantiles
  Dk_xcs_quantiles <-
    aperm(
      apply(Dk_xcs[,-1,], c(1,3), quantile, probs = quantiles,
            na.rm = TRUE, type = 1),
      c(2,1,3)
    )
  
  # P-scores by cause (mean + simulations)
  Dk_psc <- array(NA, dim = list(i = N, j = nsim+1, p = p+1))
  for (k in 1:(p+1)) {
    # note: infinities are possible here if expected deaths are 0
    Dk_psc[,,k] <- Dk_xcs[,,k]/Dk_hat[,,k]*100
  }
  # quantiles
  Dk_psc_quantiles <-
    aperm(
      apply(Dk_psc[,-1,], c(1,3), quantile, probs = quantiles,
            na.rm = TRUE, type = 7),
      c(2,1,3)
    )
  
  # bind predictions and derived stats + quantiles to input data frame
  for (k in 1:(p+1)) {
    X <- cbind(
      # observed
      Dk_obs[,k],
      # expected
      Dk_hat[,1,k], Dk_hat_quantiles[,,k],
      # excess
      Dk_xcs[,1,k], Dk_xcs_quantiles[,,k],
      # pscore
      Dk_psc[,1,k], Dk_psc_quantiles[,,k]
    )
    colnames(X) <- c(
      # observed
      paste0('OBS_', cols_prop[k]),
      # expected
      paste0('XPC_', cols_prop[k], '_AVG'), paste0('XPC_', cols_prop[k], '_Q', quantiles),
      # excess
      paste0('XCS_', cols_prop[k], '_AVG'), paste0('XCS_', cols_prop[k], '_Q', quantiles),
      # p-score
      paste0('PSC_', cols_prop[k], '_AVG'), paste0('PSC_', cols_prop[k], '_Q', quantiles)
    )
    df <- cbind(df, X)
  }
  
  return(df)
}

df <- XCOD(
  df = the_analysis_data,
  cols_prop = c('pA', 'pB', 'pC', 'pD'),
  col_total = 'deathsTotal',
  col_stratum = 'stratum',
  col_origin_time = 't',
  col_seasonal_time = 'w',
  col_cvflag = 'cv_flag',
  nsim = 1000,
  quantiles = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975),
  basis = 'ilr'
)

# Convert to long format ------------------------------------------

df_long <-
  df %>%
  pivot_longer(cols = starts_with(c('OBS', 'XPC', 'XCS', 'PSC'))) %>%
  separate(col = name, into = c('stat', 'part', 'quantile'), sep = '_')

# Visualize training and predicted proportions --------------------

df_long %>%
  ggplot() +
  geom_area(
    aes(x = t, y = value, fill = part),
    # plot observed counts by cause over training period
    data = . %>% filter(stat == 'OBS', cv_flag == 'training')
  ) +
  geom_area(
    aes(x = t, y = value, fill = part),
    # plot average expected counts by cause over test period
    data = . %>% filter(stat == 'XPC', cv_flag == 'test', quantile == 'AVG')
  ) +
  facet_grid(age ~ sex, scales = 'free_y') +
  scale_fill_brewer(type = 'div', palette = 9) +
  theme_minimal() +
  labs(
    title = 'Observed and forecasted deaths by cause',
    fill = 'Cause of death',
    y = 'Deaths',
    x = 'Months since 2015'
  )

ggsave('cover.png', path = './ass/',
       device = ragg::agg_png, width = 170, height = 150, units = 'mm')

# Visualize P-scores ----------------------------------------------

# should be mostly non-significant in this case
df %>%
  filter(cv_flag == 'test') %>%
  ggplot(aes(x = t)) +
  geom_ribbon(
    aes(ymin = PSC_pB_Q0.025, ymax = PSC_pB_Q0.975),
    color = NA, fill = 'grey80') +
  geom_hline(yintercept = 0) +
  geom_line(
    aes(y = PSC_pB_AVG), color = 'red'
  ) +
  facet_grid(age ~ sex, scales = 'free_y') +
  theme_minimal() +
  labs(
    title = 'Percent excess for cause B',
    y = 'Monthly P-score',
    x = 'Months since 2015'
  )

# Observed versus predicted versus true mean ----------------------

ggplot(df) +
  aes(x = t) +
  geom_line(aes(y = lambda), color = 'grey', data = sim$causeA$pred) +
  geom_point(aes(y = OBS_pA, color = cv_flag)) +
  geom_ribbon(aes(ymin = XPC_pA_Q0.025, ymax = XPC_pA_Q0.975, fill = cv_flag),
              color = NA, alpha = 0.1) +
  geom_line(aes(y = XPC_pA_AVG, color = cv_flag), data = df) +
  facet_grid(age ~ sex, scales = 'free_y')

ggplot(df) +
  aes(x = t) +
  geom_line(aes(y = lambda), color = 'grey', data = sim$causeB$pred) +
  geom_point(aes(y = OBS_pB, color = cv_flag)) +
  geom_ribbon(aes(ymin = XPC_pB_Q0.025, ymax = XPC_pB_Q0.975, fill = cv_flag),
              color = NA, alpha = 0.1) +
  geom_line(aes(y = XPC_pB_AVG, color = cv_flag), data = df) +
  facet_grid(age ~ sex, scales = 'free_y')

ggplot(df) +
  aes(x = t) +
  geom_line(aes(y = lambda), color = 'grey', data = sim$causeC$pred) +
  geom_point(aes(y = OBS_pC, color = cv_flag)) +
  geom_ribbon(aes(ymin = XPC_pC_Q0.025, ymax = XPC_pC_Q0.975, fill = cv_flag),
              color = NA, alpha = 0.1) +
  geom_line(aes(y = XPC_pC_AVG, color = cv_flag), data = df) +
  facet_grid(age ~ sex, scales = 'free_y')
