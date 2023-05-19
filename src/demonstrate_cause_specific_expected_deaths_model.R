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

# Functions -------------------------------------------------------

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
    nsim = 10, basis = 'ilr'
) {
  
  ## preparation
  
  require(mgcv)      # for gam()
  require(coda.base) # for compositional data analysis operations
  
  N = nrow(df)
  p = length(cols_prop)-1
  
  idx_train <- which(df[,col_cvflag] == 'training')
  vec_stratum <- df[,col_stratum]
  # predict over whole data...
  df_prediction <- df[,c(
    col_stratum, col_origin_time, col_seasonal_time, col_cvflag,
    col_total, cols_prop
  )]
  # standardize names
  colnames(df_prediction) <- c('stratum', 'origin_time', 'seasonal_time',
                               'cv_flag', 'OBS_ALLCAUSE', cols_prop)
  # ...but train over this part of input data:
  df_training <- df_prediction[idx_train,]
  
  ## model expected total deaths over time
  
  deathsTotal_form <- as.formula(paste0(
    "OBS_ALLCAUSE~",
    # log-linear time trend separate by stratum
    "origin_time*stratum",
    # cyclical spline for seasonality separate
    # by stratum
    "+s(seasonal_time, bs = 'cp', by = stratum)"
  ))
  deathsTotal_fit <- gam(
    deathsTotal_form,
    data = df_training,
    family = poisson(link = 'log')
  )
  
  ## simulate from expected totals posterior predictive dist
  
  # design matrix
  deathsTotal_Xprd <-
    predict(deathsTotal_fit, newdata = df_prediction, type = 'lpmatrix')
  # coefficients
  deathsTotal_beta <- coef(deathsTotal_fit)
  # simulate coefficients from MV-Normal distribution
  deathsTotal_beta_sim <- MASS::mvrnorm(
    nsim, deathsTotal_beta,
    vcov(deathsTotal_fit, freq = FALSE, unconditional = TRUE)
  )
  # derive simulated replications of linear predictor
  # (rate of total deaths) 
  deathsTotal_lambda_sim <- c(apply(
    deathsTotal_beta_sim, 1,
    FUN = function (b) exp(deathsTotal_Xprd%*%b)
  ))
  L <- matrix(deathsTotal_lambda_sim, nrow = N, ncol = nsim)
  # add mean linear predictor as first column to simulations
  L <- cbind(exp(deathsTotal_Xprd%*%deathsTotal_beta), L)
  
  ## model expected proportions of deaths by cause over time
  
  # transform proportions to log-ratio analysis space
  P <- as.matrix(df_training[,cols_prop])
  Plr <- coordinates(P, basis = basis)
  Plr_hat <- array(NA, dim = list(i = N, j = nsim+1, p = p))
  
  # model and extrapolate expected proportions in log-ratio space
  # separate by cause...
  props_form <- update.formula(deathsTotal_form, plr ~ .)
  for (k in 1:p) {
    the_data <- cbind(df_training, plr = Plr[,k])
    # zero proportions are excluded from fitting
    exclude_zero_props <- is.infinite(the_data$plr)
    the_data <- the_data[!exclude_zero_props,]
    prop_gam_fit <- gam(props_form,
                        family = gaussian(link = 'identity'),
                        data = the_data)
    # simulate predicted proportions from fitted model
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
  
  # convert predicted proportions from log-ratio to proportion space
  P_hat <- array(NA, dim = list(i = N, j = nsim+1, p = p+1))
  for (j in 1:(nsim+1)) {
    P_hat_j <- composition(Plr_hat[,j,], basis = basis)
    for (k in 1:(p+1))
      P_hat[,j,k] <- P_hat_j[,k]
  }
  
  ## simulate expected deaths by cause
  
  # observed deaths by cause
  Dk_obs <- round(df_prediction[,cols_prop]*df_prediction[,'OBS_ALLCAUSE'],
                  digits = 0)
  
  # expected deaths by cause (mean + simulations)
  Dk_hat <- array(NA, dim = list(i = N, j = nsim+1, p = p+1))
  for (k in 1:(p+1)) {
    Dk_hat[,-1,k] <-
      apply(P_hat[,-1,k]*L[,-1], 2, function (lambda_k) rpois(N, lambda_k))
    Dk_hat[,1,k] <- P_hat[,1,k]*L[,1]
  }
  
  # bind cause-specific predictions and simulations to input data
  for (k in 1:(p+1)) {
    X <- cbind(
      # observed
      Dk_obs[,k],
      # expected average & simulated
      Dk_hat[,,k]
    )
    colnames(X) <- c(
      # observed
      paste0('OBS_', cols_prop[k]),
      # expected average & simulated
      paste0('XPC_AVG_', cols_prop[k]),
      paste0('XPC_SIM', 1:nsim, '_', cols_prop[k])
    )
    df_prediction <- cbind(df_prediction, X)
  }
  # sum cause-specific predictions and simulations to totals and
  # bind to input data
  X <- cbind(
    # expected average & simulated
    apply(Dk_hat, 1:2, sum)
  )
  colnames(X) <- c(
    # expected average & simulated
    paste0('XPC_AVG_ALLCAUSE'),
    paste0('XPC_SIM', 1:nsim, '_ALLCAUSE')
  )
  df_prediction <- cbind(df_prediction, X)
  
  return(df_prediction)
}

EnsureSortedXCOD <- function (xcod_out) {
  xcod_out[order(xcod_out[['stratum']],
                 xcod_out[['origin_time']],
                 xcod_out[['seasonal_time']], decreasing = FALSE),]
}

# return data frame of row-wise quantiles over columns of X
Rowquantiles <- function (X, prob, type = 4, na.rm = TRUE) {
  t(apply(X, 1, quantile, prob = prob, type = type, na.rm = na.rm, names = FALSE))
}

GetExcessByCause <- function (
    xcod_out, name_parts,
    measure = 'absolute',
    quantiles = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975),
    cumulative = FALSE, origin_time_start_of_cumulation = 0
) {
  
  # forgiveness please, this info should really be explicit in xcod_out
  nsim = max(as.integer(
    sub('(^.+_SIM)([[:digit:]]+)(.*$)','\\2',
        grep('[[:digit:]]',names(xcod_out), value = TRUE))
  ))
  nrow = NROW(xcod_out)
  xcod_sorted <- EnsureSortedXCOD(xcod_out)
  
  # data columns
  Y <- xcod_sorted[,grepl('XPC_|OBS_',colnames(xcod_sorted))]
  # label columns
  X <- xcod_sorted[,c('stratum', 'origin_time', 'seasonal_time', 'cv_flag')]
  
  # aggregate parts if requested
  if (is.list(name_parts)) {
    amalgamation_names <- names(name_parts)
    amalgamation_parts <- name_parts
    Y <- Map(function (x, y) {
      name <- x
      parts <- unlist(y)
      n_parts <- length(parts)
      
      # OBS
      OBS_colnames <- grep(paste0('OBS_', parts, collapse = '|'),
                           colnames(Y), value = TRUE)
      OBS_parts <- matrix(unlist(Y[,OBS_colnames]),
                          nrow = nrow, ncol = n_parts)
      OBS_amalgamation <- rowSums(OBS_parts)
      # XPC AVG
      XPC_AVG_colnames <- grep(paste0('XPC_AVG_', parts, collapse = '|'),
                               colnames(Y), value = TRUE)
      XPC_AVG_parts <- matrix(unlist(Y[,XPC_AVG_colnames]),
                              nrow = nrow, ncol = n_parts)
      XPC_AVG_amalgamation <- rowSums(XPC_AVG_parts)
      # XPC SIM
      XPC_SIM_colnames <- grep(paste0('XPC_SIM.+_', parts, collapse = '|'),
                               colnames(Y), value = TRUE)
      XPC_SIM_parts <- array(
        unlist(Y[,XPC_SIM_colnames]),
        dim = c(nrow, nsim, n_parts)
      )
      XPC_SIM_amalgamation <- apply(XPC_SIM_parts, 1:2, sum)
      # XCOD amalgamation
      xcod_amalgamation <- cbind(
        OBS_amalgamation, XPC_AVG_amalgamation, XPC_SIM_amalgamation
      )
      colnames(xcod_amalgamation) <-
        c(paste0(c('OBS_', 'XPC_AVG_'), name),
          paste0('XPC_SIM', 1:100, '_', name))
      
      return(xcod_amalgamation)
    },
    amalgamation_names, amalgamation_parts)
    
    # reorder columns
    Y <- do.call('cbind', Y)
    Y <- Y[,c(which(grepl('OBS',colnames(Y))),
              which(grepl('XPC_AVG',colnames(Y))),
              which(grepl('XPC_SIM',colnames(Y))))]
  } else {
    amalgamation_names <- name_parts
  }
  
  # accumulate data columns if requested
  if (isTRUE(cumulative)) {
    # set data to 0 for time points prior to accumulation start
    # so that we only accumulate from the accumulation start
    vec_timeselect <- X[['origin_time']] < origin_time_start_of_cumulation
    Y[vec_timeselect,] <- 0
    # this restarts the accumulation of a data vector whenever the
    # stratum vector changes in value
    vec_stratum <- X[['stratum']]
    Y <- apply(Y, 2, function (x) ave(x, vec_stratum, FUN = cumsum))
    # set values to NA prior to accumulation start so that derived
    # values become NA as well
    Y[vec_timeselect,] <- NA
  }

  for (part in amalgamation_names) {
    OBS <- Y[,paste0('OBS_', part)]
    j <- grepl(paste0('^XPC_SIM[[:digit:]]+_',part,'$'), colnames(Y))
    if (identical(measure, 'observed')) {
      MEASURE <- as.matrix(OBS)
    }
    if (identical(measure, 'expected')) {
      MEASURE <- Y[,j]
    }
    if (identical(measure, 'absolute')) {
      MEASURE <- apply(Y[,j], 2, function (XPC_SIM) {round(OBS-XPC_SIM,0)})
    }
    if (identical(measure, 'pscore')) {
      MEASURE <- apply(Y[,j], 2, function (XPC_SIM) {(OBS-XPC_SIM)/XPC_SIM*100})
    }
    if (identical(measure, 'ratio')) {
      MEASURE <- apply(Y[,j], 2, function (XPC_SIM) {OBS/XPC_SIM})
    }
    # get quantiles
    Q <- Rowquantiles(MEASURE, quantiles, type = 1)
    colnames(Q) <-
      paste0('Q', substr(
        formatC(quantiles, format = 'f', digits = 3),
        start = 3, stop = 5
      ), '_', part)
    X <- cbind(X,Q)
  }
  return(X)
}

# Setup data grid -------------------------------------------------

sim <- list()
sim$dims <- list(
  # t: months since origin
  t   = 0:(12*3-1),
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
  ) %>%
  mutate(
    excess_start = ifelse(t < 30, 0, 1)
  )
# months into year
sim$X_df$w <- sim$X_df$t%%12+1
# we use the following model to specify the expected number of deaths
# for a given cause by time, month of year, sex, and age group
# lambda_t = exp( b0 + b1*t + b2*sin(2*pi*w/12) + b3*sex + b4*age2 +
#                b5*age3 + b6*sex*t + b7*excess_start)
sim$X <- model.matrix(~1 + t + sin(2*pi*w/12) + sex + age + sex:t, data = sim$X_df)
sim$X <- cbind(sim$X, excess_start = sim$X_df$excess_start)

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
    b4 = log(4.87), b5 = log(8.44), b6 = log(0.98), b7 = 0)
sim$causeA$pred <- data.frame(
  sim$X_df, lambda = exp(sim$X%*%sim$causeA$para)
)
sim$causeA$pred$deaths <-
  rpois(nrow(sim$X_df), sim$causeA$pred$lambda)

sim$causeB <- list()
sim$causeB$para <-
  c(b0 = 6, b1 = -0.003, b2 = 0.001, b3 = log(1.17),
    b4 = log(4.14), b5 = log(7.66), b6 = log(0.99), b7 = 0)
sim$causeB$pred <- data.frame(
  sim$X_df, lambda = exp(sim$X%*%sim$causeB$para)
)
sim$causeB$pred$deaths <-
  rpois(nrow(sim$X_df), sim$causeB$pred$lambda)

sim$causeC <- list()
sim$causeC$para <-
  c(b0 = 6, b1 = 0.003, b2 = -0.2, b3 = log(1.01),
    b4 = log(3.11), b5 = log(3.55), b6 = log(1.001), b7 = 0.2)
sim$causeC$pred <- data.frame(
  sim$X_df, lambda = exp(sim$X%*%sim$causeC$para)
)
sim$causeC$pred$deaths <-
  rpois(nrow(sim$X_df), sim$causeC$pred$lambda)

sim$causeD <- list()
sim$causeD$para <-
  c(b0 = 6, b1 = 0.003, b2 = 0.2, b3 = log(1.01),
    b4 = log(3.11), b5 = log(3.55), b6 = log(1.001), b7 = 1)
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

expected <- list()

expected$agesex <- XCOD(
  df = the_analysis_data,
  cols_prop = c('pA', 'pB', 'pC', 'pD'),
  col_total = 'deathsTotal',
  col_stratum = 'stratum',
  col_origin_time = 't',
  col_seasonal_time = 'w',
  col_cvflag = 'cv_flag',
  nsim = 100,
  basis = 'ilr'
)

# Aggregate observed and expected over age-sex --------------------

library(data.table)

expected$agesex_dt <- as.data.table(expected$agesex)

expected$strata_cols <- c('stratum', 'origin_time', 'seasonal_time', 'cv_flag')
expected$value_cols <-
  c(
    grep('XPC_|OBS_', names(expected$agesex_dt), value = TRUE)
  )

expected$total <-
  groupingsets(
    expected$agesex_dt,
    j = lapply(.SD, sum),
    by = expected$strata_cols,
    sets = list(
      c('origin_time', 'seasonal_time', 'cv_flag')
    ),
    .SDcols = expected$value_cols
  ) %>%
  mutate(stratum = 'Total') %>%
  as.data.frame()

# Demonstrate derivation of excess statistics ---------------------

# pscores by cause
GetExcessByCause(
  # output of XCOD function
  xcod_out = expected$total,
  # parts of interest
  name_parts = c('ALLCAUSE'),
  # what measure to return
  measure = 'pscore',
  # quantiles of interest
  quantiles = c(0.1, 0.5, 0.9),
  cumulative = FALSE, origin_time_start_of_cumulation = 25
)

# cumulative pscores starting from time 30
GetExcessByCause(
  xcod_out = expected$agesex,
  name_parts = c('pA', 'pB', 'pC', 'pD'),
  measure = 'pscore',
  quantiles = c(0.1, 0.5, 0.9),
  # cumulative measures
  cumulative = TRUE,
  # accumulation starts here
  origin_time_start_of_cumulation = 30
)

# absolute cumulative number of excess by cause
# starting from time 25
GetExcessByCause(
  xcod_out = expected$agesex,
  name_parts = c('pA', 'pB', 'pC', 'pD'),
  measure = 'absolute',
  quantiles = c(0.1, 0.5, 0.9),
  # cumulative measures
  cumulative = TRUE,
  # accumulation starts here
  origin_time_start_of_cumulation = 25
)

# all cause p-scores
GetExcessByCause(
  xcod_out = expected$agesex,
  name_parts = c('ALLCAUSE'),
  measure = 'pscore',
  quantiles = c(0.1, 0.5, 0.9)
)

# p-scores for causes A & B combined into AB and C
GetExcessByCause(
  xcod_out = expected$agesex,
  name_parts = list(AB = c('pA', 'pB'), C = 'pC', D = 'pD'),
  measure = 'absolute',
  quantiles = c(0.1, 0.5, 0.9)
)

# test if causes ABCD combined give same results as ALLCAUSE
GetExcessByCause(
  xcod_out = expected$agesex,
  name_parts = list(ABCD = c('pA', 'pB', 'pC', 'pD'), ALLCAUSE = 'ALLCAUSE'),
  measure = 'absolute',
  quantiles = c(0.1, 0.5, 0.9), cumulative = TRUE
)

# Visualize training and predicted proportions --------------------

observed_vs_expected <- list()
observed_vs_expected$data <- list()
observed_vs_expected$data <-
  bind_rows(
    expected = GetExcessByCause(
      expected$total, name_parts = c('pA', 'pB', 'pC', 'pD'),
      measure = 'expected'
    ),
    observed = GetExcessByCause(
      expected$total, name_parts = c('pA', 'pB', 'pC', 'pD'),
      measure = 'observed'
    ),
    .id = 'measure'
  ) %>%
  pivot_longer(cols = starts_with('Q')) %>%
  separate(col = name, into = c('quantile', 'part'), sep = '_')

observed_vs_expected$data %>%
  ggplot() +
  geom_area(
    aes(x = origin_time, y = value, fill = part),
    # plot observed counts by cause over training period
    data = . %>% filter(cv_flag == 'training', measure == 'observed')
  ) +
  geom_area(
    aes(x = origin_time, y = value, fill = part),
    # plot average expected counts by cause over test period
    data = . %>% filter(cv_flag == 'test', measure == 'expected',
                        quantile == 'Q500')
  ) +
  facet_wrap(~ stratum, ncol = 2, scales = 'free_y') +
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

expected$total %>%
  GetExcessByCause(
    name_parts = c('pA', 'pB', 'pC', 'pD'),
    measure = 'pscore'
  ) %>%
  pivot_longer(cols = starts_with('Q')) %>%
  separate(col = name, into = c('quantile', 'part'), sep = '_') %>%
  pivot_wider(names_from = quantile, values_from = value) %>%
  filter(cv_flag == 'test') %>%
  ggplot(aes(x = origin_time)) +
  geom_ribbon(
    aes(ymin = Q025, ymax = Q975),
    color = NA, fill = 'grey80') +
  geom_hline(yintercept = 0) +
  geom_line(
    aes(y = Q500), color = 'red'
  ) +
  scale_x_continuous(breaks = 0:40) +
  facet_wrap(~ part) +
  theme_minimal() +
  labs(
    title = 'Percent excess by cause',
    y = 'Monthly P-score',
    x = 'Months since 2015'
  )
