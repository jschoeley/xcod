# Misc ------------------------------------------------------------

# Ensure that the xcod object is sorted
EnsureSortedXCOD <- function (xcod_out) {
  xcod_out[order(xcod_out[['stratum']],
                 xcod_out[['origin_time']],
                 xcod_out[['seasonal_time']], decreasing = FALSE),]
}

# Return data frame of row-wise quantiles over columns of X
Rowquantiles <- function (X, prob, type = 4, na.rm = TRUE) {
  t(apply(X, 1, quantile, prob = prob, type = type, na.rm = na.rm, names = FALSE))
}

# Expected estimation ---------------------------------------------

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
  Dk_obs <- round(df_prediction[,cols_prop]*df_prediction[['OBS_ALLCAUSE']],
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

# Excess derivation -----------------------------------------------

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