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
#' @param formula_total RHS of formula passed to mgcv::gam describing
#'   time series of total deaths. Character string.
#' @param formula_prop_dense RHS of formula passed to mgcv::gam
#'   describing time series of cause of death proportions for parts with
#'   no or very little zero-shares. Character string.
#' @param formula_prop_sparse RHS of formula passed to mgcv::gam
#'   describing time series of cause of death proportions for parts with
#'   all or many zero-shares. Character string.
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
    formula_total = "s(origin_time, bs = 'tp', fx = TRUE, k = 3) + s(seasonal_time, bs = 'cp')",
    formula_prop_dense = "s(origin_time, bs = 'tp', fx = TRUE, k = 3) + s(seasonal_time, bs = 'cp')",
    formula_prop_sparse = "1",
    cols_prop, col_total, col_stratum = NULL, col_origin_time,
    col_seasonal_time, col_cvflag,
    nsim = 10, basis = 'ilr'
) {
  
  ## preparation
  
  require(mgcv)      # for gam()
  require(coda.base) # for compositional data analysis operations
  
  N = nrow(df)
  
  idx_cols_prop_sparse = grepl('^-', cols_prop)
  cols_prop_sparse = substr(cols_prop[idx_cols_prop_sparse], 2, 1000000L)
  cols_prop_dense = cols_prop[!idx_cols_prop_sparse]
  cols_prop_all = c(cols_prop_dense, cols_prop_sparse)
  # number of parts which are "sparse", i.e. mostly 0
  p_sparse = length(cols_prop_sparse)
  # number of parts which are not sparse
  p_dense = length(cols_prop_dense)
  # total number of parts
  p = p_sparse + p_dense
  
  idx_train <- which(df[,col_cvflag] == 'training')
  
  # prepare prediction data
  if (is.null(col_stratum)) {
    vec_stratum <- rep('All', N)
    # predict over whole data...
    df_prediction <- df[,c(
      col_origin_time, col_seasonal_time, col_cvflag,
      col_total, cols_prop_all
    )]
    df_prediction <- cbind(vec_stratum, df_prediction)
  } else {
    vec_stratum <- df[,col_stratum]
    # predict over whole data...
    df_prediction <- df[,c(
      col_stratum, col_origin_time, col_seasonal_time, col_cvflag,
      col_total, cols_prop_all
    )]
  }
  
  # standardize names
  colnames(df_prediction) <- c('stratum', 'origin_time', 'seasonal_time',
                               'cv_flag', 'OBS_ALLCAUSE', cols_prop_all)
  # ...but train over this part of input data:
  df_training <- df_prediction[idx_train,]
  
  ## model expected total deaths over time
  
  deathsTotal_form <- as.formula(paste0(
    "OBS_ALLCAUSE~",
    formula_total
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
  # derive simulated replications of Lambda
  # (rate of total deaths) 
  deathsTotal_lambda_sim <- c(apply(
    deathsTotal_beta_sim, 1,
    FUN = function (b) exp(deathsTotal_Xprd%*%b)
  ))
  L <- matrix(deathsTotal_lambda_sim, nrow = N, ncol = nsim)
  # add mean Lambda as first column to simulations
  L <- cbind(exp(deathsTotal_Xprd%*%deathsTotal_beta), L)
  
  ## model expected proportions of deaths by cause over time
  
  # predict proportions for sparse parts
  sparse <- list()
  if (p_sparse > 0) {
    # matrix holding training proportions of sparse parts
    sparse$P <- as.matrix(df_training[,cols_prop_sparse])
    sparse$props_form <- as.formula(paste0("p~", formula_prop_sparse))
    # matrix holding predicted proportions of sparse parts
    sparse$P_hat <- array(NA, dim = list(i = N, j = nsim+1, p = p_sparse))
    
    for (k in 1:(p_sparse)) {
      sparse$the_data <- cbind(
        df_training[,c('stratum', 'origin_time', 'seasonal_time')],
        p = sparse$P[,k]
      )
      sparse$prop_gam_fit <- gam(sparse$props_form,
                                 family = gaussian(link = 'identity'),
                                 data = sparse$the_data)
      # simulate predicted proportions from fitted model
      sparse$prop_gam_Xprd <-
        predict(sparse$prop_gam_fit, newdata = df_prediction, type = 'lpmatrix')
      sparse$prop_gam_beta <- coef(sparse$prop_gam_fit)
      sparse$prop_gam_beta_sim <- MASS::mvrnorm(
        nsim, sparse$prop_gam_beta,
        vcov(sparse$prop_gam_fit, freq = FALSE, unconditional = TRUE)
      )
      sparse$prop_gam_p_sim <- c(apply(
        sparse$prop_gam_beta_sim, 1,
        FUN = function (b) sparse$prop_gam_Xprd%*%b
      ))
      sparse$P_hat[,-1,k] <-
        matrix(sparse$prop_gam_p_sim, nrow = N, ncol = nsim)
      sparse$P_hat[,1,k] <- sparse$prop_gam_Xprd%*%sparse$prop_gam_beta
    }
    
    # closure
    for (i in 1:N) {
      for (j in nsim+1) {
        sparse$P_hat[i,j,]/sum(sparse$P_hat[i,j,])
      }
    }
    
    # clip negative predictions to 0
    sparse$P_hat[sparse$P_hat < 0] <- 0
    
    # share of sparse parts on all parts
    sparse$share_on_all <- apply(sparse$P_hat, 1:2, sum)
  }
  
  # predict proportions for dense parts
  dense <- list()
  if (p_dense > 0) {
    
    # matrix holding training proportions of dense parts
    dense$P <- as.matrix(df_training[,cols_prop_dense])
    # transform proportions to log-ratio analysis space
    dense$Plr <- coordinates(dense$P, basis = basis)
    
    # matrix holding predicted proportions of dense parts in log-ratio space
    dense$Plr_hat <- array(NA, dim = list(i = N, j = nsim+1, p = p_dense-1))
    
    # model and extrapolate expected proportions in log-ratio space
    # separate by cause...
    dense$props_form <- as.formula(paste0("plr~", formula_prop_dense))
    for (k in 1:(p_dense-1)) {
      dense$the_data <- cbind(df_training, plr = dense$Plr[,k])
      # zero proportions are excluded from fitting
      dense$exclude_zero_props <- is.infinite(dense$the_data$plr)
      dense$the_data <- dense$the_data[!dense$exclude_zero_props,]
      dense$prop_gam_fit <- gam(dense$props_form,
                                family = gaussian(link = 'identity'),
                                data = dense$the_data)
      # simulate predicted proportions from fitted model
      dense$prop_gam_Xprd <-
        predict(dense$prop_gam_fit, newdata = df_prediction, type = 'lpmatrix')
      dense$prop_gam_beta <- coef(dense$prop_gam_fit)
      dense$prop_gam_beta_sim <- MASS::mvrnorm(
        nsim, dense$prop_gam_beta,
        vcov(dense$prop_gam_fit, freq = FALSE, unconditional = TRUE)
      )
      dense$prop_gam_plr_sim <- c(apply(
        dense$prop_gam_beta_sim, 1,
        FUN = function (b) dense$prop_gam_Xprd%*%b
      ))
      dense$Plr_hat[,-1,k] <-
        matrix(dense$prop_gam_plr_sim, nrow = N, ncol = nsim)
      dense$Plr_hat[,1,k] <- dense$prop_gam_Xprd%*%dense$prop_gam_beta
    }
  }
  
  P_hat <- array(NA, dim = list(i = N, j = nsim+1, p = p))
  # if there are sparse parts of the expected composition, write those
  # into the P_hat composition array. Use the number of dense parts as
  # offset for the index as the dense parts should come first.
  if (p_sparse > 0) {
    for (k in 1:p_sparse) { P_hat[,,(k+p_dense)] <- sparse$P_hat[,,k] }
  }
  # if there are dense parts of the expected composition, write those
  # into the P_hat composition array, after converting them from
  # compositional to proportion space.
  if (p_dense > 0) {
    for (j in 1:(nsim+1)) {
      dense$P_hat_j <- composition(dense$Plr_hat[,j,], basis = basis)
      for (k in 1:p_dense) { P_hat[,j,k] <- dense$P_hat_j[,k] }
    }
  }
  # If the dense parts form a sub-composition, rescale this
  # sub-composition to its share on the total composition
  if ( p_dense > 0 & p_sparse > 0) {
    for (k in 1:p_dense) {
      P_hat[,,k] <- P_hat[,,k]*(1-sparse$share_on_all)
    }
  }
  
  ## simulate expected deaths by cause
  
  # observed deaths by cause
  Dk_obs <- round(df_prediction[,cols_prop_all]*df_prediction[['OBS_ALLCAUSE']],
                  digits = 0)
  
  # expected deaths by cause (mean + simulations)
  Dk_hat <- array(NA, dim = list(i = N, j = nsim+1, p = p))
  for (k in 1:p) {
    Dk_hat[,-1,k] <-
      apply(P_hat[,-1,k]*L[,-1], 2, function (lambda_k) rpois(N, lambda_k))
    Dk_hat[,1,k] <- P_hat[,1,k]*L[,1]
  }
  
  # bind cause-specific predictions and simulations to input data
  for (k in 1:p) {
    X <- cbind(
      # observed
      Dk_obs[,k],
      # expected average & simulated
      Dk_hat[,,k]
    )
    colnames(X) <- c(
      # observed
      paste0('OBS_', cols_prop_all[k]),
      # expected average & simulated
      paste0('XPC_AVG_', cols_prop_all[k]),
      paste0('XPC_SIM', 1:nsim, '_', cols_prop_all[k])
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
          paste0('XPC_SIM', 1:nsim, '_', name))
      
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