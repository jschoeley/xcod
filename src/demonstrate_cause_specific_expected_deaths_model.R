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

# Load XCOD functions ---------------------------------------------

source('src/XCOD.R')

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
  rnbinom(nrow(sim$X_df), mu = sim$causeD$pred$lambda, size = 5)

sim$causeE <- list()
sim$causeE$para <-
  c(b0 = 0, b1 = 0.003, b2 = 0.2, b3 = log(1.01),
    b4 = log(3.11), b5 = log(3.55), b6 = log(1.001), b7 = 0)
sim$causeE$pred <- data.frame(
  sim$X_df, lambda = exp(sim$X%*%sim$causeE$para)
)
sim$causeE$pred$deaths <-
  rpois(nrow(sim$X_df), sim$causeE$pred$lambda)

# Sum causes to total number of deaths ----------------------------

sim$total <-
  cbind(
    sim$X_df,
    deathsA = sim$causeA$pred$deaths,
    deathsB = sim$causeB$pred$deaths,
    deathsC = sim$causeC$pred$deaths,
    deathsD = sim$causeD$pred$deaths,
    deathsE = sim$causeE$pred$deaths
  )
sim$total$deathsTotal <-
  sim$total$deathsA + sim$total$deathsB + sim$total$deathsC +
  sim$total$deathsD + sim$total$deathsE

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
    pD = sim$total$deathsD / sim$total$deathsTotal,
    pE = sim$total$deathsE / sim$total$deathsTotal
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
expected$strata <- unique(the_analysis_data$stratum)

# fit separately by stratum
expected$agesex <- map(
  expected$strata, ~{
    stratum_subset <- filter(the_analysis_data, stratum == .x)
    XCOD(
      df = stratum_subset,
      formula_total = "origin_time + as.factor(seasonal_time)",
      formula_prop_dense = "origin_time + as.factor(seasonal_time)",
      formula_prop_sparse = "1",
      cols_prop = c('pA', 'pB', 'pC', 'pD', '-pE'),
      col_total = 'deathsTotal',
      col_origin_time = 't',
      col_seasonal_time = 'w',
      col_cvflag = 'cv_flag',
      nsim = 100,
      basis = 'ilr'
    ) %>%
      mutate(stratum = .x)
  }
)
# merge strata
expected$agesex <- bind_rows(expected$agesex)

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
  name_parts = c('pA', 'pB', 'pC', 'pD', 'pE'),
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
  name_parts = c('pA', 'pB', 'pC', 'pD', 'pE'),
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

# p-scores for causes A & B combined into AB, as well as C & D & E
GetExcessByCause(
  xcod_out = expected$agesex,
  name_parts = list(AB = c('pA', 'pB'), C = 'pC', D = 'pD', E = 'pE'),
  measure = 'absolute',
  quantiles = c(0.1, 0.5, 0.9)
)

# test if causes ABCDE combined give same results as ALLCAUSE
GetExcessByCause(
  xcod_out = expected$agesex,
  name_parts = list(ABCDE = c('pA', 'pB', 'pC', 'pD', 'pE'), ALLCAUSE = 'ALLCAUSE'),
  measure = 'absolute',
  quantiles = c(0.1, 0.5, 0.9), cumulative = TRUE
)

# Visualize training and predicted proportions --------------------

observed_vs_forecasted <- list()
observed_vs_forecasted$data <- list()
observed_vs_forecasted$data <-
  bind_rows(
    expected = GetExcessByCause(
      expected$total, name_parts = c('pA', 'pB', 'pC', 'pD', 'pE'),
      measure = 'expected'
    ),
    observed = GetExcessByCause(
      expected$total, name_parts = c('pA', 'pB', 'pC', 'pD', 'pE'),
      measure = 'observed'
    ),
    .id = 'measure'
  ) %>%
  pivot_longer(cols = starts_with('Q')) %>%
  separate(col = name, into = c('quantile', 'part'), sep = '_')

observed_vs_forecasted$data %>%
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

# Visualize observed vs. expected ---------------------------------

observed_vs_expected <- list()
observed_vs_expected$data <- list()
observed_vs_expected$data <-
  bind_rows(
    expected = GetExcessByCause(
      expected$total, name_parts = c('pA', 'pB', 'pC', 'pD', 'pE'),
      measure = 'expected'
    ),
    observed = GetExcessByCause(
      expected$total, name_parts = c('pA', 'pB', 'pC', 'pD', 'pE'),
      measure = 'observed'
    ),
    .id = 'measure'
  ) %>%
  pivot_longer(cols = starts_with('Q')) %>%
  separate(col = name, into = c('quantile', 'part'), sep = '_') %>%
  pivot_wider(names_from = quantile, values_from = value)

observed_vs_expected$data %>%
  ggplot() +
  geom_ribbon(
    aes(x = origin_time, ymin = Q025, ymax = Q975), fill = 'grey90',
    data = . %>% filter(measure == 'expected')
  ) +
  geom_point(
    aes(x = origin_time, y = Q500), color = 'grey50',
    # plot observed counts by cause over training period
    data = . %>% filter(measure == 'observed')
  ) +
  geom_line(
    aes(x = origin_time, y = Q500), color = 'black',
    # plot average expected counts by cause over test period
    data = . %>% filter(measure == 'expected')
  ) +
  facet_wrap(~part, ncol = 2, scales = 'free_y') +
  scale_fill_brewer(type = 'div', palette = 9) +
  theme_minimal() +
  labs(
    title = 'Observed and forecasted deaths by cause',
    fill = 'Cause of death',
    y = 'Deaths',
    x = 'Months since 2015'
  )

# Visualize P-scores ----------------------------------------------

expected$total %>%
  GetExcessByCause(
    name_parts = c('pA', 'pB', 'pC', 'pD', 'pE'),
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
