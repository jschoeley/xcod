# XCOD: Estimate expected and excess deaths by cause using coherent compositional regression

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13353995.svg)](https://doi.org/10.5281/zenodo.13353995)

[Jonas Schöley](https://orcid.org/0000-0002-3340-8518)

With XCOD you can estimate time series of expected and excess deaths by cause. After fitting XCOD to the data and learning the cause-specific trends and seasonalities as well as residual dependencies across causes you can use XCOD to extrapolate the trends and simulate time series of cause specific expected deaths. Various excess death statistics (P-score, absolute excess deaths, mortality ratio) can be calculated along with associated prediction intervals and p-values.

## Model specification

A vector of death counts by cause $j = 1, \ldots, k$ at time $t$, $\mathbf{D}^{(j)}_t = (D_t^1, \ldots, D_t^k)$, is modeled as the product of total expected deaths, $\lambda_t$, and a vector of cause-specific expected proportions on all deaths, $\mathbf{p}^{(j)}_t = (p^1_t, \ldots, p^k_t)$,

$$
\mathbf{D}^{(j)}_t = \lambda_t\mathbf{p}^{(j)}_t \times \exp \left( \mathbf{U}_t^{(j)}\right),
$$

with  $\mathbf{U}^{(j)}_t = (U^1_t, \ldots, U^k_t)$ being a vector of multivariate-normal distributed residuals, $\mathbf{U}_t \sim \mathcal{N}_k(0, \mathbf{\Sigma}_k)$, capturing the dependencies across causes of death. The exponent and multiplication are element wise.

In order to achieve robust and speedy estimates the data is fitted in a three-step procedure.

First, the expected total deaths at time t, $\lambda_t$, are fitted via Poisson regression as a log-linear function of time, featuring both a long-term trend and a cyclical seasonal component:

$$
\lambda_t = \exp \left( \beta_0 + g(t; \boldsymbol{\beta}_t) + s(t; T, \boldsymbol{\beta}_s) \right),
$$

where $g(t; \boldsymbol{\beta}_t)$ is a function of time modeling the long term trend, and $s(t; \boldsymbol{\beta}_s)$ is a cyclical function of time $t$ with period $T$ to capture seasonality.

Second, compositional regressions are fitted to estimate the expected proportion of cause-specific deaths on all deaths. Specifically, the cause-specific proportions, transformed via the centered-log-ratio, are expressed as a linear function of time,

$$
\mathrm{clr}_j(p^{j}_t) = \gamma_0^j + g(t; \boldsymbol{\gamma}_t) + s(t; T, \boldsymbol{\gamma}_s^j),
$$

with $p^{j}_t$ being the $j^\textrm{th}$ element of $\mathbf{p}_t$ and $\mathrm{clr}_j = \left( \frac {p^j_t} {G(\mathbf{p}_t)} \right)$, where $G(\mathbf{p}_t) = \sqrt[k]{\prod^k p^j_t}$ is the geometric mean over the proportions at a given time. The regressions are fitted via least-squares to the observed proportions.

Third, the residual covariance across causes of death, $\mathbf{\Sigma}_k$, is estimated from the within-sample model residuals, $\mathbf{U}_t = \log D_t - \log ( \hat\lambda_t\mathbf{\hat p}^{(j)}_t )$.

Simulated expected deaths by cause are drawn from the quantized predictive distribution of the fitted model:

$$
\mathbf{D}^{\mathrm{sim}(j)}_t = \left\lfloor \hat\lambda_t\mathbf{\hat p}^{(j)}_t \times \exp \left( \mathbf{U}_t^{\mathrm{sim}(j)}\right) \right\rfloor,
$$

where

$$
\mathbf{U}_t^{\mathrm{sim}(j)} \sim \mathcal{N}_k(0, \mathbf{\hat\Sigma}_k).
$$

These samples of expected deaths are the basis for samples of excess death estimates. We define the P-score for cause $j$ at time $t$ as the percent difference between the observed deaths and the deaths simulated under the expected distribution:

$$
\mathrm{Pscore}^{\mathrm{sim},j}_t = \frac {D^{\mathrm{obs},j}_t - D^{\mathrm{sim},j}_t} {D^{\mathrm{sim},j}_t}.
$$

The cumulative P-score over the interval $(a, b)$ is defined as

$$
\mathrm{Pscore}^{\mathrm{sim},j}_{(a,b)} = \frac {\sum_{t=a}^b (D^{\mathrm{obs},j}_t - D^{\mathrm{sim},j}_t)} {\sum_{t=a}^b D^{\mathrm{sim},j}_t}.
$$

Likewise the (cumulative) absolute excess deaths for cause $j$ at time $t$ are simulated as

$$
\mathrm{Excess}^{\mathrm{sim},j}_t = D^{\mathrm{obs},j}_t - D^{\mathrm{sim},j}_t,
$$

and

$$
\mathrm{Excess}^{\mathrm{sim},j}_{(a,b)} = \sum_{t=a}^b (D^{\mathrm{obs},j}_t - D^{\mathrm{sim},j}_t).
$$