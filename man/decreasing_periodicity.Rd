% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decreasing_periodicity.R
\name{decreasing_periodicity}
\alias{decreasing_periodicity}
\title{Periodic decreasing temperature trend}
\usage{
decreasing_periodicity(
  y_ini = c(N = 400, N = 400, N = 400),
  temp_ini = rep(20, 3),
  temp_cmin = rep(18, 3),
  temp_cmax = c(25, 28, 35),
  ro = rep(0.7, 3),
  lambda = rep(5e-05, 3),
  A = rep(5, 3),
  B = rep(0.06, 3),
  m = rep(1/5, 3),
  time_start = 2005,
  time_end = 2100,
  leap = 1/12
)
}
\arguments{
\item{y_ini}{Initial population values (must be written with its name: N).}

\item{temp_ini}{Initial temperature.}

\item{temp_cmin}{Minimum critical temperature.}

\item{temp_cmax}{Maximum critical temperature.}

\item{ro}{Population growth rate at optimum temperature.}

\item{lambda}{Marginal loss by non-thermodependent intraspecific competition.}

\item{A}{Temperature wave amplitude.}

\item{B}{Parameter specifying the period of the trend (period is (2 pi)/|B|).}

\item{m}{Decreasing slope.}

\item{time_start}{Start of time sequence.}

\item{time_end}{End of time sequence.}

\item{leap}{Time sequence step.}
}
\value{
(1) A data.frame with columns having the simulated trends.

(2) A two-panel figure in which (a) shows the population abundance curves
           represented by solid lines and the corresponding carrying capacities
           are represented by shaded areas. In (b) the temperature trend is shown.
           The three simultaneous simulations are depicted by different colors, i.e.
           1st brown, 2nd green and 3rd blue.
}
\description{
This function simulates the effect of a decrease temperature trend
             with periodic variability on the abundance of ectotherm populations.
}
\details{
Three populations and/or scenarios can be simulated simultaneously.
        The temperature trend is determined by an expression that shows a
        linear decrease with periodic variability. The amplitude, amplitude,
        period and downward speed of change can be specified. In each input
        vector, the parameters for the three simulations must be specified
        (finite numbers for the initial population abundance). The simulations
        are obtained by a model that incorporates the effects of temperature
        over time, which leads to a non-autonomous ODE approach. This is
        function uses the ODE solver implemented in the package deSolve
        (Soetaert et al., 2010).
}
\examples{

#######################################################################
  #Example 1: Different initial population abundances.
#######################################################################

decreasing_periodicity(y_ini = c(N = 100, N = 200, N = 400),
                      temp_ini = rep(30,3),
                      temp_cmin = rep(18,3),
                      temp_cmax = rep(35,3),
                      ro = rep(0.7,3),
                      lambda = rep(0.00005,3),
                      A = rep(3,3),
                      B = rep(0.6,3),
                      m = rep(1/5,3),
                      time_start = 2005,
                      time_end = 2100,
                      leap = 1/12)

#######################################################################
  #Example 2: Different thermal tolerance ranges.
#######################################################################

temp_cmin3 <- 18
temp_cmin2 <- 10/9*temp_cmin3
temp_cmin1 <- 10/9*temp_cmin2

temp_cmax1 <- 32.4
temp_cmax2 <- 10/9*temp_cmax1
temp_cmax3 <- 10/9*temp_cmax2

decreasing_periodicity(y_ini = c(N = 100, N = 100, N = 100),
                      temp_ini = rep(26,3),
                      temp_cmin = c(temp_cmin1,temp_cmin2,temp_cmin3),
                      temp_cmax = rep(temp_cmax1,temp_cmax2,temp_cmax3),
                      ro = rep(0.7,3),
                      lambda = rep(0.00005,3),
                      A = rep(2,3),
                      B = rep(0.6,3),
                      m = rep(1/5,3),
                      time_start = 2005,
                      time_end = 2100,
                      leap = 1/12)
\donttest{
#######################################################################
  #Example 3: Different relationships between initial environmental
  #           temperature and optimum temperature.
#######################################################################

temp_cmin <- 18
temp_cmax <- 35

# Temperature at which performance is at its maximum value.
temp_op <- (temp_cmax+temp_cmin)/3+sqrt(((temp_cmax+temp_cmin)/3)^2-
           (temp_cmax*temp_cmin)/3)

temp_ini1 <- (temp_cmin+temp_op)/2
temp_ini2 <- temp_op
temp_ini3 <- (temp_op+temp_cmax)/2

decreasing_periodicity(y_ini = c(N = 100, N = 100, N = 100),
                      temp_ini = c(temp_ini1,temp_ini2,temp_ini3),
                      temp_cmin = rep(temp_cmin,3),
                      temp_cmax = rep(temp_cmax,3),
                      ro = rep(0.7,3),
                      lambda = rep(0.00005,3),
                      A = rep(2,3),
                      B = rep(0.6,3),
                      m = rep(1/5,3),
                      time_start = 2005,
                      time_end = 2100,
                      leap = 1/12)

#######################################################################
  #Example 4: Different marginal losses by a non-thermodependent
  #           component of intraspecific competition.
#######################################################################

lambda3 <- 0.01
lambda2 <- 1/2*lambda3
lambda1 <- 1/2*lambda2

decreasing_periodicity(y_ini = c(N = 100, N = 100, N = 100),
                      temp_ini = rep(30,3),
                      temp_cmin = rep(18,3),
                      temp_cmax = rep(35,3),
                      ro = rep(0.7,3),
                      lambda = c(lambda1,lambda2,lambda3),
                      A = rep(2,3),
                      B = rep(0.6,3),
                      m = rep(1/5,3),
                      time_start = 2005,
                      time_end = 2100,
                      leap = 1/12)

#######################################################################
  #Example 5: Different wave amplitude.
#######################################################################

A3 <- 2
A2 <- 1/2 * A3
A1 <- 1/2 * A2

decreasing_periodicity(y_ini = c(N = 100, N = 100, N = 100),
                      temp_ini = rep(30,3),
                      temp_cmin = rep(18,3),
                      temp_cmax = rep(35,3),
                      ro = rep(0.7,3),
                      lambda = rep(0.00005,3),
                      A = c(A1,A2,A3),
                      B = rep(0.6,3),
                      m = rep(1/5,3),
                      time_start = 2005,
                      time_end = 2100,
                      leap = 1/12)

#######################################################################
  #Example 6: Different period.
#######################################################################

B3 <- pi/5
B2 <- 1/2 * B3
B1 <- 1/2 * B2

decreasing_periodicity(y_ini = c(N = 100, N = 100, N = 100),
                      temp_ini = rep(30,3),
                      temp_cmin = rep(18,3),
                      temp_cmax = rep(35,3),
                      ro = rep(0.7,3),
                      lambda = rep(0.00005,3),
                      A = rep(2,3),
                      B = c(B1,B2,B3),
                      m = rep(1/5,3),
                      time_start = 2005,
                      time_end = 2100,
                      leap = 1/12)
}
}
\references{
Soetaert, K., Petzoldt, T., & Setzer, R. (2010). Solving Differential Equations
           in R: Package deSolve. Journal of Statistical Software, 33(9), 1 - 25.
           doi:http://dx.doi.org/10.18637/jss.v033.i09
}
