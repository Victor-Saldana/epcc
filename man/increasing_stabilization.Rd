% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/increasing_stabilization.R
\name{increasing_stabilization}
\alias{increasing_stabilization}
\title{Increasing temperature and stabilization}
\usage{
increasing_stabilization(
  y_ini = c(N = 400, N = 400, N = 400),
  temp_ini = rep(20, 3),
  temp_cmin = rep(18, 3),
  temp_cmax = c(25, 28, 32),
  ro = rep(0.7, 3),
  lambda = rep(5e-05, 3),
  temp_stabilization = rep(25, 3),
  q = rep(0.03, 3),
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

\item{temp_stabilization}{Stabilization temperature.}

\item{q}{Temparature increase factor.}

\item{time_start}{Start of time sequence.}

\item{time_end}{End of time sequence.}

\item{leap}{sequence increase.}
}
\value{
(1) A data.frame with columns having the simulated trends.

(2) A two-panel figure in which (a) shows the population abundance curves represented
           by solid lines and the corresponding carrying capacities are represented by shaded
           areas. In (b) the temperature trend is shown. The three simultaneous simulations are
           depicted by different colors, i.e. 1st brown, 2nd green and 3rd blue.
}
\description{
This function allows simulating the effect of an increase in environmental temperature,
             which stabilizes at a specific value, on the abundance of ectotherm populations.
}
\details{
Three populations and/or scenarios can be simulated simultaneously. A logistic
        type function determines the temperature trend. The temperature increases and
        then stabilizes at a specific value. In each input vector, the parameters for
        the three simulations must be specified (finite numbers for the initial population
        abundance). The simulations are obtained by a model that incorporates the effects
        of temperature over time, which leads to a non-autonomous ODE approach. This is
        function uses the ODE solver implemented in the package deSolve
        (Soetaert et al., 2010).
}
\examples{

#######################################################################
  #Example 1: Different initial population abundances.
#######################################################################

increasing_stabilization(y_ini = c(N = 100, N = 200, N = 400),
                        temp_ini = rep(26,3),
                        temp_cmin = rep(18,3),
                        temp_cmax = rep(30,3),
                        ro = rep(0.7,3),
                        lambda = rep(0.00005,3),
                        temp_stabilization = rep(30,3),
                        q = rep(0.03,3),
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

increasing_stabilization(y_ini = c(N = 100, N = 100, N = 100),
                        temp_ini = rep(30,3),
                        temp_cmin = c(temp_cmin1,temp_cmin2,temp_cmin3),
                        temp_cmax = c(temp_cmax1,temp_cmax2,temp_cmax3),
                        ro = rep(0.7,3),
                        lambda = rep(0.00005,3),
                        temp_stabilization = rep(33,3),
                        q = rep(0.03,3),
                        time_start = 2005,
                        time_end = 2100,
                        leap = 1/12)
\donttest{
#######################################################################
  #Example 3: Different relationships between initial environmental
  #           temperature and optimum temperature.
#######################################################################

temp_cmin <- 18
temp_cmax <- 30

# Temperature at which performance is at its maximum value.
temp_op <- (temp_cmax+temp_cmin)/3+sqrt(((temp_cmax+temp_cmin)/3)^2-
           (temp_cmax*temp_cmin)/3)

temp_ini1 <- (temp_cmin+temp_op)/2
temp_ini2 <- temp_op
temp_ini3 <- (temp_op+temp_cmax)/2

increasing_stabilization(y_ini = c(N = 100, N = 100, N = 100),
                        temp_ini = c(temp_ini1,temp_ini2,temp_ini3),
                        temp_cmin = rep(temp_cmin,3),
                        temp_cmax = rep(temp_cmax,3),
                        ro = rep(0.7,3),
                        lambda = rep(0.00005,3),
                        temp_stabilization = rep(32,3),
                        q = rep(0.03,3),
                        time_start = 2005,
                        time_end = 2100,
                        leap = 1/12)

#######################################################################
  #Example 4: Different stabilitations temperature.
#######################################################################

temp_stabilization3 <- 42
temp_stabilization2 <- 14/13*temp_stabilization3
temp_stabilization1 <- 14/13*temp_stabilization2

increasing_stabilization(y_ini = c(N = 100, N = 100, N = 100),
                        temp_ini = rep(30,3),
                        temp_cmin = rep(18,3),
                        temp_cmax = rep(40,3),
                        ro = rep(0.7,3),
                        lambda = rep(0.00005,3),
                        temp_stabilization = c(temp_stabilization1,
                                               temp_stabilization2,
                                               temp_stabilization3),
                        q = rep(0.03,3),
                        time_start = 2005,
                        time_end = 2100,
                        leap = 1/12)

#######################################################################
  #Example 5: Different marginal losses by a non-thermodependent
  #           component of intraspecific competition.
#######################################################################

lambda3 <- 0.01
lambda2 <- 1/2*lambda3
lambda1 <- 1/2*lambda2

increasing_stabilization(y_ini = c(N = 100, N = 100, N = 100),
                        temp_ini = rep(30,3),
                        temp_cmin = rep(18,3),
                        temp_cmax = rep(40,3),
                        ro = rep(0.7,3),
                        lambda = c(lambda1,lambda2,lambda3),
                        temp_stabilization = rep(35,3),
                        q = rep(0.03,3),
                        time_start = 2005,
                        time_end = 2100,
                        leap = 1/12)
}
}
\references{
Soetaert, K., Petzoldt, T., & Setzer, R. (2010). Solving Differential Equations in R: Package
           deSolve. Journal of Statistical Software, 33(9), 1 - 25.
           doi:http://dx.doi.org/10.18637/jss.v033.i09
}
