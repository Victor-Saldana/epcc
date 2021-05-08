
# ${\huge e}_{\tiny PCC}$

<!-- badges: start -->
<!-- badges: end -->

The R package ${\Huge e}_{\tiny PCC}$ provides several functions that allow to simulate the effects of thermal sensitivity and the exposition to different environmental temperature trends on the abundance dynamics of ectothermic populations. More specifically, parameters associated with the optimum population performance ($r_{o}$) and critical thermal limits for survival ($CT_{min}$ and $CT_{max}$) can be specified. For instance, it is possible to simulate if the thermal optimum is below or above the current temperature and also to determine the potential outcome when considering that the population is constituted by a thermal specialist or generalist organisms (i.e., wider or narrower thermal limit ranges). Regarding environmental temperature, the package is intended to encompass a variety of situations, ranging from predicted scenarios proposed by the Intergovernmental Group of Experts on Climate Change (IPCC) at global and local specific areas, potential trends ranging from heating or cooling pulses, and trends with different temperature variability levels through time. Estos escenarios potenciales permiten simular tendencias disímiles que pueden ocurrir en diferentes latitudes y lapsos de tiempo. Además, los posibles efectos no térmicos intraespecíficos sobre la dinámica de la población (Svanback \& Bolnick, 2007; Rich et al., 2009) también pueden incorporarse mediante un parámetro específico (es decir, $\lambda$, la pérdida marginal por competencia).
The package also provides functions to assess the outcome of two common interspecific interactions, i.e., competence and predation, when the population growth of one of the species is affected by temperature (the prey population in the case of predation).  In addition, a function showing an ectotherm population with age structure is presented. These functions have been developed considering global warming trends as proposed by the IPCC (2014).
The package ${\Huge e}_{\tiny PCC}$ has been built upon a classical ordinary differential equation (ODE) solver, i.e., the  R package deSolve (Soetaert et al., 2010), but our approach involves incorporating temperature effects through time, which leads to a non-autonomous system ODE approach.
For each temperature trend, ${\Huge e}_{\tiny PCC}$ provides a specific function that allows visualizing variations in abundance, the corresponding carrying capacity, and temperature trends.

## Installation

You can install the released version of epcc from [GitHub](https://github.com/Victor-Saldana/epcc) with:

``` r
library(devtools)
install\_github(''Victor-Saldana/epcc'')
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(epcc)
## basic example code
```
## References
``` r
IPCC, 2014: Climate Change 2014: Synthesis Report. Contribution of Working Groups I, II and III to the Fifth Assessment Report of the Intergovernmental Panel on Climate Change [Core Writing Team, R.K. Pachauri and L.A. Meyer (eds.)]. IPCC, Geneva, Switzerland, 151 pp.

```
``` r
Rich,, H. B., Quinn, T. P., Scheuerell, M. D., \& Schindler, D. E. (2009). Climate and intraspecific competition control the growth and life history of juvenile sockeye salmon (Oncorhynchus nerka) in Iliamna Lake, Alaska. Canadian Journal of Fisheries and Aquatic Sciences, 66(2), 238-246.doi:10.1139/f08-210

```
``` r

Soetaert, K., Petzoldt, T., \& Setzer, R. (2010). Solving Differential Equations in R: Package deSolve. Journal of Statistical Software, 33(9), 1 - 25. doi:http://dx.doi.org/10.18637/jss.v033.i09
```

``` r
Svanback, R., \& Bolnick, D. I. (2007). Intraspecific competition drives increased resource use diversity within a natural population. Proceedings of the Royal Society B: Biological Sciences, 274(1611), 839-844. doi:10.1098/rspb.2006.0198 
```