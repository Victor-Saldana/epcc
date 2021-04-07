#'Projected values under IPCC RCP2.6 scenario
#'
#'@description This function allows obtaining the projected increment in environmental
#'             temperature according to the IPCC RCP2.6 scenario (2014).
#'
#'@param date A specific year or a vector of years.
#'
#'@details The temperature increment projection of the change in global mean surface
#'         temperature according to the IPCC RCP2.6 scenario. It is possible to get
#'         the value for one or various years.
#'
#'@references IPCC. (2014): Cambio climático 2014: Informe de síntesis. Contribución de los Grupos de trabajo
#'            I, II y III al Quinto Informe de Evaluación del Grupo Intergubernamental de Expertos sobre el
#'            Cambio Climático [Equipo principal de redacción, R.K. Pachauri y L.A. Meyer (eds.)]. IPCC, Ginebra,
#'            Suiza, 157 págs.
#'
#'@examples
#'
#'########################################################################
#'  #Example 1: Projection of the temperature increase for a given year.
#'########################################################################
#'
#'date <- 2050
#'temp <- get_RCP2.6(date)
#'temp
#'
#'########################################################################
#'  #Example 2: Projection of the temperature increase for a vector of years.
#'########################################################################
#'
#'date <- seq(2005,2100,1/12)
#'temp <- get_RCP2.6(date)
#'plot(date,temp, type="l")

get_RCP2.6  <- function(date){
  ifelse(( 2005 <= date & date<=2065),0.61*(100/61)^{(date-2005)/60},ifelse((2065<date & date<=2100),1,NA))
}
