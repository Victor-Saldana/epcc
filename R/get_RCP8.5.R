#'Projected values under IPCC RCP8.5 scenario
#'
#'@description This function allows obtaining the projected increment in environmental
#'             temperature according to the IPCC RCP8.5 scenario.
#'
#'@param date A specific year or a vector of years.
#'
#'@details The temperature increment projection of the change in global mean surface
#'         temperature according to the IPCC RCP8.5 scenario. It is possible to get
#'         the value for one or various years.
#'
#'@references IPCC. (2014): Cambio climático 2014: Informe de síntesis. Contribución de los Grupos de trabajo
#'            I, II y III al Quinto Informe de Evaluación del Grupo Intergubernamental de Expertos sobre el
#'            Cambio Climático [Equipo principal de redacción, R.K. Pachauri y L.A. Meyer (eds.)]. IPCC, Ginebra,
#'            Suiza, 157 págs.
#'
#'@export
#'@examples
#'
#'########################################################################
#'  #Example 1: Projection of the temperature increase for a given year.
#'########################################################################
#'
#'date <- 2050
#'temp <- get_RCP8.5(date)
#'temp
#'
#'########################################################################
#'   #Example 2: Projection of the temperature increase for a vector of years.
#'########################################################################
#'
#'date <- seq(2005,2100,1/12)
#'temp <- get_RCP8.5(date)
#'plot(date,temp,type="l")
#'



get_RCP8.5 <- function(date) {
  values <- c(0.61, 2, 3.7)
  q<- c(2005,2065,2100)
  p<- values
  df <- data.frame(q, p)

  m<- nls(p ~ exp(loga + b * q), df, start = list( loga = log(2), b = 0.005),control = list (maxiter = 500))
  y_est<-predict(m,df$q)

  k<- exp(coef(m)[1]) * exp(coef(m)[2] * date)
  return(k)
  }



