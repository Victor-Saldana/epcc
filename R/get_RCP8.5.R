#'Projected values under IPCC RCP8.5 scenario
#'
#'@description This function allows obtaining the projected increment in environmental
#'             temperature according to the IPCC RCP8.5 scenario.
#'
#'@param date A specific year or a vector of years.
#'@param a,b Parameters obtained from the adjustment of the temperature trend.
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
#'@examples
#'
#'########################################################################
#'  #Example 1: Projection of the temperature increase for a given year.
#'########################################################################
#'
#'date <- 2050
#'temp <- get_RCP8.5(date,a=exp(coef(m)[1]), b=coef(m)[2])
#'temp
#'
#'########################################################################
#'   #Example 2: Projection of the temperature increase for a vector of years.
#'########################################################################
#'
#'date <- seq(2005,2100,1/12)
#'temp <- get_RCP8.5(date,a=exp(coef(m)[1]), b=coef(m)[2])
#'plot(date,temp,type="l")
#'

get_RCP8.5 <- function(date,a,b) {a * exp(b * date)}
values <- c(0.61, 2, 3.7)
x<- c(2005,2065,2100)
y<- values
df <- data.frame(x, y)

m<- nls(y ~ exp(loga + b * x), df, start = list( loga = log(2), b = 0.005),control = list (maxiter = 500))
y_est<-predict(m,df$x)


