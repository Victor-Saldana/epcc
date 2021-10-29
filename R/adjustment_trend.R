#'Adjustment to observed population trends
#'
#' @description This function allows to fit the observed time series of abundance and ambient
#'              temperature trends together with the main descriptors of the TPCs for the
#'              salamander Desmognathus ochrophaeus.
#'
#'
#'@param y_ini Initial population for adjustment.
#'@param temp_cmin Minimum critical temperature for Desmognathus ochrophaeus (Layne & Claussen, 1987).
#'@param temp_cmax Maximum critical temperature for Desmognathus ochrophaeus (Markle & Kozak 2018).
#'@param roi Initial value for optimal intrinsic growth rate.
#'@param lambdai Initial value for the marginal loss of competition.
#'@param trend Acronym corresponding to the type of temperature trend used.
#'
#'
#'@details This function allows you to adjust the intrinsic growth rate of an ectotherm
#'         population by providing as input the values of environmental or body temperature
#'         together with the growth rate, using the nls2 function you can find the parameters
#'         temp_cmax, temp_cmax and ro, which are necessary to adjust the curve through a cubic polynomial
#'         which fulfills the essential conditions of a TPC.
#'
#'@export
#'@examples
#'
#'######################################################################
#'   #Example 1:Fitting observed time series of abundance and ambient
#'              temperature, using an increasing linear temperature trend.
#'#######################################################################
#'
#'adjustment_trend <- function(y_ini = c(N = 15),
#'                             temp_cmin = 4,
#'                             temp_cmax = 33.1,
#'                             roi = 0.3,
#'                             lambdai = 0.05,
#'                             trend = 1)
#'
#'######################################################################
#'   #Example 2:Fitting the observed time series of abundance and ambient
#'              temperature, using an periodic increasing temperature trend.
#'#######################################################################
#'
#'adjustment_trend <- function(y_ini = c(N = 15),
#'                             temp_cmin = 4,
#'                             temp_cmax = 33.1,
#'                             roi = 0.3,
#'                             lambdai = 0.05,
#'                             trend = 3)
#'
#'######################################################################
#'   #Example 3:Fitting the observed time series of abundance and ambient
#'              temperature, using a constant temperature trend.
#'#######################################################################
#'
#'adjustment_trend <- function(y_ini = c(N = 15),
#'                             temp_cmin = 4,
#'                             temp_cmax = 33.1,
#'                             roi = 0.3,
#'                             lambdai = 0.05,
#'                             trend = 5)
#'
#'
#'@references Rezende, E. L., & Bozinovic, F. (2019). Thermal performance across levels of biological
#'            organization. Philosophical Transactions of the Royal Society B: Biological Sciences,
#'            374(1778), 20180549.doi:10.1098/rstb.2018.0549
#'
#'@references Saldaña-Núñez, V.N., Córdova-Lepe, F.D. & Moreno-Gómez, F.N. (2019). Population dynamics in the face of
#'            climate change: Analysis of a cubic thermal performance curve in ectotherms. J. Phys.: Conf.
#'            Ser. 1329 012007.   doi:10.1088/1742-6596/1329/1/012007
#'

adjustment_trend <- function(y_ini = c(N = 15),
                             temp_cmin = 4,
                             temp_cmax = 33.1,
                             roi = 0.3,
                             lambdai = 0.05,
                             trend = 1){

  if(trend == 1) {

    ##############################################################################
    #Temperatura promedio, creciente lineal
    x1<- TAVG$Temperature.DATE.31.45.
    y1<- (TAVG$Temperature.TAVG.31.45.- 32)*5/9

    df1 <- data.frame(x1, y1)
    time_start<-x1[1]
    f1 <- function (times,temp_ini,m,time_start) {
      T <- temp_ini+m*(times-time_start)
    }

    m1<- nls2(y1 ~ f1(x1,temp_ini,m,time_start), df1, start = list( temp_ini=11.5,m=1/6),control = list (maxiter = 500))
    y_est1<-predict(m1,df1$x1)

    times<-seq(x1[1],x1[length(x1)])
    temp1=f1(times,coef(m1)[1],coef(m1)[2],time_start)

    dat1<-data.frame('x'=times,'y'=temp1)
    dat2<-data.frame('x'=x1,'y'=y1)
    data<-rbind(dat1,dat2)
    p1<- ggplot(data, aes(x=.data$x, y=.data$y)) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_line(data =subset(dat1,times>x1[1]  & times<x1[length(x1)]), color = "black",size=rel(1.1))+
      geom_point(data =subset(dat2,x1>x1[1]  & x1<x1[length(x1)]), color = "black",size=rel(2.5))+
      labs(x = "Time (years)",y="Temperature")+
      coord_cartesian(xlim = c(x1[1],x1[length(x1)]))+
      theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
      theme(axis.title.x = element_text(size = rel(1.2), angle = 00))


    # ABUNDANCIA
    BIO3=data.frame(BIO$YEAR[31:45],BIO$ABUNDANCE[31:45])

    xa<- BIO3$BIO.YEAR.31.45.
    ya<- BIO3$BIO.ABUNDANCE.31.45.

    dfa <- data.frame(xa, ya)

    temp_ini<-coef(m1)[1]

    fa <- function (times,y_ini,temp_cmin,temp_ini, temp_cmax,ro,lambda) {

      parms1<-c(temp_cmin,temp_ini, temp_cmax,ro,lambda)

      model1 <- function (times,y,parms1) {
        with(as.list(c(y)), {

          T <- f1(times,coef(m1)[1],coef(m1)[2],time_start)
          temp_op1<- ( temp_cmax+temp_cmin)/3+sqrt((( temp_cmax+temp_cmin)/3)^2-
                                                     ( temp_cmax*temp_cmin)/3)
          r1<- ro*((T-temp_cmin)*(temp_cmax-T)*T)/((temp_op1-temp_cmin)*(temp_cmax-temp_op1)*temp_op1)
          dN <-   r1 * N * (1 - lambda*(N / r1))
          list(dN,T,r1) })
      }

      out1 <- ode(y=y_ini, times, model1 ,parms1,method = "ode45")

      return(out1[,2])
    }

    ma<- nls2(ya ~ fa(xa,y_ini,temp_cmin,temp_ini, temp_cmax,ro,lambda), dfa, alg = "plinear", start = list(ro=0.7, lambda=0.0005),control = nls.control (maxiter = 500))
    y_esta<-predict(ma,dfa$xa)

    times<-seq(x1[1],x1[length(x1)])
    Abund=fa(times,y_ini,temp_cmin,temp_ini, temp_cmax,coef(ma)[1],coef(ma)[2])
    da1<-data.frame('x'=times,'y'=Abund )
    da2<-data.frame('x'=xa,'y'=ya )

    data<-rbind(da1,da2)

    p2<-ggplot(data, aes(x=.data$x, y=.data$y)) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_line(data =subset(da1,times>xa[1]  & times<xa[length(xa)]), color = "black",size=rel(1.1))+
      geom_point(data =subset(da2,xa>xa[1]  & xa<xa[length(xa)]), color = "black",size=rel(2))+
      labs(x = "Time (years)",y="Abundance")+
      coord_cartesian(xlim = c(xa[1],xa[length(xa)]))+
      theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
      theme(axis.title.x = element_text(size = rel(1.2), angle = 00))


    plot_grid(p1, p2)

    ##############################################################################
  } else if(trend== 2) {

    #Temperatura promedio, decreciente lineal
    x1<- TAVG$Temperature.DATE.31.45.
    y1<- (TAVG$Temperature.TAVG.31.45.- 32)*5/9

    df1 <- data.frame(x1, y1)
    time_start<-x1[1]
    f1 <- function (times,temp_ini,m,time_start) {
      T <- temp_ini-m*(times-time_start)
    }

    m1<- nls2(y1 ~ f1(x1,temp_ini,m,time_start), df1, start = list( temp_ini=13,m=1/6),control = list (maxiter = 500))
    y_est1<-predict(m1,df1$x1)

    times<-seq(x1[1],x1[length(x1)])
    temp1=f1(times,coef(m1)[1],coef(m1)[2],time_start)

    dat1<-data.frame('x'=times,'y'=temp1)
    dat2<-data.frame('x'=x1,'y'=y1)
    data<-rbind(dat1,dat2)
    p1<- ggplot(data, aes(x=.data$x, y=.data$y)) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_line(data =subset(dat1,times>x1[1]  & times<x1[length(x1)]), color = "black",size=rel(1.1))+
      geom_point(data =subset(dat2,x1>x1[1]  & x1<x1[length(x1)]), color = "black",size=rel(2.5))+
      labs(x = "Time (years)",y="Temperature")+
      coord_cartesian(xlim = c(x1[1],x1[length(x1)]))+
      theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
      theme(axis.title.x = element_text(size = rel(1.2), angle = 00))


    # ABUNDANCIA
    BIO3=data.frame(BIO$YEAR[31:45],BIO$ABUNDANCE[31:45])

    xa<- BIO3$BIO.YEAR.31.45.
    ya<- BIO3$BIO.ABUNDANCE.31.45.

    dfa <- data.frame(xa, ya)

    temp_ini<-coef(m1)[1]

    fa <- function (times,y_ini,temp_cmin,temp_ini, temp_cmax,ro,lambda) {

      parms1<-c(temp_cmin,temp_ini, temp_cmax,ro,lambda)

      model1 <- function (times,y,parms1) {
        with(as.list(c(y)), {

          T <- f1(times,coef(m1)[1],coef(m1)[2],time_start)
          temp_op1<- ( temp_cmax+temp_cmin)/3+sqrt((( temp_cmax+temp_cmin)/3)^2-
                                                     ( temp_cmax*temp_cmin)/3)
          r1<- ro*((T-temp_cmin)*(temp_cmax-T)*T)/((temp_op1-temp_cmin)*(temp_cmax-temp_op1)*temp_op1)
          dN <-   r1 * N * (1 - lambda*(N / r1))
          list(dN,T,r1) })
      }

      out1 <- ode(y=y_ini, times, model1 ,parms1,method = "ode45")

      return(out1[,2])
    }

    ma<- nls2(ya ~ fa(xa,y_ini,temp_cmin,temp_ini, temp_cmax,ro,lambda), dfa, alg = "plinear", start = list(ro=0.7, lambda=0.0005),control = nls.control (maxiter = 500))
    y_esta<-predict(ma,dfa$xa)

    times<-seq(x1[1],x1[length(x1)])
    Abund=fa(times,y_ini,temp_cmin,temp_ini, temp_cmax,coef(ma)[1],coef(ma)[2])
    da1<-data.frame('x'=times,'y'=Abund )
    da2<-data.frame('x'=xa,'y'=ya )

    data<-rbind(da1,da2)

    p2<-ggplot(data, aes(x=.data$x, y=.data$y)) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_line(data =subset(da1,times>xa[1]  & times<xa[length(xa)]), color = "black",size=rel(1.1))+
      geom_point(data =subset(da2,xa>xa[1]  & xa<xa[length(xa)]), color = "black",size=rel(2))+
      labs(x = "Time (years)",y="Abundance")+
      coord_cartesian(xlim = c(xa[1],xa[length(xa)]))+
      theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
      theme(axis.title.x = element_text(size = rel(1.2), angle = 00))


    plot_grid(p1, p2)

  } else if(trend== 3) {


    #Temperatura promedio, creciente con periocidad
    x1<- TAVG$Temperature.DATE.31.45.
    y1<- (TAVG$Temperature.TAVG.31.45.- 32)*5/9

    df1 <- data.frame(x1, y1)

    f1 <- function (times,temp_ini,mp, A,B) {
      T <- temp_ini+A*sin(B*(times-x1[1]))+ mp*(times-x1[1])
    }

    m1<- nls2(y1 ~ f1(x1,temp_ini,mp, A,B), df1, start = list( temp_ini=11.5,mp=1/6, A=1,B=0.2),control = list (maxiter = 500))
    y_est1<-predict(m1,df1$x1)

    times<-seq(x1[1],x1[length(x1)])
    temp1=f1(times,coef(m1)[1],coef(m1)[2],coef(m1)[3],coef(m1)[4])

    dat1<-data.frame('x'=times,'y'=temp1)
    dat2<-data.frame('x'=x1,'y'=y1)
    data<-rbind(dat1,dat2)
    p1<- ggplot(data, aes(x=.data$x, y=.data$y)) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_line(data =subset(dat1,times>x1[1]  & times<x1[length(x1)]), color = "black",size=rel(1.1))+
      geom_point(data =subset(dat2,x1>x1[1]  & x1<x1[length(x1)]), color = "black",size=rel(2.5))+
      labs(x = "Time (years)",y="Temperature")+
      coord_cartesian(xlim = c(x1[1],x1[length(x1)]))+
      theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
      theme(axis.title.x = element_text(size = rel(1.2), angle = 00))

    # ABUNDANCIA
    BIO3=data.frame(BIO$YEAR[31:45],BIO$ABUNDANCE[31:45])

    xa<- BIO3$BIO.YEAR.31.45.
    ya<- BIO3$BIO.ABUNDANCE.31.45.

    dfa <- data.frame(xa, ya)

    temp_ini<-coef(m1)[1]

    fa <- function (times,y_ini,temp_cmin,temp_ini, temp_cmax,ro,lambda) {

      parms1<-c(temp_cmin,temp_ini, temp_cmax,ro,lambda)

      model1 <- function (times,y,parms1) {
        with(as.list(c(y)), {

          T <- f1(times,coef(m1)[1],coef(m1)[2],coef(m1)[3],coef(m1)[4])
          temp_op1<- ( temp_cmax+temp_cmin)/3+sqrt((( temp_cmax+temp_cmin)/3)^2-
                                                     ( temp_cmax*temp_cmin)/3)
          r1<- ro*((T-temp_cmin)*(temp_cmax-T)*T)/((temp_op1-temp_cmin)*(temp_cmax-temp_op1)*temp_op1)
          dN <-   r1 * N * (1 - lambda*(N / r1))
          list(dN,T,r1) })
      }

      out1 <- ode(y=y_ini, times, model1 ,parms1,method = "ode45")

      return(out1[,2])
    }

    ma<- nls2(ya ~ fa(xa,y_ini,temp_cmin,temp_ini, temp_cmax,ro,lambda), dfa, alg = "plinear", start = list(ro=0.7, lambda=0.0005),control = nls.control (maxiter = 500))
    y_esta<-predict(ma,dfa$xa)

    times<-seq(x1[1],x1[length(x1)])
    Abund=fa(times,y_ini,temp_cmin,temp_ini, temp_cmax,coef(ma)[1],coef(ma)[2])
    da1<-data.frame('x'=times,'y'=Abund )
    da2<-data.frame('x'=xa,'y'=ya )

    data<-rbind(da1,da2)

    p2<-ggplot(data, aes(x=.data$x, y=.data$y)) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_line(data =subset(da1,times>xa[1]  & times<xa[length(xa)]), color = "black",size=rel(1.1))+
      geom_point(data =subset(da2,xa>xa[1]  & xa<xa[length(xa)]), color = "black",size=rel(2))+
      labs(x = "Time (years)",y="Abundance")+
      coord_cartesian(xlim = c(xa[1],xa[length(xa)]))+
      theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
      theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
    plot_grid(p1, p2)
  } else if(trend== 4) {

    #Temperatura promedio, decreciente con periocidad
    x1<- TAVG$Temperature.DATE.31.45.
    y1<- (TAVG$Temperature.TAVG.31.45.- 32)*5/9

    df1 <- data.frame(x1, y1)

    f1 <- function (times,temp_ini,mp, A,B) {
      T <- temp_ini+A*sin(B*(times-x1[1]))+ mp*(times-x1[1])
    }

    m1<- nls2(y1 ~ f1(x1,temp_ini,mp, A,B), df1, start = list( temp_ini=11.5,mp=1/6, A=1,B=0.2),control = list (maxiter = 500))
    y_est1<-predict(m1,df1$x1)

    times<-seq(x1[1],x1[length(x1)])
    temp1=f1(times,coef(m1)[1],coef(m1)[2],coef(m1)[3],coef(m1)[4])

    dat1<-data.frame('x'=times,'y'=temp1)
    dat2<-data.frame('x'=x1,'y'=y1)
    data<-rbind(dat1,dat2)
    p1<- ggplot(data, aes(x=.data$x, y=.data$y)) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_line(data =subset(dat1,times>x1[1]  & times<x1[length(x1)]), color = "black",size=rel(1.1))+
      geom_point(data =subset(dat2,x1>x1[1]  & x1<x1[length(x1)]), color = "black",size=rel(2.5))+
      labs(x = "Time (years)",y="Temperature")+
      coord_cartesian(xlim = c(x1[1],x1[length(x1)]))+
      theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
      theme(axis.title.x = element_text(size = rel(1.2), angle = 00))

    # ABUNDANCIA
    BIO3=data.frame(BIO$YEAR[31:45],BIO$ABUNDANCE[31:45])

    xa<- BIO3$BIO.YEAR.31.45.
    ya<- BIO3$BIO.ABUNDANCE.31.45.

    dfa <- data.frame(xa, ya)

    temp_ini<-coef(m1)[1]

    fa <- function (times,y_ini,temp_cmin,temp_ini, temp_cmax,ro,lambda) {

      parms1<-c(temp_cmin,temp_ini, temp_cmax,ro,lambda)

      model1 <- function (times,y,parms1) {
        with(as.list(c(y)), {

          T <- f1(times,coef(m1)[1],coef(m1)[2],coef(m1)[3],coef(m1)[4])
          temp_op1<- ( temp_cmax+temp_cmin)/3+sqrt((( temp_cmax+temp_cmin)/3)^2-
                                                     ( temp_cmax*temp_cmin)/3)
          r1<- ro*((T-temp_cmin)*(temp_cmax-T)*T)/((temp_op1-temp_cmin)*(temp_cmax-temp_op1)*temp_op1)
          dN <-   r1 * N * (1 - lambda*(N / r1))
          list(dN,T,r1) })
      }

      out1 <- ode(y=y_ini, times, model1 ,parms1,method = "ode45")

      return(out1[,2])
    }

    ma<- nls2(ya ~ fa(xa,y_ini,temp_cmin,temp_ini, temp_cmax,ro,lambda), dfa, alg = "plinear", start = list(ro=0.7, lambda=0.0005),control = nls.control (maxiter = 500))
    y_esta<-predict(ma,dfa$xa)

    times<-seq(x1[1],x1[length(x1)])
    Abund=fa(times,y_ini,temp_cmin,temp_ini, temp_cmax,coef(ma)[1],coef(ma)[2])
    da1<-data.frame('x'=times,'y'=Abund )
    da2<-data.frame('x'=xa,'y'=ya )

    data<-rbind(da1,da2)

    p2<-ggplot(data, aes(x=.data$x, y=.data$y)) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_line(data =subset(da1,times>xa[1]  & times<xa[length(xa)]), color = "black",size=rel(1.1))+
      geom_point(data =subset(da2,xa>xa[1]  & xa<xa[length(xa)]), color = "black",size=rel(2))+
      labs(x = "Time (years)",y="Abundance")+
      coord_cartesian(xlim = c(xa[1],xa[length(xa)]))+
      theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
      theme(axis.title.x = element_text(size = rel(1.2), angle = 00))


    plot_grid(p1, p2)

  } else if(trend== 5) {

    ##############################################################################

    #Temperatura promedio, constante
    x1<- TAVG$Temperature.DATE.31.45.
    y1<- (TAVG$Temperature.TAVG.31.45.- 32)*5/9

    df1 <- data.frame(x1, y1)
    time_start<-x1[1]
    f1 <- function (times,temp_ini) {
      T <- temp_ini
    }

    m1<- nls2(y1 ~ f1(x1,temp_ini), df1, start = list( temp_ini=11.5),control = list (maxiter = 500))
    y_est1<-predict(m1,df1$x1)

    times<-seq(x1[1],x1[length(x1)])
    temp1=rep(coef(m1)[1],length(times))

    dat1<-data.frame('x'=times,'y'=temp1)
    dat2<-data.frame('x'=x1,'y'=y1)
    data<-rbind(dat1,dat2)
    p1<- ggplot(data, aes(x=.data$x, y=.data$y)) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_line(data =subset(dat1,times>x1[1]  & times<x1[length(x1)]), color = "black",size=rel(1.1))+
      geom_point(data =subset(dat2,x1>x1[1]  & x1<x1[length(x1)]), color = "black",size=rel(2.5))+
      labs(x = "Time (years)",y="Temperature")+
      coord_cartesian(xlim = c(x1[1],x1[length(x1)]))+
      theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
      theme(axis.title.x = element_text(size = rel(1.2), angle = 00))

    # ABUNDANCIA
    BIO3=data.frame(BIO$YEAR[31:45],BIO$ABUNDANCE[31:45])

    xa<- BIO3$BIO.YEAR.31.45.
    ya<- BIO3$BIO.ABUNDANCE.31.45.

    dfa <- data.frame(xa, ya)

    temp_ini<-coef(m1)[1]

    fa <- function (times,y_ini,temp_cmin,temp_ini, temp_cmax,ro,lambda) {

      parms1<-c(temp_cmin,temp_ini, temp_cmax,ro,lambda)

      model1 <- function (times,y,parms1) {
        with(as.list(c(y)), {

          T <- rep(coef(m1)[1],length(times))
          temp_op1<- ( temp_cmax+temp_cmin)/3+sqrt((( temp_cmax+temp_cmin)/3)^2-
                                                     ( temp_cmax*temp_cmin)/3)
          r1<- ro*((T-temp_cmin)*(temp_cmax-T)*T)/((temp_op1-temp_cmin)*(temp_cmax-temp_op1)*temp_op1)
          dN <-   r1 * N * (1 - lambda*(N / r1))
          list(dN,T,r1) })
      }

      out1 <- ode(y=y_ini, times, model1 ,parms1,method = "ode45")

      return(out1[,2])
    }

    ma<- nls2(ya ~ fa(xa,y_ini,temp_cmin,temp_ini, temp_cmax,ro,lambda), dfa, alg = "plinear", start = list(ro=0.7, lambda=0.0005),control = nls.control (maxiter = 500))
    y_esta<-predict(ma,dfa$xa)

    times<-seq(x1[1],x1[length(x1)])
    Abund=fa(times,y_ini,temp_cmin,temp_ini, temp_cmax,coef(ma)[1],coef(ma)[2])
    da1<-data.frame('x'=times,'y'=Abund )
    da2<-data.frame('x'=xa,'y'=ya )

    data<-rbind(da1,da2)

    p2<-ggplot(data, aes(x=.data$x, y=.data$y)) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_line(data =subset(da1,times>xa[1]  & times<xa[length(xa)]), color = "black",size=rel(1.1))+
      geom_point(data =subset(da2,xa>xa[1]  & xa<xa[length(xa)]), color = "black",size=rel(2))+
      labs(x = "Time (years)",y="Abundance")+
      coord_cartesian(xlim = c(xa[1],xa[length(xa)]))+
      theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
      theme(axis.title.x = element_text(size = rel(1.2), angle = 00))


    plot_grid(p1, p2)

    #############################################################################



  }  else {

    print("The inserted code is not associated with an established temperature trend, please enter a valid code")

  }
}
