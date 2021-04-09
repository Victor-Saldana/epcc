#'Age structure under IPCC RCP2.6 or RCP8.5 scenarios
#'
#'@description This function allows simulating the effect of the IPCC RCP2.6 or RCP8.5 scenarios
#'             (2014) on the abundances of two species interacting through predation. The prey is an
#'             ectotherm population, and the predator is not affected by temperature.
#'
#'
#'
#'@param y_ini Initial population values (must be written with its name: N).
#'@param temp_ini Initial temperature.
#'@param temp_cmin Minimum critical temperature (prey).
#'@param temp_cmax Maximum critical temperature (prey).
#'@param ro Population growth rate at optimum temperature (prey).
#'@param lambda1 Marginal loss by non-thermodependent intraspecific competition (age1).
#'@param lambda2 Marginal loss by non-thermodependent intraspecific competition (age2).
#'@param lambda3 Marginal loss by non-thermodependent intraspecific competition (age3).
#'@param alpha1 Stage 1 to 2 transition coefficient.
#'@param alpha2 Stage 2 to 3 transition coefficient.
#'@param d2 Mortality rate (age2).
#'@param d3 Mortality rate (age3).
#'@param RCP Representative concentration trajectories (RCP2.6 and RCP8.5 scenarios).
#'@param time_start Start of time sequence.
#'@param time_end End of time sequence.
#'@param leap Time sequence step.
#'
#'@details Three scenarios can be evaluated for a predation interaction where the prey is an ectotherm
#'         population. The temperature trends correspond to IPCC projections under the RCP2.6 or RCP8.5
#'         scenarios. In each input vector, the parameters for the three simulations must be specified.
#'
#'
#'@return (1) A data.frame with columns having the simulated trends.
#'@return (2) A four-panel figure where (a), (b), and (c) show the population abundance curves for each simulation.
#'        The brown curve corresponds to the abundance of prey and the green curve to predators. Panel (d)
#'        shows the temperature trend curves used for each simulation, green, blue, and black, respectively.
#'
#'@references IPCC. (2014): Cambio climático 2014: Informe de síntesis. Contribución de los Grupos de trabajo
#'            I, II y III al Quinto Informe de Evaluación del Grupo Intergubernamental de Expertos sobre el
#'            Cambio Climático [Equipo principal de redacción, R.K. Pachauri y L.A. Meyer (eds.)]. IPCC, Ginebra,
#'            Suiza, 157 págs.
#'
#'@examples
#'
#'#######################################################################
#'   #Example 1: Different thermal tolerance ranges (scenario RCP2.6).
#'#######################################################################
#'
#'temp_cmin <- 18
#'
#'# Temperature that occurs before the minimum simulation time.
#'temp_i <- 22
#'
#'time_end <- 2100
#'
#'# Temperature that occurs in the maximum time ofthe simulation.
#'temp_max <- get_RCP2.6(time_end)+temp_i
#'
#'# Simulation thermal range.
#'RS <- temp_max-temp_cmin
#'
#'temp_cmax1 <- 4/3*RS+temp_cmin
#'temp_cmax2 <- 2/3*RS+temp_cmin
#'temp_cmax3 <- 1/3*RS+temp_cmin
#'temp_ini <- (temp_cmin+temp_cmax3)/2
#'
#'age_structure(y_ini = c(N1 = 800, N1 = 800, N1 = 800,
#'                        N2 = 600, N2 = 600, N2 = 600,
#'                        N3 = 400, N3 = 400, N3 = 400),
#'              temp_ini = rep(temp_ini,3),
#'              temp_cmin = rep(temp_cmin,3),
#'              temp_cmax = c(temp_cmax1,temp_cmax2,temp_cmax3),
#'              ro = rep(0.7,3),
#'              lambda1 = rep(0.0004,3),
#'              lambda2 = rep(0.0006,3),
#'              lambda3 = rep(0.0005,3),
#'              alpha1 = 0.3,
#'              alpha2 = 0.4,
#'              d2 = 0.2,
#'              d3 = 0.3,
#'              RCP = 2.6,
#'              time_start = 2005,
#'              time_end = time_end,
#'              leap = 1/50)
#'
#'#######################################################################
#'   #Example 2: Different thermal tolerance ranges (scenario RCP8.5).
#'#######################################################################
#'
#'temp_cmin <- 18
#'
#'# Temperature that occurs before the minimum simulation time.
#'temp_i <- 22
#'
#'time_end <- 2100
#'
#'# Temperature that occurs in the maximum time ofthe simulation.
#'temp_max <- get_RCP8.5(time_end,a=exp(coef(m)[1]), b=coef(m)[2])+temp_i
#'
#'# Simulation thermal range.
#'RS <- temp_max-temp_cmin
#'
#'temp_cmax1 <- 4/3*RS+temp_cmin
#'temp_cmax2 <- 2/3*RS+temp_cmin
#'temp_cmax3 <- 1/3*RS+temp_cmin
#'temp_ini <- (temp_cmin+temp_cmax3)/2
#'
#'age_structure(y_ini = c(N1 = 800, N1 = 800, N1 = 800,
#'                        N2 = 600, N2 = 600, N2 = 600,
#'                        N3 = 400, N3 = 400, N3 = 400),
#'              temp_ini = rep(temp_ini,3),
#'              temp_cmin = rep(temp_cmin,3),
#'              temp_cmax = c(temp_cmax1,temp_cmax2,temp_cmax3),
#'              ro = rep(0.7,3),
#'              lambda1 = rep(0.0004,3),
#'              lambda2 = rep(0.0006,3),
#'              lambda3 = rep(0.0005,3),
#'              alpha1 = 0.3,
#'              alpha2 = 0.4,
#'              d2 = 0.3,
#'              d3 = 0.4,
#'              RCP = 8.5,
#'              time_start = 2005,
#'              time_end = time_end,
#'              leap = 1/50)
#'
#'#######################################################################
#'   #Example 3: Different conversion efficiencies (scenario RCP2.6).
#'#######################################################################
#'
#'e1 <- 0.2
#'e2 <- 2*e1
#'e3 <- 2*e2
#'
#'age_structure(y_ini = c(N1 = 800, N1 = 800, N1 = 800,
#'                        N2 = 600, N2 = 600, N2 = 600,
#'                        N3 = 400, N3 = 400, N3 = 400),
#'              temp_ini = rep(22,3),
#'              temp_cmin = rep(20,3),
#'              temp_cmax = rep(35,3),
#'              ro = rep(0.9,3),
#'              lambda1 = rep(0.0004,3),
#'              lambda2 = rep(0.0006,3),
#'              lambda3 = rep(0.0005,3),
#'              alpha1 = 0.3,
#'              alpha2 = 0.4,
#'              d2 = 0.3,
#'              d3 = 0.4,
#'              RCP = 2.6,
#'              time_start = 2005,
#'              time_end = 2100,
#'              leap = 1/50)
#'
#'#######################################################################
#'   #Example 4: Different conversion efficiencies (scenario RCP8.5).
#'#######################################################################
#'
#'e1 <- 0.2
#'e2 <- 2*e1
#'e3 <- 2*e2
#'
#'age_structure(y_ini = c(N1 = 800, N1 = 800, N1 = 800,
#'                        N2 = 600, N2 = 600, N2 = 600,
#'                        N3 = 400, N3 = 400, N3 = 400),
#'              temp_ini = rep(22,3),
#'              temp_cmin = rep(20,3),
#'              temp_cmax = rep(35,3),
#'              ro = rep(0.9,3),
#'              lambda1 = rep(0.0004,3),
#'              lambda2 = rep(0.0006,3),
#'              lambda3 = rep(0.0005,3),
#'              alpha1 = 0.3,
#'              alpha2 = 0.4,
#'              d2 = 0.3,
#'              d3 = 0.4,
#'              RCP = 8.5,
#'              time_start = 2005,
#'              time_end = 2100,
#'              leap = 1/50)



age_structure<- function(y_ini = c(N1 = 800, N1 = 800, N1 = 800,
                        N2 = 600, N2 = 600, N2 = 600,
                        N3 = 400, N3 = 400, N3 = 400),
              temp_ini = rep(25,3),
              temp_cmin = rep(18,3),
              temp_cmax = c(25,28,35),
              ro = rep(0.7,3),
              lambda1 = rep(0.0004,3),
              lambda2 = rep(0.0004,3),
              lambda3 = rep(0.0004,3),
              alpha1 = 0.1,
              alpha2 = 0.7,
              d2 = 0.4,
              d3 = 0.4,
              RCP = 2.6,
              time_start = 2005,
              time_end = 2100,
              leap = 1/50){

  times<- seq(time_start, time_end, leap)

  if(time_end<=2100){
    if(time_start<=time_end){

      if(temp_cmin[1]<temp_cmax[1] && temp_cmin[2]<temp_cmax[2] && temp_cmin[3]<temp_cmax[3] ){


        if(temp_cmin[1]<=temp_ini[1] && temp_ini[1]<=temp_cmax[1] && temp_cmin[2]<=temp_ini[2] && temp_ini[2]<=temp_cmax[2] && temp_cmin[3]<=temp_ini[3] && temp_ini[3]<=temp_cmax[3]){

          options(warn = - 1)

          ##########################################################
          # Optimum growing temperature
          ##########################################################

          temp_op1<- (temp_cmax[1]+temp_cmin[1])/3+sqrt(((temp_cmax[1]+temp_cmin[1])/3)^2-(temp_cmax[1]*temp_cmin[1])/3)

          temp_op2<- (temp_cmax[2]+temp_cmin[2])/3+sqrt(((temp_cmax[2]+temp_cmin[2])/3)^2-(temp_cmax[2]*temp_cmin[2])/3)

          temp_op3<- (temp_cmax[3]+temp_cmin[3])/3+sqrt(((temp_cmax[3]+temp_cmin[3])/3)^2-(temp_cmax[3]*temp_cmin[3])/3)

          ##########################################################
          # Parameters
          ##########################################################

          parms1<-c(temp_cmin[1],temp_ini[1],temp_cmax[1],temp_op1,ro[1], lambda1[1],lambda1[2],lambda1[3],lambda2[1],lambda2[2],lambda2[3],lambda3[1],lambda3[2],lambda3[3],alpha1,alpha2,d2,d3)
          parms2<-c(temp_cmin[2],temp_ini[2],temp_cmax[2],temp_op2,ro[2], lambda1[1],lambda1[2],lambda1[3],lambda2[1],lambda2[2],lambda2[3],lambda3[1],lambda3[2],lambda3[3],alpha1,alpha2,d2,d3)
          parms3<-c(temp_cmin[3],temp_ini[3],temp_cmax[3],temp_op3,ro[3], lambda1[1],lambda1[2],lambda1[3],lambda2[1],lambda2[2],lambda2[3],lambda3[1],lambda3[2],lambda3[3],alpha1,alpha2,d2,d3)

          ##############################################



          if(RCP==2.6) {
            ##########################################################
            # Model for each trend
            ##########################################################
            model1 <- function (times, y,parms1) {
              with(as.list(c(y)), {
                T1 <- get_RCP2.6(times)+temp_ini[1]  # IPCC1
                r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
                dN1 <-  r1*N1-lambda1[1]*N1*N1-lambda1[2]*N1*N2-lambda1[3]*N1*N3-alpha1*N1
                dN2 <-  alpha1*N1-d2*N2-lambda2[1]*N2*N1-lambda2[2]*N2*N2-lambda2[3]*N2*N3-alpha2*N2
                dN3 <-  alpha2*N2-d3*N3-lambda3[1]*N3*N1-lambda3[2]*N3*N2-lambda3[3]*N3*N3

                return(list(c(dN1,dN2,dN3)))
              })
            }
            ###############################################################

            model2 <- function (times, y,parms2) {
              with(as.list(c(y)), {
                T2 <- get_RCP2.6(times)+temp_ini[2]  # IPCC1
                r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
                dN1 <-  r2*N1-lambda1[1]*N1*N1-lambda1[2]*N1*N2-lambda1[3]*N1*N3-alpha1*N1
                dN2 <-  alpha1*N1-d2*N2-lambda2[1]*N2*N1-lambda2[2]*N2*N2-lambda2[3]*N2*N3-alpha2*N2
                dN3 <-  alpha2*N2-d3*N3-lambda3[1]*N3*N1-lambda3[2]*N3*N2-lambda3[3]*N3*N3

                return(list(c(dN1,dN2,dN3)))
              })
            }
            ###############################################################

            model3 <- function (times, y,parms3) {
              with(as.list(c(y)), {
                T3 <- get_RCP2.6(times)+temp_ini[3]  # IPCC1
                r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)
                dN1 <-  r3*N1-lambda1[1]*N1*N1-lambda1[2]*N1*N2-lambda1[3]*N1*N3-alpha1*N1
                dN2 <-  alpha1*N1-d2*N2-lambda2[1]*N2*N1-lambda2[2]*N2*N2-lambda2[3]*N2*N3-alpha2*N2
                dN3 <-  alpha2*N2-d3*N3-lambda3[1]*N3*N1-lambda3[2]*N3*N2-lambda3[3]*N3*N3

                return(list(c(dN1,dN2,dN3)))
              })
            }
            ###############################################################
            y_ini1<-c(y_ini[1],y_ini[4],y_ini[7])
            y_ini2<-c(y_ini[2],y_ini[5],y_ini[8])
            y_ini3<-c(y_ini[3],y_ini[6],y_ini[9])

            ###############################################################
            # Solution
            ##############################################################

            out1 <- ode(y=y_ini1, times, model1, parms1, method = "ode45")
            out2 <- ode(y=y_ini2, times, model2, parms2, method = "ode45")
            out3 <- ode(y=y_ini3, times, model3, parms3, method = "ode45")
            #############################################################


            ###############################################################
            # Abundance
            ##############################################################

            data1<-data.frame('x'=times,'y'=out1[,2])
            data2<-data.frame('x'=times,'y'=out1[,3])
            data3<-data.frame('x'=times,'y'=out1[,4])

            # data1$group<-"Pop1"
            # data2$group<-"Pop2"
            # data3$group<-"Pop3"


            dat1<-data.frame('x'=times,'y'=out2[,2])
            dat2<-data.frame('x'=times,'y'=out2[,3])
            dat3<-data.frame('x'=times,'y'=out2[,4])


            # dat1$group<-"Pop1"
            # dat2$group<-"Pop2"
            # dat3$group<-"Pop3"


            da1<-data.frame('x'=times,'y'=out3[,2])
            da2<-data.frame('x'=times,'y'=out3[,3])
            da3<-data.frame('x'=times,'y'=out3[,4])

            # da1$group<-"Pop1"
            # da2$group<-"Pop2"
            # da3$group<-"Pop3"

            T1 <- get_RCP2.6(times)+temp_ini[1]
            T2 <- get_RCP2.6(times)+temp_ini[2]
            T3 <- get_RCP2.6(times)+temp_ini[3]

            d1<-data.frame('x'=times,'y'=T1)
            d2<-data.frame('x'=times,'y'=T2)
            d3<-data.frame('x'=times,'y'=T3)

            # d1$group<-"Pop1"
            # d2$group<-"Pop2"
            # d3$group<-"Pop3"

            r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
            r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
            r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)

            # K1=r1/lambda[1]
            # K2=r2/lambda[2]
            # K3=r3/lambda[3]
            #
            #
            # cap1<-data.frame('x'=times,'y'=K1 )
            # cap2<-data.frame('x'=times,'y'=K2 )
            # cap3<-data.frame('x'=times,'y'=K3 )
            #
            # cap1$group<-"Pop1"
            # cap2$group<-"Pop2"
            # cap3$group<-"Pop3"


            ###############################################################
            # Data
            ###############################################################

            Data<- data.frame(times,out1[,2],out1[,3], out1[,4],out2[,2],out2[,3],out2[,4],out3[,2],out3[,3],out3[,4])
            names(Data)<- c("Time","Abundance species-1","Abundance species-2","Abundance species-3","Abundance species-1","Abundance species-2","Abundance species-3","Abundance species-1","Abundance species-2","Abundance species-3")
            View(Data)


            times_new1<-vector(mode = "numeric", length = 0)
            times_new2<-vector(mode = "numeric", length = 0)
            times_new3<-vector(mode = "numeric", length = 0)

            for (i in 2: length(times)){

              if(out1[i-1,2]>=0 && ( out1[i,2])<0){
                times_new1[i-1]<- times[i-1]
              }else{
                times_new1[i-1]<- 0
              }
            }


            for (i in 2: length(times)){

              if(out2[i-1,2]>=0 && ( out2[i,2])<0){
                times_new2[i-1]<- times[i-1]
              }else{
                times_new2[i-1]<- 0
              }
            }


            for (i in 2: length(times)){

              if(out3[i-1,2]>=0 && ( out3[i,2])<0){
                times_new3[i-1]<- times[i-1]
              }else{
                times_new3[i-1]<- 0
              }
            }

            index1<- which(times_new1!=0)[1]
            index2<- which(times_new2!=0)[1]
            index3<- which(times_new3!=0)[1]

            index1<- as.integer(index1)
            index2<- as.integer(index2)
            index3<- as.integer(index3)

            if(!is.na(as.integer(index1))== FALSE){
              times_sup11<- times[length(times)]
            }else{
              times_sup11<- times[index1]
            }
            if(!is.na(as.integer(index2))== FALSE){
              times_sup21<- times[length(times)]
            }else{
              times_sup21<- times[index2]
            }

            if(!is.na(as.integer(index3))== FALSE){
              times_sup31<- times[length(times)]
            }else{
              times_sup31<- times[index3]
            }

            times_new7<-vector(mode = "numeric", length = 0)
            times_new8<-vector(mode = "numeric", length = 0)
            times_new9<-vector(mode = "numeric", length = 0)

            for (i in 2: length(times)){

              if(( temp_cmax[1]-T1[i-1])>=0 && ( temp_cmax[1]-T1[i])<0){
                times_new7[i-1]<- times[i-1]

              }else if(( temp_cmax[1]-T1[i-1])<=0 && ( temp_cmax[1]-T1[i])>0){

                times_new7[i-1]<- times[i-1]
              }else{
                times_new7[i-1]<- 0
              }
            }


            for (i in 2: length(times)){

              if(( temp_cmax[2]-T2[i-1])>=0 && ( temp_cmax[2]-T2[i])<0){
                times_new8[i-1]<- times[i-1]

              }else if(( temp_cmax[2]-T2[i-1])<=0 && ( temp_cmax[2]-T2[i])>0){

                times_new8[i-1]<- times[i-1]
              }else{
                times_new8[i-1]<- 0
              }
            }


            for (i in 2: length(times)){

              if(( temp_cmax[3]-T3[i-1])>=0 && ( temp_cmax[3]-T3[i])<0){
                times_new9[i-1]<- times[i-1]

              }else if(( temp_cmax[3]-T3[i-1])<=0 && ( temp_cmax[3]-T3[i])>0){

                times_new9[i-1]<- times[i-1]
              }else{
                times_new9[i-1]<- 0
              }
            }

            index7<- which(times_new7!=0)[1]
            index8<- which(times_new8!=0)[1]
            index9<- which(times_new9!=0)[1]

            index7<- as.integer(index7)
            index8<- as.integer(index8)
            index9<- as.integer(index9)

            if(!is.na(as.integer(index7))== FALSE){
              times_sup12<- times[length(times)]
            }else{
              times_sup12<- times[index7]
            }
            if(!is.na(as.integer(index8))== FALSE){
              times_sup22<- times[length(times)]
            }else{
              times_sup22<- times[index8]
            }

            if(!is.na(as.integer(index9))== FALSE){
              times_sup32<- times[length(times)]
            }else{
              times_sup32<- times[index9]
            }


            if(times_sup11<= times_sup12){
              times_sup1<-times_sup11
            }else{
              times_sup1<-times_sup12
            }

            if(times_sup21<= times_sup22){
              times_sup2<-times_sup21
            }else{
              times_sup2<-times_sup22
            }

            if(times_sup31<= times_sup32){
              times_sup3<-times_sup31
            }else{
              times_sup3<-times_sup32
            }
            ###############################################################
            # Carrying capacity
            ##############################################################

            # times_new4<-vector(mode = "numeric", length = 0)
            # times_new5<-vector(mode = "numeric", length = 0)
            # times_new6<-vector(mode = "numeric", length = 0)
            #
            # for (i in 2: length(times)){
            #
            #   if(K1[i-1]>=0 && ( K1[i])<0){
            #     times_new4[i-1]<- times[i-1]
            #   }else{
            #     times_new4[i-1]<- 0
            #   }
            # }
            #
            #
            # for (i in 2: length(times)){
            #
            #   if(K2[i-1]>=0 && ( K2[i])<0){
            #     times_new5[i-1]<- times[i-1]
            #   }else{
            #     times_new5[i-1]<- 0
            #   }
            # }
            #
            #
            # for (i in 2: length(times)){
            #
            #   if(K3[i-1]>=0 && ( K3[i])<0){
            #     times_new6[i-1]<- times[i-1]
            #   }else{
            #     times_new6[i-1]<- 0
            #   }
            # }
            #
            # index4<- which(times_new4!=0)[1]
            # index5<- which(times_new5!=0)[1]
            # index6<- which(times_new6!=0)[1]
            #
            # index4<- as.integer(index4)
            # index5<- as.integer(index5)
            # index6<- as.integer(index6)
            #
            # if(!is.na(as.integer(index4))== FALSE){
            #   times_sup4<- times[length(times)]
            # }else{
            #   times_sup4<- times[index4]
            # }
            # if(!is.na(as.integer(index5))== FALSE){
            #   times_sup5<- times[length(times)]
            # }else{
            #   times_sup5<- times[index5]
            # }
            #
            # if(!is.na(as.integer(index6))== FALSE){
            #   times_sup6<- times[length(times)]
            # }else{
            #   times_sup6<- times[index6]
            # }


            ###############################################################
            # Plots
            ##############################################################

            data<-rbind(data1, data2, data3)

            p1<- ggplot(data, aes(x=x, y=y)) +
              theme_bw()+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
              # geom_ribbon(data=subset(cap1,x>times[1] & x<times_sup4),aes(x=x,ymax=y),ymin=0,alpha=0.3, fill="brown") +
              geom_vline(xintercept = times_sup1, size=.5, color="brown",linetype="dashed")+
              geom_line(data =subset(data1,x>times[1] & x<times_sup1), color = "brown")+
              geom_line(data =subset(data2,x>times[1] & x<times_sup1), color = "green4")+
              geom_line(data =subset(data3,x>times[1] & x<times_sup1), color = "blue")+
              labs(x = "Time",y="Abundance")+
              theme(plot.title = element_text(size=40))+
              theme(plot.title = element_text(hjust = 0.5))+
              theme(axis.title.y = element_text(size = rel(1), angle = 90))+
              theme(axis.title.x = element_text(size = rel(1), angle = 00))+
              labs(tag = "(a)")



            dat<-rbind(dat1, dat2, dat3)

            p2<- ggplot(dat, aes(x=x, y=y)) +
              theme_bw()+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
              # geom_ribbon(data=subset(cap2,x>times[1] & x<times_sup5),aes(x=x,ymax=y),ymin=0,alpha=0.3, fill="brown") +
              geom_vline(xintercept = times_sup2, size=.5, color="brown",linetype="dashed")+
              geom_line(data =subset(dat1,x>times[1] & x<times_sup2), color = "brown")+
              geom_line(data =subset(dat2,x>times[1] & x<times_sup2), color = "green4")+
              geom_line(data =subset(dat3,x>times[1] & x<times_sup2), color = "blue")+
              labs(x = "Time",y="Abundance")+
              theme(plot.title = element_text(size=40))+
              theme(plot.title = element_text(hjust = 0.5))+
              theme(axis.title.y = element_text(size = rel(1), angle = 90))+
              theme(axis.title.x = element_text(size = rel(1), angle = 00))+
              labs(tag = "(b)")


            da<-rbind(da1, da2, da3)

            p3<- ggplot(da, aes(x=x, y=y)) +
              theme_bw()+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
              # geom_ribbon(data=subset(cap3,x>times[1] & x<times_sup6),aes(x=x,ymax=y),ymin=0,alpha=0.3, fill="brown") +
              geom_vline(xintercept = times_sup3, size=.5, color="brown",linetype="dashed")+
              geom_line(data =subset(da1,x>times[1] & x<times_sup3), color = "brown")+
              geom_line(data =subset(da2,x>times[1] & x<times_sup3), color = "green4")+
              geom_line(data =subset(da3,x>times[1] & x<times_sup3), color = "blue")+
              labs(x = "Time",y="Abundance")+
              theme(plot.title = element_text(size=40))+
              theme(plot.title = element_text(hjust = 0.5))+
              theme(axis.title.y = element_text(size = rel(1), angle = 90))+
              theme(axis.title.x = element_text(size = rel(1), angle = 00))+
              labs(tag = "(c)")


            d<-rbind(d1, d2, d3)

            p4<- ggplot(d, aes(x=x, y=y)) +
              theme_bw()+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
              geom_vline(xintercept = times_sup1, size=.5, color="green",linetype="dashed")+
              geom_vline(xintercept = times_sup2, size=.5, color="blue",linetype="dashed")+
              geom_vline(xintercept = times_sup3, size=.5, color="black",linetype="dashed")+
              geom_line(data =subset(d1,x>times[1] & x<times_sup1), color = "green")+
              geom_line(data =subset(d2,x>times[1] & x<times_sup2), color = "blue")+
              geom_line(data =subset(d3,x>times[1] & x<times_sup3), color = "black")+
              labs(x = "Time",y="Temperature")+
              theme(plot.title = element_text(size=40))+
              theme(plot.title = element_text(hjust = 0.5))+
              theme(axis.title.y = element_text(size = rel(1), angle = 90))+
              theme(axis.title.x = element_text(size = rel(1), angle = 00))+
              labs(tag = "(d)")

            plot_grid(p1, p2,p3,p4)

          } else if(RCP==8.5) {

            ##########################################################
            # Model for each trend
            ##########################################################
            model1 <- function (times, y,parms1) {
              with(as.list(c(y)), {
                T1<-  get_RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[1]    #IPCC2
                r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
                dN1 <-  r1*N1-lambda1[1]*N1*N1-lambda1[2]*N1*N2-lambda1[3]*N1*N3-alpha1*N1
                dN2 <-  alpha1*N1-d2*N2-lambda2[1]*N2*N1-lambda2[2]*N2*N2-lambda2[3]*N2*N3-alpha2*N2
                dN3 <-  alpha2*N2-d3*N3-lambda3[1]*N3*N1-lambda3[2]*N3*N2-lambda3[3]*N3*N3

                return(list(c(dN1,dN2,dN3)))
              })
            }
            ###############################################################

            model2 <- function (times, y,parms2) {
              with(as.list(c(y)), {
                T2<-  get_RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[2]    #IPCC2
                r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
                dN1 <-  r2*N1-lambda1[1]*N1*N1-lambda1[2]*N1*N2-lambda1[3]*N1*N3-alpha1*N1
                dN2 <-  alpha1*N1-d2*N2-lambda2[1]*N2*N1-lambda2[2]*N2*N2-lambda2[3]*N2*N3-alpha2*N2
                dN3 <-  alpha2*N2-d3*N3-lambda3[1]*N3*N1-lambda3[2]*N3*N2-lambda3[3]*N3*N3

                return(list(c(dN1,dN2,dN3)))
              })
            }
            ###############################################################

            model3 <- function (times, y,parms3) {
              with(as.list(c(y)), {
                T3<-  get_RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[3]    #IPCC2
                r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)
                dN1 <-  r3*N1-lambda1[1]*N1*N1-lambda1[2]*N1*N2-lambda1[3]*N1*N3-alpha1*N1
                dN2 <-  alpha1*N1-d2*N2-lambda2[1]*N2*N1-lambda2[2]*N2*N2-lambda2[3]*N2*N3-alpha2*N2
                dN3 <-  alpha2*N2-d3*N3-lambda3[1]*N3*N1-lambda3[2]*N3*N2-lambda3[3]*N3*N3

                return(list(c(dN1,dN2,dN3)))
              })
            }
            ###############################################################
            y_ini1<-c(y_ini[1],y_ini[4],y_ini[7])
            y_ini2<-c(y_ini[2],y_ini[5],y_ini[8])
            y_ini3<-c(y_ini[3],y_ini[6],y_ini[9])

            ###############################################################
            # Solution
            ##############################################################

            out1 <- ode(y=y_ini1, times, model1, parms1,method = "ode45")
            out2 <- ode(y=y_ini2, times, model2, parms2,method = "ode45")
            out3 <- ode(y=y_ini3, times, model3, parms3,method = "ode45")
            #############################################################


            ###############################################################
            # Abundance
            ##############################################################

            data1<-data.frame('x'=times,'y'=out1[,2] )
            data2<-data.frame('x'=times,'y'=out1[,3] )
            data3<-data.frame('x'=times,'y'=out1[,4] )

            # data1$group<-"Pop1"
            # data2$group<-"Pop2"
            # data3$group<-"Pop3"


            dat1<-data.frame('x'=times,'y'=out2[,2] )
            dat2<-data.frame('x'=times,'y'=out2[,3] )
            dat3<-data.frame('x'=times,'y'=out2[,4] )


            # dat1$group<-"Pop1"
            # dat2$group<-"Pop2"
            # dat3$group<-"Pop3"

            da1<-data.frame('x'=times,'y'=out3[,2] )
            da2<-data.frame('x'=times,'y'=out3[,3] )
            da3<-data.frame('x'=times,'y'=out3[,4] )

            # da1$group<-"Pop1"
            # da2$group<-"Pop2"
            # da3$group<-"Pop3"

            T1<-  get_RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[1]
            T2<-  get_RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[2]
            T3<-  get_RCP8.5(times,a=exp(coef(m)[1]), b=coef(m)[2])+temp_ini[3]

            d1<-data.frame('x'=times,'y'=T1)
            d2<-data.frame('x'=times,'y'=T2)
            d3<-data.frame('x'=times,'y'=T3)

            # d1$group<-"Pop1"
            # d2$group<-"Pop2"
            # d3$group<-"Pop3"

            r1<- rate_TPC(T1,ro[1],temp_cmin[1],temp_cmax[1],temp_op1)
            r2<- rate_TPC(T2,ro[2],temp_cmin[2],temp_cmax[2],temp_op2)
            r3<- rate_TPC(T3,ro[3],temp_cmin[3],temp_cmax[3],temp_op3)

            # K1=r1/lambda[1]
            # K2=r2/lambda[2]
            # K3=r3/lambda[3]
            #
            #
            # cap1<-data.frame('x'=times,'y'=K1 )
            # cap2<-data.frame('x'=times,'y'=K2 )
            # cap3<-data.frame('x'=times,'y'=K3 )
            #
            # cap1$group<-"Pop1"
            # cap2$group<-"Pop2"
            # cap3$group<-"Pop3"

            ###############################################################
            # Data
            ###############################################################

            Data<- data.frame(times,out1[,2],out1[,3],out1[,4],out2[,2],out2[,3],out2[,4],out3[,2],out3[,3],out3[,4])
            names(Data)<- c("Time","Abundance species-1","Abundance species-2","Abundance species-3","Abundance species-1","Abundance species-2","Abundance species-3","Abundance species-1","Abundance species-2","Abundance species-3")
            View(Data)

            times_new1<-vector(mode = "numeric", length = 0)
            times_new2<-vector(mode = "numeric", length = 0)
            times_new3<-vector(mode = "numeric", length = 0)

            for (i in 2: length(times)){

              if(out1[i-1,2]>=0 && ( out1[i,2])<0){
                times_new1[i-1]<- times[i-1]
              }else{
                times_new1[i-1]<- 0
              }
            }


            for (i in 2: length(times)){

              if(out2[i-1,2]>=0 && ( out2[i,2])<0){
                times_new2[i-1]<- times[i-1]
              }else{
                times_new2[i-1]<- 0
              }
            }


            for (i in 2: length(times)){

              if(out3[i-1,2]>=0 && ( out3[i,2])<0){
                times_new3[i-1]<- times[i-1]
              }else{
                times_new3[i-1]<- 0
              }
            }

            index1<- which(times_new1!=0)[1]
            index2<- which(times_new2!=0)[1]
            index3<- which(times_new3!=0)[1]

            index1<- as.integer(index1)
            index2<- as.integer(index2)
            index3<- as.integer(index3)

            if(!is.na(as.integer(index1))== FALSE){
              times_sup11<- times[length(times)]
            }else{
              times_sup11<- times[index1]
            }
            if(!is.na(as.integer(index2))== FALSE){
              times_sup21<- times[length(times)]
            }else{
              times_sup21<- times[index2]
            }

            if(!is.na(as.integer(index3))== FALSE){
              times_sup31<- times[length(times)]
            }else{
              times_sup31<- times[index3]
            }


            times_new7<-vector(mode = "numeric", length = 0)
            times_new8<-vector(mode = "numeric", length = 0)
            times_new9<-vector(mode = "numeric", length = 0)

            for (i in 2: length(times)){

              if(( temp_cmax[1]-T1[i-1])>=0 && ( temp_cmax[1]-T1[i])<0){
                times_new7[i-1]<- times[i-1]

              }else if(( temp_cmax[1]-T1[i-1])<=0 && ( temp_cmax[1]-T1[i])>0){

                times_new7[i-1]<- times[i-1]
              }else{
                times_new7[i-1]<- 0
              }
            }


            for (i in 2: length(times)){

              if(( temp_cmax[2]-T2[i-1])>=0 && ( temp_cmax[2]-T2[i])<0){
                times_new8[i-1]<- times[i-1]

              }else if(( temp_cmax[2]-T2[i-1])<=0 && ( temp_cmax[2]-T2[i])>0){

                times_new8[i-1]<- times[i-1]
              }else{
                times_new8[i-1]<- 0
              }
            }


            for (i in 2: length(times)){

              if(( temp_cmax[3]-T3[i-1])>=0 && ( temp_cmax[3]-T3[i])<0){
                times_new9[i-1]<- times[i-1]

              }else if(( temp_cmax[3]-T3[i-1])<=0 && ( temp_cmax[3]-T3[i])>0){

                times_new9[i-1]<- times[i-1]
              }else{
                times_new9[i-1]<- 0
              }
            }

            index7<- which(times_new7!=0)[1]
            index8<- which(times_new8!=0)[1]
            index9<- which(times_new9!=0)[1]

            index7<- as.integer(index7)
            index8<- as.integer(index8)
            index9<- as.integer(index9)

            if(!is.na(as.integer(index7))== FALSE){
              times_sup12<- times[length(times)]
            }else{
              times_sup12<- times[index7]
            }
            if(!is.na(as.integer(index8))== FALSE){
              times_sup22<- times[length(times)]
            }else{
              times_sup22<- times[index8]
            }

            if(!is.na(as.integer(index9))== FALSE){
              times_sup32<- times[length(times)]
            }else{
              times_sup32<- times[index9]
            }


            if(times_sup11<= times_sup12){
              times_sup1<-times_sup11
            }else{
              times_sup1<-times_sup12
            }

            if(times_sup21<= times_sup22){
              times_sup2<-times_sup21
            }else{
              times_sup2<-times_sup22
            }

            if(times_sup31<= times_sup32){
              times_sup3<-times_sup31
            }else{
              times_sup3<-times_sup32
            }


            ###############################################################
            # Carrying capacity
            ##############################################################

            # times_new4<-vector(mode = "numeric", length = 0)
            # times_new5<-vector(mode = "numeric", length = 0)
            # times_new6<-vector(mode = "numeric", length = 0)
            #
            # for (i in 2: length(times)){
            #
            #   if(K1[i-1]>=0 && ( K1[i])<0){
            #     times_new4[i-1]<- times[i-1]
            #   }else{
            #     times_new4[i-1]<- 0
            #   }
            # }
            #
            #
            # for (i in 2: length(times)){
            #
            #   if(K2[i-1]>=0 && ( K2[i])<0){
            #     times_new5[i-1]<- times[i-1]
            #   }else{
            #     times_new5[i-1]<- 0
            #   }
            # }
            #
            #
            # for (i in 2: length(times)){
            #
            #   if(K3[i-1]>=0 && ( K3[i])<0){
            #     times_new6[i-1]<- times[i-1]
            #   }else{
            #     times_new6[i-1]<- 0
            #   }
            # }
            #
            # index4<- which(times_new4!=0)[1]
            # index5<- which(times_new5!=0)[1]
            # index6<- which(times_new6!=0)[1]
            #
            # index4<- as.integer(index4)
            # index5<- as.integer(index5)
            # index6<- as.integer(index6)
            #
            # if(!is.na(as.integer(index4))== FALSE){
            #   times_sup4<- times[length(times)]
            # }else{
            #   times_sup4<- times[index4]
            # }
            # if(!is.na(as.integer(index5))== FALSE){
            #   times_sup5<- times[length(times)]
            # }else{
            #   times_sup5<- times[index5]
            # }
            #
            # if(!is.na(as.integer(index6))== FALSE){
            #   times_sup6<- times[length(times)]
            # }else{
            #   times_sup6<- times[index6]
            # }
            #

            ###############################################################
            # Plots
            ##############################################################

            data<-rbind(data1, data2, data3)

            p1<- ggplot(data, aes(x=x, y=y)) +
              theme_bw()+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
              # geom_ribbon(data=subset(cap1,x>times[1] & x<times_sup4),aes(x=x,ymax=y),ymin=0,alpha=0.3, fill="brown") +
              geom_vline(xintercept = times_sup1, size=.5, color="brown",linetype="dashed")+
              geom_line(data =subset(data1,x>times[1] & x<times_sup1), color = "brown")+
              geom_line(data =subset(data2,x>times[1] & x<times_sup1), color = "green4")+
              geom_line(data =subset(data3,x>times[1] & x<times_sup1), color = "blue")+
              labs(x = "Time",y="Abundance")+
              theme(plot.title = element_text(size=40))+
              theme(plot.title = element_text(hjust = 0.5))+
              theme(axis.title.y = element_text(size = rel(1), angle = 90))+
              theme(axis.title.x = element_text(size = rel(1), angle = 00))+
              labs(tag = "(a)")



            dat<-rbind(dat1, dat2, dat3)

            p2<- ggplot(dat, aes(x=x, y=y)) +
              theme_bw()+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
              # geom_ribbon(data=subset(cap2,x>times[1] & x<times_sup5),aes(x=x,ymax=y),ymin=0,alpha=0.3, fill="brown") +
              geom_vline(xintercept = times_sup2, size=.5, color="brown",linetype="dashed")+
              geom_line(data =subset(dat1,x>times[1] & x<times_sup2), color = "brown")+
              geom_line(data =subset(dat2,x>times[1] & x<times_sup2), color = "green4")+
              geom_line(data =subset(dat3,x>times[1] & x<times_sup2), color = "blue")+
              labs(x = "Time",y="Abundance")+
              theme(plot.title = element_text(size=40))+
              theme(plot.title = element_text(hjust = 0.5))+
              theme(axis.title.y = element_text(size = rel(1), angle = 90))+
              theme(axis.title.x = element_text(size = rel(1), angle = 00))+
              labs(tag = "(b)")


            da<-rbind(da1, da2, da3)

            p3<- ggplot(da, aes(x=x, y=y)) +
              theme_bw()+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
              # geom_ribbon(data=subset(cap3,x>times[1] & x<times_sup6),aes(x=x,ymax=y),ymin=0,alpha=0.3, fill="brown") +
              geom_vline(xintercept = times_sup3, size=.5, color="brown",linetype="dashed")+
              geom_line(data =subset(da1,x>times[1] & x<times_sup3), color = "brown")+
              geom_line(data =subset(da2,x>times[1] & x<times_sup3), color = "green4")+
              geom_line(data =subset(da3,x>times[1] & x<times_sup3), color = "blue")+
              labs(x = "Time",y="Abundance")+
              theme(plot.title = element_text(size=40))+
              theme(plot.title = element_text(hjust = 0.5))+
              theme(axis.title.y = element_text(size = rel(1), angle = 90))+
              theme(axis.title.x = element_text(size = rel(1), angle = 00))+
              labs(tag = "(c)")


            d<-rbind(d1, d2, d3)

            p4<- ggplot(d, aes(x=x, y=y)) +
              theme_bw()+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
              geom_vline(xintercept = times_sup1, size=.5, color="green",linetype="dashed")+
              geom_vline(xintercept = times_sup2, size=.5, color="blue",linetype="dashed")+
              geom_vline(xintercept = times_sup3, size=.5, color="black",linetype="dashed")+
              geom_line(data =subset(d1,x>times[1] & x<times_sup1), color = "green")+
              geom_line(data =subset(d2,x>times[1] & x<times_sup2), color = "blue")+
              geom_line(data =subset(d3,x>times[1] & x<times_sup3), color = "black")+
              labs(x = "Time",y="Temperature")+
              theme(plot.title = element_text(size=40))+
              theme(plot.title = element_text(hjust = 0.5))+
              theme(axis.title.y = element_text(size = rel(1), angle = 90))+
              theme(axis.title.x = element_text(size = rel(1), angle = 00))+
              labs(tag = "(d)")

            plot_grid(p1, p2,p3,p4)

          } else {

            stop("The available trends are associated with the RCP2.6 and RCP8.5 scenarios.")


          }

        }else{
          stop("The initial study temperature must be within the thermal tolerance range")
        }

      }else{

        stop("The minimum critical temperature must be less than the maximum critical temperature")
      }

    }else{

      stop("time_start must be less than time_end ")
    }

  }else{

    stop("The maximum simulation time is the year 2100 ")
  }

}




