library(shiny)
library(ggplot2)
library(reshape2)
library(Rcpp)
library(gridExtra)
library(data.table)
library(deSolve)
library(driftSim)
library(plyr)

saveDir <- "~/tmp/"
options(shiny.maxRequestSize=1000*1024^2)

shinyServer(
    function(inputs, output, session){
        
        values <- reactiveValues()
        
        create_startingK <- observeEvent(inputs$iniCalc,{
            print("Generating hostK...")
            N <- isolate(inputs$s0) + isolate(inputs$i0) +  isolate(inputs$r0)
            if(is.null(isolate(inputs$hostInput))){
                inputFile <- "hosts.csv"
            }
            else {
                inputFile <- isolate(inputs$hostInput$datapath)
            }
            iniK <- generateHostKDist(inputFile,N)
            write.table(iniK, file=paste(saveDir, "outputs/iniK.csv",sep=""),row.names=FALSE, col.names=FALSE, sep=",")
        })
        
        plot_SIR <- function(filename){
            N <- isolate(inputs$s0) + isolate(inputs$i0) +  isolate(inputs$r0) + 500
            dat <- read.csv(filename,header=0)
            dat <- cbind(seq(1,nrow(dat),by=1),dat)
            colnames(dat) <- c("t","S","I","R")
            dat <- melt(dat,id="t")
            colnames(dat) <- c("t","Population","value")
            SIR_plot <- ggplot() +
                geom_line(data=dat,aes(x=t,y=value,colour=Population,group=Population)) +
                xlab("Time (days)") +
                ylab("Number Individuals") +
                scale_y_continuous(limits=c(0,N),expand=c(0,0))+
                scale_x_continuous(limits=c(0,inputs$dur+1),expand=c(0,0))+
                theme(
                    text=element_text(colour="gray20",size=14),
                    plot.title=element_text(size=28),
                    legend.text=element_text(size=20,colour="gray20"),
                    legend.title=element_text(size=20,colour="gray20"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x = element_line(colour = "gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    axis.text.x=element_text(colour="gray20"),
                    panel.background=element_blank(),
                    axis.text.y=element_text(colour="gray20"))
            return(SIR_plot)
        }

        output$incidenceText1 <- renderText({
            N <- isolate(inputs$s0) + isolate(inputs$i0) +  isolate(inputs$r0)
            paste(round((values$infections/N)*100,3),"%",sep="")
        })

        
        
        calculate_deltaVMat <- observeEvent(inputs$dVcalc,{
            print("Calculating deltaV matrix...")
            maxV = 3
            maxK = 80
            time_step = 1
            p = as.numeric(inputs$p) 
            r = as.numeric(inputs$r)
            q = as.numeric(inputs$q)
            a =as.numeric(inputs$a)
            b = as.numeric(inputs$b)
            kc = as.numeric(inputs$kc)
            pars <- c(p,r,b,a,kc,q)
            
            difeq <- function(t, V, params){
                x <- V
                j <- params[1]
                p <- params[2]
                r <- params[3]
                b <- params[4]
                a <- params[5]
                kc <- params[6]
                q <- params[7]
                
                immK <- r*j
                if(immK <0) immK <- 0
                
                f_x <- (1-exp(-p*(x+q)))^(immK)
                g_x <- exp(-a*(x^b))
                f_dx <- p*(immK)*((1-exp(-p*(V+q)))^(immK-1))*(exp(-p*(x+q)))
                g_dx <- -a*b*x^(b-1)*exp(-a*(x^b))
                
                dV <- f_x*g_dx + g_x*f_dx
                return(list(dV))
            }
            V = seq(0,maxV,by=0.01)
            immKs = seq(0,maxK,by=0.1)
            allV <- matrix(nrow=length(immKs),ncol=length(V))
            for(j in 1:length(immKs)){
                print(j)
                for(i in 1:length(V)){
                    deltaV <- ode(y=c(V[i]),seq(0,time_step,by=1/40),difeq,c(immKs[j],pars))
                    allV[j,i] <- deltaV[length(deltaV)] - V[i]
                }
            }
            write.table(allV, file=paste(saveDir,"outputs/deltaVMat.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",")
        })
        
        run_sim <- observeEvent(inputs$run, {
            print("Running simulation")
            #' Make sure Rcpp will compile and run
               
            SIR_flag <- 1 %in% inputs$flags #' Flag to save SIR dynamics
            voutput1_flag <- 2 %in% inputs$flags #' Flag to save virus information for Sean's phylogenetic tree
            voutput2_flag <- 3 %in% inputs$flags #' Flag to save pairwise distance matrix
            time_flag <- 4 %in% inputs$flags #' Flag to record time taken for simulation
            VERBOSE <- 9 %in% inputs$flags #' Outputs in simulation
            save_state <- 5 %in% inputs$flags #' Flag to save the final state of the simulation
            input_flagA <- 6 %in% inputs$flags #' Flag to use specified file as input for simulation
            input_flagB <- 7 %in% inputs$flags #' Flag to use specified file as input for simulation
            input_flag <- input_flagA || input_flagB
            save_k <- 8 %in% inputs$flags
            
            flags <- c(SIR_flag, voutput1_flag, voutput2_flag, time_flag, save_state, input_flagA, input_flagB, save_k)
            flags <- as.numeric(flags)
            
            S0 <- as.numeric(inputs$s0)
            I0 <- as.numeric(inputs$i0)
            R0 <- as.numeric(inputs$r0)

            contactRate <- as.numeric(inputs$contact)
            mu <- 1/(as.numeric(inputs$mu)*365)
            wane <- 1/as.numeric(inputs$wane)
            gamma <- 1/as.numeric(inputs$gamma)
            iniBind <- as.numeric(inputs$iniBinding)
            meanBoost = as.numeric(inputs$boost)
            maxTitre = as.numeric(inputs$maxTitre)
            iniDist = as.numeric(inputs$iniDist)
            saveFreq = as.numeric(inputs$kSaveFreq)
            
            hostpars <- c(S0,I0, R0,contactRate,mu,wane,gamma,iniBind, meanBoost, iniDist,saveFreq,maxTitre)

            deltaVMat <- unname(as.matrix(read.csv(paste(saveDir,"outputs/deltaVMat.csv",sep=""),header=FALSE)))

            p = as.numeric(inputs$p) 
            r = as.numeric(inputs$r)
            q = as.numeric(inputs$q)
            a =as.numeric(inputs$a)
            b = as.numeric(inputs$b)
            n = as.numeric(inputs$n)
            v = as.numeric(inputs$v)

            probMut = inputs$probMut
            expDist = inputs$expDist
            kc = inputs$kc
            VtoD = inputs$VtoD

            viruspars <- c(p,r,q,a,b,n,v,probMut,expDist,kc,VtoD)
            
            progress_within<- shiny::Progress$new()
            on.exit(progress_within$close())
            callback <- function(x) {
                progress_within$set(value=x[[1]]/inputs$dur,detail= x[[1]])
            }
            if(is.null(inputs$hostInput)){
                inputFiles <- c("hosts.csv")
            }
            else {
                inputFiles <- c(inputs$hostInput$datapath)
            }
            print("Working directory: " )
            print(saveDir)

            iniK <- NULL
            if(input_flagA){
                print(paste("Generating starting K from... ",inputFiles[1],sep=""))
                iniK <- generateHostKDist(inputFiles[1],S0+I0+R0)
            }
            if(input_flagB){
                print(paste("Using saved starting K from ",inputFiles[1],sep=""))
                iniK <- read.csv(inputFiles[1],header=FALSE)[,1]
            }

            if(length(inputs$scenarios) > 0){
                infections <- NULL
                withProgress(message="Simulation number", value=0, detail=1, {
                    index <- 1
                    for(i in inputs$scenarios){
                        print(paste("Scenario number: ",as.numeric(i),sep=""))
                        progress_within$set(message = "Day", value = 0)
                        Sys.sleep(0.1)
                        filename1 <- paste("scenario_",i,"_SIR.csv",sep="")
                        filename2 <- paste("voutput1_",i,".csv",sep="")
                        filename3 <- paste("voutput2_",i,".csv",sep="")
                        filename4 <- paste("hosts_",i,".csv",sep="")
                        filename5 <- paste("hostKs_",i,".csv",sep="")
                        filenames <- c(filename1, filename2, filename3, filename4, filename5)
                        infections[[index]] <- run_simulation(flags,hostpars,viruspars,deltaVMat,iniK,0,inputs$dur,inputFiles,filenames,VERBOSE, as.numeric(i),callback)
                        print(paste("Number of new viruses: ", infections[[index]], sep=""))
                        incProgress(1/length(inputs$scenarios),detail=(as.numeric(i)+1))
                        index <- index + 1
                        for(j in filenames){
                            if(file.exists(j)) file.rename(from=j,to = paste(saveDir,"/outputs/",j,sep=""))
                        }
                    }                   
                })
                values$infections <- infections
            }
            print("Saving simulation parameters used...")
            vParNames <- c("p","r","q","a","b","n","v","probMut","expDist","kc","VtoD")
            hParNames <- c("s0","i0","r0","c","mu","w","g","iniBind","meanBoost","iniDist","saveFreq","maxTitre")
            allNames <- c(vParNames, hParNames)
            allPars <- c(viruspars,hostpars)
            names(allPars) <- allNames
            write.table(allPars,file=paste(saveDir,"outputs/parameters.csv",sep=""))
        })

        output$sim_main_1<- renderPlot({
            inputs$run
            if(1 %in% inputs$scenarios){
                g <- plot_SIR(paste(saveDir,"outputs/scenario_1_SIR.csv",sep=""))
                g
            } else{
                NULL
            }
        })
        output$sim_main_2<- renderPlot({
            inputs$run
            if(2 %in% inputs$scenarios){
                g <- plot_SIR(paste(saveDir,"outputs/scenario_2_SIR.csv",sep=""))
                g
            }else{
                NULL
            }
        })
        output$sim_main_3<- renderPlot({
            inputs$run
            if(3 %in% inputs$scenarios){
                g <- plot_SIR(paste(saveDir,"outputs/scenario_3_SIR.csv",sep=""))
                g
            } else{
                NULL
            }
        })
        output$sim_main_4<- renderPlot({
            inputs$run
            if(4 %in% inputs$scenarios){
                g <- plot_SIR(paste(saveDir,"outputs/scenario_4_SIR.csv",sep=""))
                g
            } else{
                NULL
            }
        })
        output$parameterPlots <- renderPlot({
            base <- ggplot(data.frame(x=c(0,inputs$N_reinfect)),aes(x)) + stat_function(fun=dexp,geom="line",colour="red",args=list(rate=inputs$expDist)) +
             xlab("Size of mutation") +
                ylab("Probability (given a mutation did occur)") +
                scale_y_continuous(expand=c(0,0))+
                scale_x_continuous(expand=c(0,0))+
                theme(
                    text=element_text(colour="gray20",size=14),
                    plot.title=element_text(size=28),
                    legend.text=element_text(size=20,colour="gray20"),
                    legend.title=element_text(size=20,colour="gray20"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x = element_line(colour = "gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    axis.text.x=element_text(colour="gray20"),
                    panel.background=element_blank(),
                    axis.text.y=element_text(colour="gray20"))
            base
        })
        output$boostPlot <- renderPlot({
            base <- ggplot(data.frame(x=c(0,3*inputs$boost)),aes(x)) + stat_function(fun=dpois,geom="bar",colour="black",n=(3*inputs$boost + 1),args=list(lambda=inputs$boost))+
                xlab("Magnitude of boost following infection") +
                ylab("Probability") +
                scale_y_continuous(expand=c(0,0))+
                scale_x_continuous(expand=c(0,0))+
                theme(
                    text=element_text(colour="gray20",size=14),
                    plot.title=element_text(size=28),
                    legend.text=element_text(size=20,colour="gray20"),
                    legend.title=element_text(size=20,colour="gray20"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x = element_line(colour = "gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    axis.text.x=element_text(colour="gray20"),
                    panel.background=element_blank(),
                    axis.text.y=element_text(colour="gray20"))
            base
        })
        output$dlPlots <- downloadHandler(
            filename = "allSIR.png",
            content = function(file){
                ggsave(file,plot=plotAllSIR(),device="png",width=8)
            }
        )
        output$dlPlot1 <- downloadHandler(
            filename = "scenario1SIR.png",
            content = function(file){
                ggsave(file,plot=plot_SIR(paste(saveDir,"outputs/scenario_1_SIR.csv",sep="")),device="png", width=10,height=6)
            }
        )
        output$dlPlot2 <- downloadHandler(
            filename = "scenario2SIR.png",
            content = function(file){
                 ggsave(file,plot=plot_SIR(paste(saveDir,"outputs/scenario_2_SIR.csv",sep="")),device="png", width=10,height=6)
            }
        )
        output$dlPlot3 <- downloadHandler(
            filename = "scenario3SIR.png",
            content = function(file){
                ggsave(file,plot=plot_SIR(paste(saveDir,"outputs/scenario4_3_SIR.csv",sep="")),device="png", width=10,height=6)
            }
        )
        output$dlPlot4 <- downloadHandler(
            filename = "scenario4SIR.png",
            content = function(file){
                ggsave(file,plot=plot_SIR(paste(saveDir,"outputs/scenario_4_SIR.csv",sep="")),device="png", width=10,height=6)
            }
        )

        
        output$dlParPlots <- downloadHandler(
            filename = "bindingAvidPlots.png",
            content = function(file){
                ggsave(file=file,plotDynamics(),device="png")
            }
        )
        output$dlOutputsSIR <- downloadHandler(
            filename = function(){paste("outputs",".tar",sep="")},
            content=function(file){
                tar(file,paste(saveDir,"outputs",sep=""))
            }
        )
        
        output$dlPars <- downloadHandler(
            filename = "parameters.csv",
            content = function(file){
                p = as.numeric(inputs$p) 
                r = as.numeric(inputs$r) 
                b = as.numeric(inputs$b)
                a =as.numeric(inputs$a) 
                n = as.numeric(inputs$n)
                v = as.numeric(inputs$v)
                q = as.numeric(inputs$q)
                N_reinfect = as.numeric(inputs$N_reinfect)
                max_reinfect = N_reinfect
                delta = as.numeric(inputs$delta)
                allPars <- list("p"=p,"r"=r,"b"=b,"a"=a,"n"=n,"v"=v,"q"=q,"max_titre"=N_reinfect,"delta"=delta)
                write.csv(allPars, file,row.names=FALSE)
            }
        )
        
        plotAllSIR <- function(){
                p1 <- plot_SIR(paste(saveDir,"outputs/scenario_1_SIR.csv",sep=""))
                p2 <- plot_SIR(paste(saveDir,"outputs/scenario_2_SIR.csv",sep=""))
                p3 <- plot_SIR(paste(saveDir,"outputs/scenario_3_SIR.csv",sep=""))
                p4 <- plot_SIR(paste(saveDir,"outputs/scenario_4_SIR.csv",sep=""))
                plot.list <- list(p1, p2, p3, p4, ncol=1)
                do.call(arrangeGrob, plot.list)
        }
        
        plotDynamics <- function(){
            p = as.numeric(inputs$p) #' parameter to control degree by which changes in binding avidity affect probability of escape from immune response
            r = as.numeric(inputs$r) #' parameter to control degree by which previous exposure reduce probability of immune escape
            b = as.numeric(inputs$b) #' parameter to control the shape of the relationship between probability of successful replication and changes in binding avidity
            a =as.numeric(inputs$a) #' controls rate of changes of relationship between probability of successful replication and change in binding avidity.
            n = as.numeric(inputs$n) #' average number of virus copies: n is number of offspring per virus replication
            v = as.numeric(inputs$v) #' v is number of virions initialyl transmitted
            nv = n*v
            q = as.numeric(inputs$q) #' parameter to control the shape of the relationship between binding avidity and immune escape (shift on the x-axis)
            V = seq(0, 2, by = 0.005) #' Binding avidity
            N_reinfect = as.numeric(inputs$N_reinfect)
            max_reinfect = N_reinfect
            delta = as.numeric(inputs$delta)
            

            #' Get a colour scale gradient blue to red

            f_array <- NULL #' Survival prob
            df_array <- NULL #'  Derivative of survival prob
            for(k in 0:(N_reinfect-1)){
                probTarget = exp(-p*(V+q)) #' probability of being targetted by immune system. As binding avidity increases, this probability decreases
                probEscape = 1-probTarget #' probability of escape
                immK = r*(k- delta) #' Strength of host immune respone. As k increases, virus must escape more antibodies. As delta increases, this effect diminishes as infecting virus is further from host immunity.
                if(immK < 0) immK = 0
                
                f = (probEscape)^(immK) #' probability of escape from immunity
                
                if(k >= 1) f_dash= immK*p*((1-probTarget)^(immK-1))*(probTarget) #' derivative of this. ie. rate of change of relationship between binding avidity and probability of immune escape
                else f_dash= rep(0, length(V))
                f_array[[k+1]] <- f  
                df_array[[k+1]] <- f_dash
            }

            
            probInf <- exp(-a*(V^b)) #' probability of successful infection within a host. as binding avidity increases in a naive host, chance of successfully replicating decreases
            probInf_dash= -a*b*(V^(b-1))*exp(-a*(V^b)) #' rate of change of this relationship

            rho_array <- NULL
            dV_array <- NULL
            probRep_array <- NULL
            for(i in 1:length(df_array)){
                R0 = f_array[[i]]*probInf*n
                rho = 1 - R0^-v
                rho[rho < 0] = 0
                rho_array[[i]] <- rho
                dV = df_array[[i]]*probInf + f_array[[i]]*probInf_dash
                dV_array[[i]] <- dV
                
                probReplication = f_array[[i]]*probInf
                probRep_array[[i]] = probReplication
            }

            rho0 = max(rho_array[[1]])
            rho1 = max(rho_array[[2]])
            rho2 = max(rho_array[[3]])

            probRep_data <- NULL
            for(i in 1:length(probRep_array)){
                data <- data.frame(x=V,y=probRep_array[[i]],z=as.character(i))
                probRep_data <- rbind(probRep_data,data)
            }

            colours <- NULL
            blue <- rev(seq(0,1,by=1/N_reinfect))
            red <- seq(0,1,by=1/N_reinfect)
            for(i in 1:N_reinfect){
                colours[i] <- rgb((i-1)/(max_reinfect-1),0,(N_reinfect-i)/(max_reinfect-1))
            }

            #' Replication
            A <- ggplot() + geom_line(data=probRep_data,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=14),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
                ylab("Probability of Successful Replication Within a Host (theta(V))") + 
                xlab("Binding Avidity")


            rho_data<- NULL
            for(i in 1:length(rho_array)){
                data <- data.frame(x=V,y=rho_array[[i]],z=as.character(i))
                rho_data<- rbind(rho_data,data)
            }

            #' Infection (Rho)
            B <- ggplot() + geom_line(data=rho_data,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=14),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
                ylab("Probability of Infection Between Hosts (rho)") + 
                xlab("Binding Avidity") +
                geom_hline(yintercept = rho0,linetype="longdash",colour="dodgerblue4")+
                geom_hline(yintercept = rho1,linetype="longdash",colour="dodgerblue4")+
                geom_hline(yintercept = rho2,linetype="longdash",colour="dodgerblue4")


            df_data<- NULL
            for(i in 1:length(df_array)){
                data <- data.frame(x=V,y=df_array[[i]],z=as.character(i))
                df_data<- rbind(df_data,data)
            }

            #' Derivative of f
            C <- ggplot() + geom_line(data=df_data,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=14),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0)) +
                ylab(expression(d*f/d*V)) + 
                xlab("Binding Avidity")


            dV_data<- NULL
            for(i in 1:length(dV_array)){
                data <- data.frame(x=V,y=dV_array[[i]],z=as.character(i))
                dV_data<- rbind(dV_data,data)
            }

            #' Infection
            D <- ggplot() + geom_line(data=dV_data,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=12),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0)) +
                ylab(expression(d*beta/d*V)) + 
                xlab("Binding Avidity") +
                geom_hline(yintercept=0,linetype='longdash',colour="gray20")


            f_data<- NULL
            for(i in 1:length(f_array)){
                data <- data.frame(x=V,y=f_array[[i]],z=as.character(i))
                f_data <- rbind(f_data,data)
            }

            E <- ggplot() + geom_line(data=f_data,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=14),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0)) +
                ylab("Probability of Evading Immune System, f(k,V)") + 
                xlab("Binding Avidity")

            probRep_dat <- data.frame(x=V,y=probInf)
            F <- ggplot() + geom_line(data=probRep_dat,aes(x=x,y=y,colour="red")) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=14),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0)) +
                ylab("Probability of Successful Replication, g(V)") + 
                xlab("Binding Avidity")

            #' Plot antigenic distance against j at a given binding avidity.
            d <- seq(0,N_reinfect,by=0.1)
            fixedV <- inputs$bindAvid
            delta_array <- NULL
            for(k in 0:(N_reinfect-1)){
                probT = exp(-p*(fixedV+q))
                probE = 1- probT
                immK1 = r*(k-d)
                immK1[immK1 < 0] <- 0
                probSurvival1 = 1 - (n*exp(-a*(fixedV^b))*(probE^immK1))^-v
                probSurvival1[probSurvival1 < 0] <- 0
                delta_array[[k+1]] <- probSurvival1
            }
            delta_data<- NULL
            for(i in 1:length(delta_array)){
                data <- data.frame(x=d,y=delta_array[[i]],z=as.character(i))
                delta_data <- rbind(delta_data,data)
            }
            G <- ggplot() + geom_line(data=delta_data,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=14),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
                ylab("Probability of Infection Between Hosts") + 
                xlab("Antigenic Distance to Host Immunity")

            immK_array <- NULL
            immK2 <- seq(0,r*N_reinfect,by=1)
            for(k in 0:(N_reinfect-1)){
                probT1 = exp(-p*(fixedV+q))
                probE1 = 1- probT1
                probSurvival2 = 1 - (n*exp(-a*(fixedV^b))*(probE1^immK2))^-v
                probSurvival2[probSurvival2 < 0] <- 0
                immK_array[[k+1]] <- probSurvival2
            }
            immK_data<- NULL
            for(i in 1:length(immK_array)){
                data <- data.frame(x=immK2/r,y=immK_array[[i]],z=as.character(i))
                immK_data <- rbind(immK_data,data)
            }
            H <- ggplot() + geom_line(data=immK_data,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=14),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
                ylab("Probability of Infection Between Hosts") + 
                xlab("ImmK")
            
            g <- grid.arrange(A,B,E,F,C,D,G,H,ncol=2)
        }
        output$Main <- renderPlot({
            plotDynamics()
        },height=1000,width=1000)

        output$antigenicDistanceTime <- renderPlot({
            cutOff <- 3
            noSamples <- 1000
            scenario <- inputs$scenarioPlot
            filename <- paste("outputs/voutput1_",scenario,".csv",sep="")
            if(file.exists(filename)){
                dat <- fread(filename,header=T,data.table=FALSE)
                tmp <- dat[sample(nrow(dat),noSamples,replace=FALSE),c("birth","distRoot","distance_to_parent")]
                tmp$class[tmp$distance_to_parent >= cutOff] <- paste("Antigenic change >= ", cutOff,sep="")
                tmp$class[tmp$distance_to_parent < cutOff] <- paste("Antigenic change <",cutOff,sep="")
                p <- ggplot(data=tmp,aes(x=birth,y=distRoot,colour=class,group=class)) + geom_point() +
                    theme(
                       # panel.grid.major=element_blank(),
                       # panel.grid.minor=element_blank(),
                        #panel.background=element_blank(),
                        axis.text.x=element_text(colour="gray20",size=12),
                        axis.text.y = element_text(colour="gray20",size=12),
                        text = element_text(colour="gray20",size=14),
                        axis.line=element_line(colour="gray20"),
                        axis.line.x =element_line(colour="gray20"),
                        axis.line.y=element_line(colour="gray20"),
                                        # legend.position = "none",
                        axis.title=element_text(size=12)
                    ) +
                    ylab("Antigenic Distance to Root") + 
                    xlab("Day of birth")
            }
            else {
                p <- NULL
            }
            return(p)
        })
                
        output$hostImmunityHist <- renderPlot({

        })
        output$immKTime <- renderPlot({
            
        })
        output$virusPairwiseDist <- renderPlot({

        })
        output$deltaRBP <- renderPlot({

        })

    })
