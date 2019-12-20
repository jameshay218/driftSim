## ----setup, include=FALSE------------------------------------------------
library(driftSim)
library(ggplot2)
attach(exampleParameters)
knitr::opts_chunk$set(echo = FALSE, fig.width=7,fig.height=6,
                      message=FALSE, results=FALSE,warning=FALSE)
p = viruspars["p"] #' parameter to control degree by which changes in binding avidity affect probability of escape from immune response
r = viruspars["r"] #' parameter to control degree by which previous exposure reduce probability of immune escape
b = viruspars["b"] #' parameter to control the shape of the relationship between probability of successful replication and changes in binding avidity
a =viruspars["a"] #' controls rate of changes of relationship between probability of successful replication and change in binding avidity.
n = viruspars["n"] #' average number of virus copies: n is number of offspring per virus replication
v = viruspars["v"] #' v is number of virions initialyl transmitted
nv = n*v
q = viruspars["q"] #' parameter to control the shape of the relationship between binding avidity and immune escape (shift on the x-axis)
V = seq(0, 3, by = 0.005) #' Binding avidity
N_reinfect = as.numeric(25)
max_reinfect = N_reinfect
delta = as.numeric(5)


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
A <- ggplot() + geom_line(data=probRep_data,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
    theme(plot.title=element_text(hjust=0.5),
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
    ylab("Probability of Successful\n Replication Within a Host (g(V))") + 
    xlab("Binding Avidity") + ggtitle("Plot A")

 rho_data<- NULL
for(i in 1:length(rho_array)){
    data <- data.frame(x=V,y=rho_array[[i]],z=as.character(i))
    rho_data<- rbind(rho_data,data)
}
B <- ggplot() + geom_line(data=rho_data,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(plot.title=element_text(hjust=0.5),
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
                ggtitle("Plot B") +
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
                theme(plot.title=element_text(hjust=0.5),
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
                xlab("Binding Avidity") + ggtitle("Plot C")


            dV_data<- NULL
            for(i in 1:length(dV_array)){
                data <- data.frame(x=V,y=dV_array[[i]],z=as.character(i))
                dV_data<- rbind(dV_data,data)
            }

            #' Infection
            D <- ggplot() + geom_line(data=dV_data,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(plot.title=element_text(hjust=0.5),
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
                geom_hline(yintercept=0,linetype='longdash',colour="gray20") + ggtitle("Plot D")


            f_data<- NULL
            for(i in 1:length(f_array)){
                data <- data.frame(x=V,y=f_array[[i]],z=as.character(i))
                f_data <- rbind(f_data,data)
            }

            E <- ggplot() + geom_line(data=f_data,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(plot.title=element_text(hjust=0.5),
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
                xlab("Binding Avidity")  + ggtitle("Plot E")

            probRep_dat <- data.frame(x=V,y=probInf)
            F <- ggplot() + geom_line(data=probRep_dat,aes(x=x,y=y,colour="red")) +
                theme(plot.title=element_text(hjust=0.5),
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
                xlab("Binding Avidity")  + ggtitle("Plot F")
 
            #' Plot antigenic distance against j at a given binding avidity.
            d <- seq(0,N_reinfect,by=0.1)
            fixedV <- 0.5
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
                    plot.title=element_text(hjust=0.5),
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
                xlab("Antigenic Distance to Host Immunity")  + ggtitle("Plot G")
 
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
                theme(plot.title=element_text(hjust=0.5),
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
                xlab("ImmK")  + ggtitle("Plot H")


## ---- fig.width=5,fig.height=3.5, fig.align="center"---------------------
E + ggtitle("Plot A")

## ---- fig.width=5,fig.height=3.5, fig.align="center"---------------------
C + ggtitle("Plot B")

## ---- fig.width=8,fig.height=4, fig.align="center"-----------------------
cowplot::plot_grid(A + ggtitle("Plot C"), F + ggtitle("Plot D"),ncol=2)

## ---- fig.width=5,fig.height=3.5, fig.align="center"---------------------
B + ggtitle("Plot E")

## ---- fig.width=5,fig.height=3.5, fig.align="center"---------------------
D + ggtitle("Plot F")

