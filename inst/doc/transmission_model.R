## ----setup, include=FALSE------------------------------------------------
library(driftSim)
library(ggplot2)
knitr::opts_chunk$set(echo = TRUE, fig.width=7,fig.height=6,
                      message=FALSE, results=FALSE,warning=FALSE)

## ----fig.width = 4, fig.height=3-----------------------------------------
 base <- ggplot(data.frame(x=c(0,3*3)),aes(x)) + stat_function(fun=dpois,geom="bar",colour="black",n=(3*3 + 1),args=list(lambda=3))+
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

## ----fig.width=7,fig.height=5--------------------------------------------
attach(exampleParameters)
 base <- ggplot(data.frame(x=c(0,10)),aes(x)) + stat_function(fun=dexp,geom="line",colour="red",args=list(rate=1/viruspars["expDist"])) +
             xlab("Size of mutation") +
                ylab("Probability (given a mutation did occur)") +
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

