library(shiny)
library(ggplot2)
library(gridExtra)


shinyUI(
    navbarPage("Antigenic Drift Simulation",
               #' Parameter exploration panel
            
               #' Simulation exploration panel
               tabPanel("Simulation",
                        sidebarPanel(
                            fluidRow(
                                column(4,actionButton("run",(strong("Run Simulation")))),
                                column(4,actionButton("dVcalc",(strong("Recalcuate deltaV")))),
                                column(4,actionButton("iniCalc",(strong("Create iniK"))))
                            ),
                            checkboxGroupInput("flags",
                                               label=h3(strong("Flags")),
                                               choices=list(
                                                   "Save SIR output"=1,
                                                   "Save virus characteristics"=2,
                                                   "Save virus pairwise distances"=3,
                                                   "Record run time"=4,
                                                   "Save final state"=5,
                                                   "Use input starting conditions (new)"=6,
                                                   "Use input starting conditions (saved)"=7,
                                                   "Save hostK"=8,
                                                   "Verbose"=9
                                               ),
                                               selected=1),
                            h3(strong("Host Population Parameters")),
                          
                            fluidRow(
                                column(4,
                                       numericInput("s0", "S0", 999900)
                                       ),
                                column(4,
                                numericInput("i0", "I0", 100)
                               ),
                               column(4,
                                      numericInput("r0", "R0", 0)
                                      )
                            ),
                            fluidRow(
                                column(4, numericInput("dur","Duration",400)),
                                column(4,numericInput("contact", "Contact Rate", 0.7)),
                                column(4,numericInput("mu", "Birth/Death Rate", 70))
                            ),
                            fluidRow(
                                column(4,numericInput("wane", "Waning Rate", 25)),
                                column(4,numericInput("gamma", "Recovery Time", 3.3)),
                                column(4,numericInput("iniBinding", "Ini Binding Avidity ", 0.67))
                            ),
                            fluidRow(
                                column(4,numericInput("boost","Mean Boost",6)),
                                column(4, numericInput("iniDist","Ini Distance", 0)),
                                column(4, numericInput("kSaveFreq","K Save Freq",5))
                            ),
                            fluidRow(
                                column(4,numericInput("maxTitre","Maximum Titre",10))
                            ),
                            h3(strong("Virus Parameters")),
                            fluidRow(
                                column(6,numericInput("probMut","Probability of a Mutation",0.1)),
                                column(6,numericInput("expDist","Shape of mutation distribution",1))
                            ),
                            fluidRow(
                                column(6,numericInput("kc","Rate of binding avidity change",0.5)),
                                column(6,numericInput("VtoD","Impact of dV on delta",0))
                            ),
                            checkboxGroupInput(
                                "scenarios",
                                h3(strong("Scenarios")),
                                choices=list(
                                    "Random drift; fixed binding avidity" = 1,
                                    "No drift; adaptive binding avidity; adaptive antigenic change" = 2,
                                    "Random drift; adaptive binding avidity; adaptive antigenic change" = 3,
                                    "Random drift; adaptive binding avidity; no adaptive antigenic change" = 4
                                ),
                                selected=1
                            ),
                            h3(strong("Filenames")),
                            fluidRow(
                                column(6,fileInput("hostInput","Host input file"))
                            ),
                            h3(strong("Actions")),
                            
                            fluidRow(
                                column(4,downloadButton("dlPlot1","Download plot 1")),
                                column(4,downloadButton("dlPlot2","Download plot 2")),
                                column(4,downloadButton("dlPlot3","Download plot 3"))
                            ),
                            fluidRow(
                                column(4,downloadButton("dlPlot4","Download plot 4")),
                                column(4,downloadButton("dlPlots","Download all plots")),
                                column(4,downloadButton("dlOutputsSIR","Download outputs"))
                            )
                        ),
                        mainPanel(
                            fluidRow(
                                column(12,
                                       h4(strong("Incidence for each scenario:")),
                                       textOutput("incidenceText1")
                                       )
                            ),
                            fluidRow(
                                column(12,
                                       h4(strong("Scenario 1: Random drift; fixed binding avidity")),
                                       plotOutput("sim_main_1")
                                       ),
                                column(12,
                                       h4(strong("Scenario 2: No drift; adaptive binding avidity; adaptive antigenic change")),
                                       plotOutput("sim_main_2")
                                       )
                            ),
                            fluidRow(
                                column(12,
                                       h4(strong("Scenario 3: Random drift; adaptive binding avidity; adaptive antigenic change")),
                                       plotOutput("sim_main_3")
                                       ),
                                column(12,
                                       h4(strong("Scenario 4: Random drift; adaptive binding avidity; no adaptive antigenic change")),
                                       plotOutput("sim_main_4")
                                       )
                            ),
                            fluidRow(
                                column(6,
                                       h4(strong("a) Mutation probability distribution, given that a mutation occured")),
                                       plotOutput("parameterPlots")
                                       ),
                                column(6,
                                       h4(strong("b) Immune boosting following recovery")),
                                       plotOutput("boostPlot")
                                       )
                            ))
                            
                        ),
               tabPanel("Parameters",
                        sidebarPanel(
                            h3(strong("Input Parameters")),
                            br(),
                            h4(strong("Probability of Survival Parameters (f)")),
                            numericInput("p",label="p - change in binding avidity affecting probability of immune escape:",value=4),
                            numericInput("q",label="q - shift on x-axis:",value=1),
                            numericInput("r",label="r - impact of previous exposure on immunity:",value=70),
                            numericInput("delta",label="delta - distance between virus and host immunity:",value=0),
                            br(),
                            h4(strong("Probability of Within Host Infection Parameters (g)")),
                            numericInput("b",label="b - shape of V/survival relationship:",value=3),
                            numericInput("a",label="a - rate of change of V/survival relationship:",value=0.7),
                            br(),
                            h4(strong("Other Parameters")),
                            numericInput("n",label="n - virus offspring per replication:",value=4),
                            numericInput("v",label="v - number of virions transmitted:",value=1),
                            numericInput("N_reinfect",label="Max titre value:",value=30),
                            numericInput("bindAvid",label="Point binding avidity:",value=0.8),
                            br(),
                            h4(strong("Downloads")),
                            fluidRow(
                                downloadButton("dlParPlots","Download plots"),
                                downloadButton("dlPars","Download parameters")
                            )
                        ),
                        
                        mainPanel(
                            plotOutput("Main")
                        )
                        ),
               tabPanel("Outputs",
                        sidebarPanel(
                            h3(strong("Simulation Output Plots")),
                            selectInput("scenarioPlot","Choose a scenario to plot",
                                        c("Scenario 1"=1,
                                          "Scenario 2"=2,
                                          "Scenario 3" =3,
                                          "Scenario 4"=4
                                          ),
                                        selected=1
                                        ),
                            actionButton("updateOutputPlots","Update plots")
                        ),
                        mainPanel(
                            plotOutput("antigenicDistanceTime")
                                        #         plotOutput("hostImmunityHist"),
                                        #        plotOutput("immKTime"),
                                        #        plotOutput("virusPairwiseDist"),
                                        #        plotOutput("deltaRBP")                            
                            )),
               
               tabPanel("Phylogeny")
               )
)




