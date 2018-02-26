# Rcode_utils is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
# FluSubsample_Rshiny is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Rcode_utils.  
# If not, see <http://www.gnu.org/licenses/>.

# FluSubsample - user interface
# S J Lycett
# 17 Feb 2018
#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#


library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("FluSubsampling"),
  
  # Output: Tabset w/ plot, summary, and table ----
  #tabsetPanel(type = "tabs",
  #            tabPanel("Plot", plotOutput("HNtypeBarPlot")),
  #            tabPanel("Summary", textOutput("dataType")),
  #            tabPanel("Table", tableOutput("segmentTable"))
  #),
  
  "Use Flu Fasta sequences from NCBI Influenza Virus Resource, FASTA Definition line MUST BE:",
  p(">{serotype}_{host}_{accession}_{strain}_{country}_{year}/{month}/{day}_{segname}", style = "color:magenta"),
  a(href="https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database", "Link for Individual Segments"),
  " or ",
  a(href="https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=genomeset", "Link for Complete Genomes"),
  br(),
  br(),
  
  wellPanel(
    # Input: Select fasta sequences file
    fileInput("seqfile", "Choose Fasta File",
              multiple = FALSE,
              accept = c(".fas",".fasta",".FAS",".FASTA")),
    
    fluidRow(
      column(4,
          p(textOutput("numSeqs")),
          p(textOutput("dataType")),
          p(textOutput("numWholeGenomes"))
      ),
      column(8,
             tableOutput("segmentTable")
             )
    ),
    
    fluidRow(
      column(12,
             radioButtons("selectSegForDisplay","Select segments to display",
                          choiceNames=c("WholeGenome","PB2","PB1","PA","HA","NP","NA","MP","NS","All"),
                          choiceValues=c(0,1,2,3,4,5,6,7,8,9),
                          inline=TRUE,
                          selected=4))
    ),
    
    fluidRow(
      column(4,
          plotOutput("dateContinentDistribution")
      ),
      column(4,
             plotOutput("HNtypeBarPlot")
             
      ),
      column(4,
             plotOutput("HostBarPlot")
      )
    )
  ),
  
#  wellPanel(
#    h4("Select categorisation scheme"),
#    fluidRow(
#      column(4,
#             radioButtons("jointTime","Select timescale",
#                          choiceNames=c("Year","Year-Month"),
#                          choiceValues=c("Year","Month"),
#                          selected="Year"),
             
             #renderText("timeSelectionDetails"),
             
#             radioButtons("jointSpace","Select spatial scale",
#                   choiceNames=c("Continent","Region","Country or State"),
#                   choiceValues=c("Continent","Region","Country"),
#                   selected="Country")
#      ),
#      column(4,
#             radioButtons("jointSubtype","Subtype choices",
#                          choiceNames=c("All","Exclude Mixed"),
#                          choiceValues=c("All","NoMix"),
#                          selected="NoMix"),
#             radioButtons("jointHost","Main Host choices",
#                          choiceNames=c("All","All known","Major only"),
#                          choiceValues=c("All","Known","Major"),
#                          selected="Known"
#                         )
#             ),
#      column(4,
#             
#             radioButtons("jointAvian","Split Avian host into",
#                          choiceNames=c("Avian","BirdOrder","3 Types","4 Types","5 Types"),
#                          choiceValues=c("Avian","BirdOrder","AvianHost3","AvianHost4","AvianHost5"),
#                          selected="AvianHost4")
#             
#             )
#      
#    ),
#    
#    hr(),
#    
#    h4("Summary of selection before subsampling"),
#    fluidRow(
#      column(12,
#             textOutput("timeSelectionDetails"),
#             textOutput("spaceSelectionDetails"),
#             textOutput("subtypeSelectionDetails"),
#             textOutput("mainHostSelectionDetails"),
#             textOutput("birdSelectionDetails")
#      )
#    ),
#    
#    h4("Detailed Host table before subsampling"),
#    fluidRow(
#      column(12,
#            tableOutput("birdTable")
#      )
#    )
#    
#    
#  ),
  
  wellPanel(
    h4("Select categorisation scheme"),
    fluidRow(
      column(3,
             radioButtons("jointTime","Select timescale",
                          choiceNames=c("Year","Year-Month"),
                          choiceValues=c("Year","Month"),
                          selected="Year")
             
             ),
      column(2,
             radioButtons("jointSpace","Select spatial scale",
                          choiceNames=c("Continent","Region","Country only","Country or State"),
                          choiceValues=c("Continent","Region","Country","State"),
                          selected="State")
             
             ),
      column(2,
             radioButtons("jointSubtype","Subtype choices",
                          choiceNames=c("All","Exclude Mixed"),
                          choiceValues=c("All","NoMix"),
                          selected="NoMix")
             
             ),
      column(2,
             radioButtons("jointHost","Main Host choices",
                          choiceNames=c("All","All known","Major only"),
                          choiceValues=c("All","Known","Major"),
                          selected="Known")
             
             ),
      column(3,
             radioButtons("jointAvian","Split Avian host into",
                          choiceNames=c("Avian","BirdOrder","3 Types","4 Types","5 Types"),
                          choiceValues=c("Avian","BirdOrder","AvianHost3","AvianHost4","AvianHost5"),
                          selected="AvianHost4")
             
             )
      
    ),
    
    hr(),
    
    fluidRow(
      column(3, textOutput("timeSelectionDetails")),
      column(2, textOutput("spaceSelectionDetails")),
      column(2, textOutput("subtypeSelectionDetails")),
      column(2, textOutput("mainHostSelectionDetails")),
      column(3, textOutput("birdSelectionDetails"))
    ),
    
    hr(),
    
    h4("Detailed Host table before subsampling"),
    fluidRow(
      column(12,
             tableOutput("birdTable")
      )
    )
    
    
  ),

  wellPanel(
    h4("Subsampling control"),
    fluidRow(
      column(4,
        sliderInput("nper", "Number per category", min=0, max=20, value=3),
        sliderInput("nreps", "Number of replicates", min=1, max=10, value=1),
        radioButtons("selectSegForSub","Select segment for subsampling",
                     choiceNames=c("WholeGenome","PB2","PB1","PA","HA","NP","NA","MP","NS"),
                     choiceValues=c(0,1,2,3,4,5,6,7,8),
                     inline=TRUE,
                     selected=0),
        actionButton("doSubSampling","Do subsampling")
      ),
      
      column(4,
             plotOutput("jointCategoryBarPlot"),
             textOutput("numJointCategories"),
             textOutput("numSubsampledSeqs"),
             textOutput("minMaxPerJointCategory")#,
             #plotOutput("jointCategoryBarPlot")
             
             ),
      column(4,
             plotOutput("subsampleBarPlot"),
             radioButtons("displaySubPlot","Select category to display",
                    choices=c("Time","Space","Subtype","Host"),
                    inline=TRUE, selected="Time")#,
             #plotOutput("subsampleBarPlot")
             )
      )#,
    
    #fluidRow(
    #  column(4,
    #         radioButtons("selectSegForSub","Select segment for subsampling",
    #                      choiceNames=c("WholeGenome","PB2","PB1","PA","HA","NP","NA","MP","NS"),
    #                      choiceValues=c(0,1,2,3,4,5,6,7,8),
    #                      inline=TRUE,
    #                      selected=0)
    #         ),
    #  column(4,
    #         plotOutput("jointCategoryBarPlot")
    #         ),
    #  column(4,
    #         plotOutput("subsampleBarPlot")
    #         )
    #)
    
    
    
  )
  
  
  
)

)
