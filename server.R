# Rcode_utils is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
# FluSubsample_Rshiny is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Rcode_utils.  
# If not, see <http://www.gnu.org/licenses/>.


# FluSubsample - server
# S J Lycett
# 17 Feb 2018
#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.


library(shiny)
library(ape)

############################################################
# compile custom functions
# only run once, e.g. read in data files etc

source("getEl.R")
source("calcDecimalDate.R")
source("birdSpecies_and_continents.R")
source("helper.R")

options(shiny.maxRequestSize=100*1024^2) 

# global variables

# filename for example sequences
defaultSeqFileName <- "example_data//example_seqs.fas"

# internal values for App
values <- reactiveValues()
values$seqs <- c()
values$taxa <- c()
values$tbl  <- c()
values$ntaxa<- 0
values$name <- ""
values$random_id <- ""
values$dataType <- ""
values$genome_tbl <- c()

subsampleValues <- reactiveValues()
subsampleValues$valid_time_inds  <- c()
subsampleValues$time             <- c()
subsampleValues$valid_space_inds <- c()
subsampleValues$space            <- c()
subsampleValues$valid_st_inds    <- c()
subsampleValues$subtype          <- c()
subsampleValues$valid_host_inds  <- c()
subsampleValues$host             <- c()
subsampleValues$valid_bird_inds  <- c()
subsampleValues$bird_host        <- c()
subsampleValues$joint_inds       <- c()
subsampleValues$jointTrait       <- c()
subsampleValues$ujointTrait      <- c()
subsampleValues$num_ujoint       <- 0
subsampleValues$nper             <- 0
subsampleValues$numSub           <- 0
subsampleValues$clickNumber      <- 1


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
   
  #######################################################################
  # Outputs
  
  #######################################################################
  # Well 1
  
  output$numSeqs <- renderText({
    paste("Total number of Sequences =",values$ntaxa)
  })
  
  output$segmentTable <- renderTable({
    cn      <- colnames(values$tbl)
    #Segment <- as.integer(values$tbl[,which(cn=="segNumber")])
    Segment <- values$tbl[,which(cn=="seg")]
    if (length(Segment) > 0) {
      stbl <- as.data.frame(t(as.matrix(table(Segment))))
    }
  })
  
  output$dataType <- renderText({
    paste("Data type =",values$dataType)
  })
  
  output$numWholeGenomes <- renderText({
    cn     <- colnames(values$tbl)
    if (length(cn) > 0) {
      inds   <- getSegInds(0,tbl=values$tbl,genome_tbl=values$genome_tbl)
      nwg    <- length(inds)
    } else {
      nwg    <- 0
    }
    paste("Number of whole genomes =",nwg)
  })
  
  output$dateContinentDistribution <- renderPlot({
    inds   <- getSegInds(as.integer(input$selectSegForDisplay),tbl=values$tbl,genome_tbl=values$genome_tbl)
    plot_dateContinentDistribution(values$tbl[inds,])
  })

  output$HNtypeBarPlot <- renderPlot({
    inds   <- getSegInds(as.integer(input$selectSegForDisplay),tbl=values$tbl,genome_tbl=values$genome_tbl)
    plot_HNtypeBarPlot( values$tbl[inds,] )
  })
  
  output$HostBarPlot <- renderPlot({
    inds   <- getSegInds(as.integer(input$selectSegForDisplay),tbl=values$tbl,genome_tbl=values$genome_tbl)
    plot_HostBarPlot( values$tbl[inds,] )
  })
  
  #################################################################################
  # Well 2
  
  output$timeSelectionDetails <- renderText({
    
    cn       <- colnames(values$tbl)
    if (length(cn) > 0) {
      decDates <- as.numeric(values$tbl[,which(cn=="decDate")])
      dateTxt  <- values$tbl[,which(cn=="dateTxt")]
      years    <- floor(decDates)
      ym       <- apply(as.matrix(dateTxt), 1, getEl, ex=3, sep="/", reconstruct=TRUE)
    
      if (input$jointTime=="Year") {
        time_pts  <- years
        time_inds <- which(is.finite(years))
      } else {
        time_pts <- ym
        time_inds<- setdiff(which(is.finite(years)), which(apply(as.matrix(ym), 1, nchar) < 6) )
      }
      
      utime <- sort(unique(time_pts[time_inds]))
      minTime<- utime[1]
      maxTime<- utime[length(utime)]
      
      subsampleValues$valid_time_inds <- time_inds
      subsampleValues$time <- time_pts
      
      paste("Number of unique time pts=",length(utime),"from",minTime,"to",maxTime)
    } else {
      paste("Number of unique time pts=0")
    }
    
  })
  
  output$spaceSelectionDetails <- renderText({
    cn <- colnames(values$tbl)
    if (length(cn) > 0) {
      
      if (input$jointSpace=="Continent") {
        space <- values$tbl[,which(cn=="continent")]
      } else if (input$jointSpace=="Region") {
        space <- values$tbl[,which(cn=="region")]
      } else if (input$jointSpace=="Country") {
        space <- values$tbl[,which(cn=="country")]
      } else if (input$jointSpace=="State") {
        bigCountry <- c("USA","Canada","China","Russia")
        space  <- values$tbl[,which(cn=="country")]
        country<- values$tbl[,which(cn=="country")] 
        place  <- values$tbl[,which(cn=="place")]
        for (j in 1:length(bigCountry)) {
          kk <- which(country==bigCountry[j])
          if (length(kk) > 0) {
            space[kk] <- place[kk]
          }
        }
      }
      
      
      space_inds <- which(space != "-")
      uspace <- sort(unique(space[space_inds]))
      subsampleValues$valid_space_inds <- space_inds
      subsampleValues$space <- space
      
      paste("Number of unique locations=",length(uspace))
      
    } else {
      paste("Number of unique locations=0")
    }
    
  })
  
  output$subtypeSelectionDetails <- renderText({
    cn <- colnames(values$tbl)
    if (length(cn) > 0) {
      hst     <- values$tbl[,which(cn=="Htype")]
      nst     <- values$tbl[,which(cn=="Ntype")]
      st      <- paste(hst,nst,sep="")
      
      if (input$jointSubtype=="All") {
        st_inds <- 1:length(st)
      } else {
        st_inds <- which((hst != "HX") & (nst != "NX"))
      }
      ust <- sort(unique(st[st_inds]))
      subsampleValues$valid_st_inds <- st_inds
      subsampleValues$subtype <- st
      
      paste("Number of unique subtypes=",length(ust))
        
    } else {
      paste("Number of unique subtypes=0")
    }
    
    
  })
  
  output$mainHostSelectionDetails <- renderText({
    cn <- colnames(values$tbl)
    if (length(cn) > 0) {
      host <- values$tbl[,which(cn=="host")]
      
      if (input$jointHost=="All") {
        host_inds <- 1:length(host)
      } else if (input$jointHost=="Known") {
        host_inds <- setdiff(1:length(host), c(which(host=="-"),which(host=="Environment")))
      } else {
        majorHosts <- c("Avian","Swine","Human","Equine","Canine")
        host_inds  <- c()
        for (j in 1:length(majorHosts)) {
          host_inds <- c(host_inds, which(host==majorHosts[j]))
        }
      }
      subsampleValues$valid_host_inds <- host_inds
      subsampleValues$host <- host
      
      uhost <- sort(unique(host[host_inds]))
      paste("Number of unique main hosts=",length(uhost))
      
    } else {
      paste("Number of unique main hosts=0")
    }
  })
  
  output$birdTable <- renderTable({
    cn <- colnames(values$tbl)
    if (length(cn) > 0) {
      selcol <- switch(input$jointAvian,
                      "Avian" = which(cn=="host"),
                      "BirdOrder" = which(cn=="birdOrder"),
                      "AvianHost3" = which(cn=="avHost3"),
                      "AvianHost4" = which(cn=="avHost4"),
                      "AvianHost5" = which(cn=="avHost5")
                     )
      bird_host <- values$tbl[,selcol]
      #bird_inds <- which(bird_host != "-")
      
      #subsampleValues$valid_host_inds <- bird_inds
      #subsampleValues$host            <- bird_host
      
      bird_inds <- subsampleValues$valid_host_inds
      if (input$jointAvian != "Avian") {
        bird_inds <- setdiff(bird_inds,c(which(bird_host=="-"),which(bird_host=="Avian")))
      }
      subsampleValues$valid_bird_inds <- bird_inds
      subsampleValues$bird_host       <- bird_host
      
      #inds <- getSegInds(as.integer(input$selectSegForDisplay),values$tbl)
      inds <- getSegInds(as.integer(input$selectSegForDisplay),tbl=values$tbl,genome_tbl=values$genome_tbl)
      inds <- intersect(bird_inds,inds)
      
      btbl <- as.data.frame(t(as.matrix(table(values$tbl[inds,selcol]))))
    }
  })
  
  output$birdSelectionDetails <- renderText({
    cn <- colnames(values$tbl)
    if (length(cn) > 0) {
      if (length(subsampleValues$valid_bird_inds) > 0) {
        bird_host <- subsampleValues$bird_host[subsampleValues$valid_bird_inds]
        ubirds    <- unique(bird_host)
        paste("Number of unique hosts=",length(ubirds))
      } else {
        paste("Number of unique hosts=0")
      }
    } else {
      paste("Number of unique hosts=0")
    }
  })
  
  #################################################################################
  # Well 3
  
  output$numJointCategories <- renderText({
    cn  <- colnames(values$tbl)
    if (length(cn) > 0) {
      time_inds <- subsampleValues$valid_time_inds
      space_inds<- subsampleValues$valid_space_inds
      st_inds   <- subsampleValues$valid_st_inds
      host_inds <- subsampleValues$valid_bird_inds
      joint_inds<- intersect(time_inds,space_inds)
      joint_inds<- intersect(joint_inds, st_inds)
      joint_inds<- intersect(joint_inds, host_inds)
      if (length(joint_inds) > 0) {
        jointTrait<- paste(subsampleValues$time, subsampleValues$space, subsampleValues$subtype, subsampleValues$bird_host)
        inds  <- getSegInds(as.integer(input$selectSegForDisplay),tbl=values$tbl,genome_tbl=values$genome_tbl)
        inds2 <- intersect(joint_inds,inds)
        #ujointTrait <- sort(unique(jointTrait[joint_inds]))
        ujointTrait <- sort(unique(jointTrait[inds2]))
        num_ujoint  <- length(ujointTrait)
      
        subsampleValues$joint_inds <- joint_inds
        subsampleValues$jointTrait <- jointTrait
        subsampleValues$ujointTrait<- ujointTrait
        subsampleValues$num_ujoint <- num_ujoint
        paste("Number of unique traits=",num_ujoint)
      } else {
        paste("Number of unique traits=0")
      }
    } else {
      paste("Number of unique traits=0")
    }
  })
  
  output$minMaxPerJointCategory <- renderText({
    if (subsampleValues$num_ujoint > 0) {
      
      #inds       <- getSegInds(as.integer(input$selectSegForDisplay),values$tbl)
      inds <- getSegInds(as.integer(input$selectSegForDisplay),tbl=values$tbl,genome_tbl=values$genome_tbl)
      
      joint_inds <- intersect(subsampleValues$joint_inds, inds)
      jointTrait <- subsampleValues$jointTrait[joint_inds]
      jtbl       <- table( jointTrait )
      
      min_per_cat<- min(jtbl)
      max_per_cat<- max(jtbl)
      av_per_cat <- mean(jtbl)
      paste("Per category: min=",min_per_cat,", max=",max_per_cat,", mean=",format(av_per_cat,digits=2),sep="")
    } else {
      paste("Per category: (no information)")
    }
  })
  
  output$jointCategoryBarPlot <- renderPlot({
    if (subsampleValues$num_ujoint > 0) {
      subsampleValues$nper <- input$nper
      #inds       <- getSegInds(as.integer(input$selectSegForDisplay),values$tbl)
      inds <- getSegInds(as.integer(input$selectSegForSub),tbl=values$tbl,genome_tbl=values$genome_tbl)
      
      joint_inds <- intersect(subsampleValues$joint_inds, inds)
      jointTrait <- subsampleValues$jointTrait[joint_inds]
      jtbl       <- table( jointTrait )
      sub_jtbl   <- jtbl
      kk <- which(sub_jtbl > subsampleValues$nper)
      sub_jtbl[kk] <- subsampleValues$nper
      numSub <- sum(sub_jtbl)
      subsampleValues$numSub <- numSub
      #print(paste("numSub=",numSub))
      
      #temp <- getSubsampledInds(subsampleValues$jointTrait, joint_inds, subsampleValues$nper)
      #print(paste("actual numSub=",length(temp)))
      
      barplot(jtbl,names=1:length(jtbl),main="Distribution of traits (before)")
      lines(c(0,length(jtbl)*1.25),c(subsampleValues$nper,subsampleValues$nper),col="red")
    }
  })
  
  output$numSubsampledSeqs <- renderText({
    paste("Number of subsampled sequences=",subsampleValues$numSub)
  })
  
  output$subsampleBarPlot <- renderPlot({
    if (subsampleValues$num_ujoint > 0) {
      inds       <- getSegInds(as.integer(input$selectSegForSub),tbl=values$tbl,genome_tbl=values$genome_tbl)
      joint_inds <- intersect(subsampleValues$joint_inds, inds)
      
      sub_inds <- getSubsampledInds(subsampleValues$jointTrait, joint_inds, subsampleValues$nper)
      #print(paste("actual numSub=",length(sub_inds)))
      
      #cn <- colnames(values$tbl)
      #print( table(values$tbl[sub_inds,which(cn=="seg")]) )
      #print( table(values$tbl[sub_inds,which(cn=="host")]))
      #print( table(values$tbl[sub_inds,which(cn=="region")]))
      
      sub_values <- switch(input$displaySubPlot,
                           "Time" = subsampleValues$time[sub_inds],
                           "Space" = subsampleValues$space[sub_inds],
                           "Subtype" = subsampleValues$subtype[sub_inds],
                           "Host" = subsampleValues$bird_host[sub_inds])
      plot_generalBarPlot(sub_values, title=paste("Distribution of ",input$displaySubPlot," (after)",sep=""))
    }
  })
  
  #  output$dateDistribution <- renderPlot({
  #    plot_dateDistribution(values$tbl)
  #  })
  
  #  output$HtypeBarPlot <- renderPlot({
  #    plot_HtypeBarPlot(values$tbl)
  #  })
  
  #  output$NtypeBarPlot <- renderPlot({
  #    plot_NtypeBarPlot(values$tbl)
  #  })
  

  ################################################################################################
  # reactives
  ################################################################################################
  
  ###########################################################
  # LOAD SEQUENCES AND PROCESS SEQUENCE NAMES
  # load the fasta format sequences
  # note the sequences names must be defined on the NCBI fasta def line like this:
  # >{serotype}_{host}_{accession}_{strain}_{country}_{year}/{month}/{day}_{segname}
  
  ##############################################
  # Reacting to file inputs
  observe({
   sname <- input$seqfile$datapath
   if (length(sname) > 0) {
    # read sequences
    seqs <- read.dna( sname, format="fasta", as.matrix=FALSE)
    taxa <- as.matrix(attributes(seqs)$names)
    
    # process sequence names
    tbl  <- taxaToTbl( taxa )                 # defined in helper.R
    tbl  <- cbind(tbl, getSeqLengths(seqs))   # add sequence lengths to tbl
    colnames(tbl)[1] <- "seqName"
    colnames(tbl)[length(tbl[1,])] <- "len"
    
    # determine dataType from segments (Whole Genome etc)
    dataType <- findDataType(tbl)
    
    genome_tbl <- getWholeGenomeTable(tbl)
   
    # update the internal values object
    values$seqs 	<- seqs
    values$taxa   <- taxa
    values$tbl		<- tbl
    values$ntaxa  <- length(taxa)
    values$name   <- input$seqfile$name
    values$dataType<- dataType
    values$genome_tbl <- genome_tbl
    values$random_id <- floor(runif(1, min=1000, max=1000000000000))
   }
  })
  
  ##############################################
  # reacting to click
  observeEvent(input$doSubSampling, {
    #session$sendCustomMessage(type = 'testmessage',
    #                          message = 'Thank you for clicking')
    
    
    if (subsampleValues$num_ujoint > 0) {
      config_id <- paste(values$random_id,"-",subsampleValues$clickNumber,sep="")
      
      print(paste("Subsampling clicked for id = ",config_id," data name = ",values$name))
      
      renamedSeqs <- renameSeqs(values$seqs, values$tbl)
    
      for (r in 1:input$nreps) {
        print(paste("Replicate",r))
        
        seg_selection <- as.integer(input$selectSegForSub)
        if (seg_selection > 0) {
          # just one segment
          inds       <- getSegInds(seg_selection,tbl=values$tbl,genome_tbl=values$genome_tbl)
          joint_inds <- intersect(subsampleValues$joint_inds, inds)
          sub_inds   <- getSubsampledInds(subsampleValues$jointTrait, joint_inds, subsampleValues$nper)
          rootname   <- paste( gsub("\\.","_",values$name),"_id",config_id,"_seg",seg_selection,sep="")
          writeSeqs(renamedSeqs, sub_inds, r, rootname=rootname)
          writeTraits(renamedSeqs, subsampleValues$time, subsampleValues$space, subsampleValues$subtype, subsampleValues$bird_host,
                      inds=sub_inds, rep=r, rootname=rootname)
          
        } else {
          # whole genome
          # for the default segment
          inds       <- getSegInds(seg_selection,tbl=values$tbl,genome_tbl=values$genome_tbl)
          joint_inds <- intersect(subsampleValues$joint_inds, inds)
          sub_inds   <- getSubsampledInds(subsampleValues$jointTrait, joint_inds, subsampleValues$nper)
          cn         <- colnames(values$tbl)
          selected_isNames <- values$tbl[sub_inds, which(cn=="isName")]
          wg_inds    <- match(selected_isNames, values$genome_tbl[,1])
          
          #print( which(!is.finite(wg_inds)) )
          
          cn <- colnames(values$genome_tbl)
          for (s in 1:8) {
            sub_inds   <- as.integer( values$genome_tbl[wg_inds,which(cn==paste("Index_",s,sep=""))] )
            rootname   <- paste( gsub("\\.","_",values$name),"_id",config_id,"_seg",s,sep="")
            writeSeqs(renamedSeqs, sub_inds, r, rootname=rootname)
            writeTraits(renamedSeqs, subsampleValues$time, subsampleValues$space, subsampleValues$subtype, subsampleValues$bird_host,
                        inds=sub_inds, rep=r, rootname=rootname)
          }
          
        }
        
      }
      
      rootname   <- paste( gsub("\\.","_",values$name),"_id",config_id,sep="")
      details    <- list(
      title="FluSubsample Rshiny App",
      processing_date=paste(Sys.time()),
      upload_id=values$random_id,
      click_id=subsampleValues$clickNumber,
      input_filename=values$name,
      dataType=values$dataType,
      ntaxa=values$ntaxa,
      nreps=input$nreps,
      nper=input$nper,
      segment_subsamped=input$selectSegForSub,
      time_category=input$jointTime,
      space_category=input$jointSpace,
      subtype_category=input$jointSubtype,
      host_processing=input$jointHost,
      avian_host_processing=input$jointAvian,
      number_categories=subsampleValues$num_ujoint,
      number_taxa_in_subsample=subsampleValues$numSub
      )
      #temp <- t(t((unlist(details))))
      #write.table(temp, row.names=TRUE, col.names=FALSE)
      
      writeSubsamplingDetails(details, rootname=rootname)
      writeFullTbl(values$tbl, rootname=rootname)
      writeGenomeTbl(values$genome_tbl, rootname=rootname)
      
      print("Done")
      
      subsampleValues$clickNumber <- as.integer(details$click_id)+1
      
    } else {
      print("Select some data first")
    }
    
    
  })
  
})
