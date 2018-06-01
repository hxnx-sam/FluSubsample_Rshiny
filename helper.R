# Rcode_utils is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
# FluSubsample_Rshiny is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Rcode_utils.  
# If not, see <http://www.gnu.org/licenses/>.

# helper function for FluSubsample server
# S J Lycett
# 17 Feb 2018
# 1 June 2018 - nper = 0 means dont do any subsampling

# load other R functions (Sam standard functions)
#source("getEl.R")
#source("calcDecimalDate.R")
#source("birdSpecies_and_continents.R")

# helper function for taxaToTbl
sequenceLength	<- function( seq ) {
  #return( length(as.character(seq)[[1]]) )
  return( length(unlist(seq)) )
}

# helper function to classify avian sequences
classifyAvian <- function( tbl ) {
  cn    <- colnames(tbl)
  host  <- tbl[,which(cn=="host")]
  host2 <- tbl[,which(cn=="host2")]
  
  # correction for condor
  jj <- which(host2=="condor")
  if (length(jj) >0) {
    host[jj] <- "Avian"
  }
  
  #jj <- which(host=="Avian" & host2=="environment")
  #if (length(jj) >0) {
  #  print(tbl[jj,1])
  #}
  
  ainds     <- which(host=="Avian")
  birdOrder <- host
  wildDom          <- array("-",length(host))
  shortLong        <- array("-",length(host))
  
  if (length(ainds) > 0) {
    birdOrder[ainds] <- birdClass( host2[ainds] )
    wildDom[ainds]   <- getWildDomestic(host2[ainds])
    wildDom          <- gsub("Domestic","Dom",wildDom)
    shortLong[ainds] <- getLongShort(host2[ainds])
  }
  
    domans_inds         <- which((wildDom=="Dom") & (birdOrder=="ans"))
    domgal_inds         <- which((wildDom=="Dom") & (birdOrder=="gal"))
    wild_inds           <- which((wildDom=="Wild"))
    wildans_inds        <- which((wildDom=="Wild") & (birdOrder=="ans"))
    wildother_inds      <- which((wildDom=="Wild") & (birdOrder != "ans"))
    wildlongans_inds    <- which((wildDom=="Wild") & (birdOrder=="ans") & shortLong=="Long")
    wildshortans_inds   <- which((wildDom=="Wild") & (birdOrder=="ans") & shortLong=="Short")
  
    avHost3 <- host
    if (length(domans_inds) > 0) { avHost3[domans_inds] <- "Dom-ans" }
    if (length(domgal_inds) > 0) { avHost3[domgal_inds] <- "Dom-gal" }
    if (length(wild_inds) > 0)   { avHost3[wild_inds]   <- "Wild-bird" }
  
    avHost4 <- avHost3
    if (length(wildans_inds) > 0)   { avHost4[wildans_inds]  <- "Wild-ans" }
    if (length(wildother_inds) > 0) { avHost4[wildother_inds]<- "Wild-other" }
  
    avHost5 <- avHost4
    if (length(wildlongans_inds) > 0)  { avHost5[wildlongans_inds]  <- "Wild-ans-long"  }
    if (length(wildshortans_inds) > 0) { avHost5[wildshortans_inds] <- "Wild-ans-short" }

  
  #print(table(birdOrder,avHost5))
  jj <- which(birdOrder=="-")
  if (length(jj)>0) {
    print("Unknown birds")
    print(table(host2[jj]))
  }

    
  avtbl <- cbind(birdOrder,wildDom,shortLong,avHost3,avHost4,avHost5)
  

  
  return(avtbl)
  
}

taxaToBeastNames <- function(tbl) {
    cn      <- colnames(tbl)
    subtype <- paste(tbl[,which(cn=="Htype")],tbl[,which(cn=="Ntype")],sep="")
    host    <- tbl[,which(cn=="host")]
    accn    <- tbl[,which(cn=="accn")]
    isName  <- tbl[,which(cn=="isName")]
    isName  <- gsub(" ","_",isName)
    isName  <- gsub("\\(","-",isName)
    isName  <- gsub("\\)","-",isName)
    isName  <- gsub("\\.","-",isName)
    isName  <- gsub("'","",isName)
    isName  <- gsub(",","",isName)
    decDates<- as.numeric(tbl[,which(cn=="decDate")])
    beastNames <- paste(paste(subtype,host,accn,isName,sep="_"),format(decDates,digits=7),sep="|")
    return( beastNames )
}

###########################################################
# PROCESS SEQUENCE NAMES TO TABLE

# note the sequences names (taxa) must be defined on the NCBI fasta def line like this:
# >{serotype}_{host}_{accession}_{strain}_{country}_{year}/{month}/{day}_{segname}
taxaToTbl <- function( taxa ) {
  
  # break the sequence names to component parts
  subtype <- apply(taxa, 1, getEl, ind=1, sep="_")
  host    <- apply(taxa, 1, getEl, ind=2, sep="_")
  accn	  <- apply(taxa, 1, getEl, ind=3, sep="_")

  # do this way
  isName	<- apply(taxa, 1, getEl, ex=c(1,2,3), sep="_", reconstruct=TRUE)
  isName	<- apply(as.matrix(isName), 1, getEl, ex=c(1,2,3), sep="_", reconstruct=TRUE, fromEnd=TRUE)

  country <- apply(taxa, 1, getEl, ind=3, sep="_", fromEnd=TRUE)
  dateTxt <- apply(taxa, 1, getEl, ind=2, sep="_", fromEnd=TRUE)

  seg	  	<- apply(taxa, 1, getEl, ind=1, sep="_", fromEnd=TRUE)
  # seg is of form "4 (HA)"
  segNumber 	<- as.integer( apply(as.matrix(seg), 1, getEl, ind=1, sep=" ") )

  # calculate decimal date
  decDate <- apply(as.matrix(dateTxt), 1, calcDecimalDate_fromTxt)

  # classify to UN GEO SCHEME world regions
  region <- unGeoScheme(country)
  
  # summarise to continents
  continent <- array("-",length(region))
  ucontinent<- c("Africa","Asia","Europe","America","Oceania","Antarctia")
  for (j in 1:length(ucontinent)) {
    jj <- grep(ucontinent[j],region)
    if (length(jj) > 0) {
      continent[jj] <- ucontinent[j]
    }
  }
  #print(table(continent))
  jj <- which(continent=="-")
  if (length(jj) > 0) {
    print("Unknown location")
    print(table(country[jj]))
    print(taxa[jj])
  }
  

  # extract place from isolate name (might not be very good)
  place2 <- apply(as.matrix(isName), 1, getEl, ind=2, sep="/")
  place3 <- apply(as.matrix(isName), 1, getEl, ind=3, sep="/")
  
  place  <- place3
  hinds <- which(host=="Human")
  if (length(hinds) > 0) {
    place[hinds] <- place2[hinds]
  }
  
  # more detailed host processing
  jj <- which( (host=="") | (host==" "))
  if (length(jj)>0) {
    host[jj] <- "-"
  }
  host2 <- host
  ainds <- which(host != "Human")
  if (length(ainds) > 0) {
    host2[ainds] <- place2[ainds]
  }
  
  # more detailed subtype processing
  st    <- as.matrix(subtype)
  Htype <- apply(st, 1, getEl, ind=1, sep="N")
  Ntype <- paste("N",apply(st, 1, getEl, ind=2, sep="N"),sep="")
  
  
  #print(table(st))
  #print(table(Htype))
  #print(table(Ntype))
  
  ok_inds <- grep("H[1-9][0-9]?",Htype)
  x_inds  <- setdiff(1:length(Htype),ok_inds)
  if (length(x_inds) > 0) {
    Htype[x_inds] <- "HX"
  }
  Htype <- gsub("mixed","",Htype)
  Htype <- gsub(",","",Htype)
  Htype <- gsub(" ","",Htype)
  
  
  ok_inds <- grep("N[1-9][0-9]?",Ntype)
  x_inds  <- setdiff(1:length(Ntype),ok_inds)
  if (length(x_inds) > 0) {
    Ntype[x_inds] <- "NX"
  }
  
  #print(table(st,Htype))
  #print(table(st,Ntype))

  # make info table
  info_tbl <- cbind(taxa, subtype, Htype, Ntype, host, host2, accn, isName, place, country, region, continent, dateTxt, decDate, seg, segNumber)
  
  # add the avian sequences classification
  av_tbl   <- classifyAvian(info_tbl)
  info_tbl <- cbind(info_tbl,av_tbl)
  
  # add the beast names
  beastNames <- taxaToBeastNames(info_tbl)
  info_tbl <- cbind(info_tbl,beastNames)
  
  return( info_tbl )
}

getSeqLengths <- function( seqs ) {
  len   <- unlist(lapply(seqs, sequenceLength)) 
  return( len )
}

findDataType <- function( tbl ) {
  cn <- colnames(tbl)
  segtbl <- table(as.integer(tbl[,which(cn=="segNumber")]))
  minSeg <- min(segtbl)
  maxSeg <- max(segtbl)
  avSeg  <- mean(segtbl)
  if ( (avSeg/maxSeg > 0.9) & (minSeg > 0) ) {
    dataType <- "Whole Genome"
  } else {
    jj <- which(segtbl != 0)
    if (length(jj) == 1) {
      dataType <- paste("Segment",j,sep=" ")
    } else {
      if (  all(segtbl > 0) ) {
        dataType <- "Segments 1-8"
      } else {
        dataType <- "Selection"
      }
    }
  }
  return( dataType ) 
}



getWholeGenomeTable <- function( tbl ) {
  cn <- colnames(tbl)
  isName  <- tbl[,which(cn=="isName")]
  segs    <- as.integer(tbl[,which(cn=="segNumber")])
  accn    <- tbl[,which(cn=="accn")]
  len     <- as.integer(tbl[,which(cn=="len")])
  uisName <- unique(isName)
  genome_numSegs <- array(0, length(uisName))
  genome_inds    <- matrix(0, length(uisName), 8)
  genome_accns   <- matrix(0, length(uisName), 8)
  
  for (j in 1:length(uisName)) {
    pos <- which(isName==uisName[j])
    temp<- sort(segs[pos])
    genome_numSegs[j] <- length(temp)
    if ((genome_numSegs[j]==8) & (all(temp==1:genome_numSegs[j]))) {
      ordered_inds <- pos[match(1:8, segs[pos])]
      genome_inds[j,] <- ordered_inds
      genome_accns[j,]<- accn[ordered_inds]
    } else {
      for (s in 1:8) {
        kk <- which(segs[pos]==s)
        
        if (length(kk)>1) {
          kk <- kk[which.max(len[pos[kk]])]
        }
        
        if (length(kk)==1) {
          genome_inds[j,s] <- pos[kk]
          genome_accns[j,s]<- accn[pos[kk]]
        }
      }
      genome_numSegs[j] <- length(which(genome_inds[j,]>0))
    }
    
  }
  
  #all( isName[genome_inds[1:100,8]]==uisName[1:100] )
  genome_tbl <- cbind(uisName,genome_accns,genome_inds,genome_numSegs)
  colnames(genome_tbl) <- c("isName",paste("Accn_",1:8,sep=""), paste("Index_",1:8,sep=""), "NumSegs")
  return( genome_tbl )
  
}

getWholeGenomeInds <- function(segNumber=4, genome_tbl=genome_tbl) {
  cn      <- colnames(genome_tbl)
  wg_inds <- which( genome_tbl[,which(cn=="NumSegs")] == 8 )
  seg_cols<- match(paste("Index_",1:8,sep=""), cn)
  if ((segNumber > 0) & (segNumber < 9)) {
    segInds <- as.integer( genome_tbl[wg_inds,seg_cols[segNumber]] )
  } else {
    segInds <- as.integer( genome_tbl[wg_inds,seg_cols[4]] )
  }
  return( segInds )
}

getSegInds <- function( segNumber=1, tbl=tbl, genome_tbl=c(-1)) {
  if ((segNumber > 0) & (segNumber < 9)) {
    cn   <- colnames(tbl)
    segs <- as.integer(tbl[,which(cn=="segNumber")])
    inds <- which(segs==segNumber)
    return(inds)
  } else {
    if ( (segNumber==0) & (length(genome_tbl)>2) ) {
      return( getWholeGenomeInds(genome_tbl=genome_tbl) )
    } else {
      return ( 1:length(tbl[,1]) )
    }
  }
}

getSubsampledInds <- function(jointTrait, joint_inds, nper) {
  if (nper == 0) {
    print("Nper=0 - not doing any subsampling")
    return ( joint_inds )
  } else {
    ujointTrait <- unique(jointTrait[joint_inds])
    num_uj      <- length(ujointTrait)
    selinds     <- c()
    for (j in 1:num_uj) {
      kk <- which(jointTrait[joint_inds]==ujointTrait[j])
      if (length(kk) > 0) {
        kk <- joint_inds[kk]
        if (length(kk) > nper) {
          kk <- sample(kk, nper)
        }
        selinds <- c(selinds,kk)
      }
    }
    return( selinds )
  }
}

renameSeqs <- function(seqs=seqs, tbl=tbl) {
  taxa <- as.matrix(attributes(seqs)$names)
  cn   <- colnames(tbl)
  tbl_taxa <- tbl[,which(cn=="seqName")]
  
  minds<- match(taxa, tbl_taxa)
  jj   <- which(is.finite(minds))
  print(length(jj))
  
  if (length(jj) >= (length(taxa)/2) ) {
    if ( all( taxa[jj] == tbl_taxa[minds[jj]] ) ) {
      newTaxa     <- taxa
      newTaxa[jj] <- tbl[minds[jj],which(cn=="beastNames")]
      attributes(seqs)$names <- newTaxa
    } else {
      print("Warning sequence names do not correspond, cant do renaming")
    }
  } else {
    print(paste("Warning not enough sequences match to do renaming - only",length(jj),"match"))
    print(taxa[1:10])
    print(tbl_taxa[1:10])
  }
  return( seqs )
}

# sequences output (local to server)
writeSeqs <- function(seqs, inds=inds, rep=1, localPath="userGenerated//", rootname="seqs") {
  sname <- paste(localPath,rootname,"_rep",rep,".fas",sep="")
  write.dna( seqs[inds], file=sname, format="fasta", nbcol=-1, colsep="")
  print(paste("Written",sname))
  return( sname )
}

# traits table output (local to server)
writeTraits <- function( seqs, time, space, subtype, bird_host, inds=inds, rep=1, localPath="userGenerated//", rootname="seqs") {
  tname <- paste(localPath,rootname,"_rep",rep,"_traits.txt",sep="")
  taxa  <- attributes(seqs)$names
  traitsTbl <- cbind(taxa, time, space, subtype, bird_host)
  traitsTbl <- traitsTbl[inds,]
  colnames(traitsTbl) <- c("traits","Year","Location","Subtype","Host")
  write.table(traitsTbl, file=tname, quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")
  print(paste("Written",tname))
  return( tname )
}

# settings output
writeSubsamplingDetails <- function( details, localPath="userGenerated//", rootname="seqs") {
  sname <- paste(localPath,rootname,"_settings.txt",sep="")
  temp <- t(t((unlist(details))))
  write.table(temp, row.names=TRUE, col.names=FALSE, file=sname, quote=FALSE, sep="=")
  print(paste("Written",sname))
  return( sname )
}

# full table output
writeFullTbl <- function( tbl, localPath="userGenerated//", rootname="seqs") {
  tname <- paste(localPath,rootname,"_info_tbl.csv",sep="")
  write.csv(tbl, file=tname, row.names=FALSE)
  print(paste("Written",tname))
  return( tname )
}

# genome table output
writeGenomeTbl <- function(genome_tbl, localPath="userGenerated//", rootname="seqs") {
  tname <- paste(localPath,rootname,"_genome_tbl.csv",sep="")
  write.csv(genome_tbl, file=tname)
  print(paste("Written",tname))
  return( tname )
}

###########################################################################
# plotting (used in server.R)

plot_dateDistribution <- function( tbl ) {
  cn       <- colnames(tbl)
  decDates <- as.numeric(tbl[,which(cn=="decDate")])

  if (length(decDates) > 0) {
    years    <- floor(decDates)
    minDate  <- min(years)
    maxDate  <- max(years)
    if ( (maxDate-minDate) < 3) {
      breaks <- seq(minDate,(maxDate+1),1/12)
    } else {
      breaks <- seq(minDate,(maxDate+1),1)
    }
    hist(decDates,breaks=breaks,xlab="Date",main="Date Distribution",col=hsv(0.6,0.8,0.8))
  }
}

plot_dateContinentDistribution <- function( tbl ) {
  cn       <- colnames(tbl)
  decDates <- as.numeric(tbl[,which(cn=="decDate")])
  
  ucontinent<- c("Asia","Africa","Europe","America","Oceania","Antarctia","-")
  ccols     <- c(hsv( (0:4)/5, 0.8, 0.8),"black","grey80")
  continents<- factor(tbl[,which(cn=="continent")], levels=ucontinent)
  
  if (length(decDates) > 0) {
    years    <- floor(decDates)
    jj <- which(is.finite(years))
    minDate  <- min(years[jj])
    maxDate  <- max(years[jj])
    
    years    <- factor(years[jj], levels=minDate:maxDate)
    yctbl    <- table(continents[jj],years[jj])
    barplot(yctbl,col=ccols,xlab="Year",main="Date Distribution")
    ctbl <- table(continents)
    ctxt <- paste(rownames(ctbl),"=",ctbl)
    legend("topleft",ctxt,pch=22,pt.bg=ccols,bty="n")
    
    print(paste("Min date=",min(decDates[jj])))
    print(paste("Max date=",max(decDates[jj])))
    print(ctbl)
  }
}

plot_HtypeBarPlot <- function( tbl ) {
  cn       <- colnames(tbl)
  Htype    <- tbl[,which(cn=="Htype")]
  if (length(Htype) > 0) {
    htbl <- table(factor(Htype,levels=c(paste("H",1:17,sep=""),"HX")))
    hcols<- c(hsv( (0:16)/17, 0.8, 0.8), "grey80")
    barplot(htbl,main="H subtype",col=hcols)
    print(htbl)
  }
}


plot_NtypeBarPlot <- function( tbl ) {
  cn       <- colnames(tbl)
  Ntype    <- tbl[,which(cn=="Ntype")]
  if (length(Ntype) > 0) {
    ntbl <- table(factor(Ntype,levels=c(paste("N",1:11,sep=""),"NX")))
    ncols<- c(hsv( (0:10)/11, 0.8, 0.8), "grey80")
    barplot(ntbl,main="N subtype",col=ncols)
    print(ntbl)
  }
}


plot_HNtypeBarPlot <- function( tbl ) {
  cn       <- colnames(tbl)
  Htype    <- tbl[,which(cn=="Htype")]
  Ntype    <- tbl[,which(cn=="Ntype")]
  if (length(Htype) > 0) {
    htbl <- table(factor(Htype,levels=c(paste("H",1:17,sep=""),"HX")))
    hcols<- c(hsv( (0:16)/17, 0.8, 0.8), "grey80")
    max_h<- paste(rownames(htbl)[which.max(htbl)],"=",format(100*(max(htbl)/sum(htbl)),digits=2),"%",sep="")
    
    ntbl <- table(factor(Ntype,levels=c(paste("N",1:11,sep=""),"NX")))
    ncols<- c(hsv( (0:10)/11, 0.8, 0.8), "grey80")
    max_n<- paste(rownames(ntbl)[which.max(ntbl)],"=",format(100*(max(ntbl)/sum(ntbl)),digits=2),"%",sep="")
    
    op <- par(mfrow=c(2,1))
    barplot(htbl,main="H subtype",col=hcols)
    legend("topright",max_h,pch=NA,bty="n")
    
    barplot(ntbl,main="N subtype",col=ncols)
    legend("topright",max_n,pch=NA,bty="n")
    par(op)
    
    print(htbl)
    print(ntbl)
  }
}

plot_HostBarPlot <- function( tbl ) {
  cn       <- colnames(tbl)
  host     <- tbl[,which(cn=="host")]
  if (length(host) > 0) {
    host_tbl <- table(host)
    nhost    <- length(host_tbl)
    
    host_cols<- c(hsv( (0:(nhost-1)/nhost), 0.8, 0.8), "grey80")
    max_host<- paste(rownames(host_tbl)[which.max(host_tbl)],"=",format(100*(max(host_tbl)/sum(host_tbl)),digits=2),"%",sep="")
    
    if (nhost==1) {
      pie(host_tbl,main="Host Species",col=host_cols)
      legend("top",max_host,pch=NA,bty="n")
    } else {
      barplot(host_tbl,main="Host Species",col=host_cols,horiz=TRUE)
      #legend("topright",max_host,pch=NA,bty="n")
      
      htxt <- paste(rownames(host_tbl),"=",host_tbl)
      legend("topright",htxt,pch=22,pt.bg=host_cols,bty="n")
    }
    
    print(host_tbl)
  }
}

plot_generalBarPlot <- function( values, titleTxt="", legpos="topright" ) {
  vtbl <- table(values)
  nu   <- length(vtbl)
  vcols<- hsv( 0:(nu-1)/nu, 0.8, 0.8)
  barplot(vtbl, col=vcols)
  if (titleTxt != "") {
    title(titleTxt)
  }
  if (legpos != "x") {
    if (length(vtbl) < 20) {
      legend(legpos, rownames(vtbl), pch=22, pt.bg=vcols, bty="n")
    }
  }
}