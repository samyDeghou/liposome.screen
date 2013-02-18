#!/usr/bin/r





loadSecondaryStructure <- function(protein,method="psi.pred"){
  return(read.)
  
}

loadProtRunMap <- function(){
  return(read.table("/Users/deghou/remote/projects/liposomes/data/proteins_restored/prot.run.mapping"))
}


getProteinConcentration <- function(run){
 f <- read.table("/Users/deghou/remote/projects/liposomes/data/proteins_restored/protein.concentration",sep="\t")
 return(as.numeric(subset(f$V2, subset = f$V1 == run)))
}



loadSecondaryStructurePrediction <- function(protein){
  fileName <- paste('/Users/deghou/remote/projects/liposomes/data/proteins_restored/sequence/predicted.secondary.structure/Jnet/from.hack/',protein,sep="")
  f = readChar(fileName, file.info(fileName)$size)
  structure = strsplit(f,fixed = FALSE, split="\n| ")[[1]]
  protein <- list()
  protein["name"] = structure[1]
  protein["sequence"] = structure[2]
  protein["structure"] = structure[3]
  return(protein)
}



loadProtToWorkWith <- function(){
  return(as.character(read.table("/Users/deghou/remote/projects/liposomes/newdata/motif.search/to.work.with")$V1))
}

loadMultipleAlignment <- function(file="",method){
#  if(file == ""){
#    return(seqinr::read.alignment("/Users/deghou/remote/projects/liposomes/data/proteins_restored/sequence/fasta.JUST.PH.DOMAINS.ALIGNMENT.MUSCLE.txt",format="fasta"))    
# }else{
#    return(seqinr::read.alignment(paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/grouping.auto.generated/",file,sep=""),format = "fasta"))
#  }
  dir = "/Users/deghou/remote/projects/liposomes/newdata/alignments/"
  if(method == "mafft"){
    return(seqinr::read.alignment(paste(dir,"mafft.fa",sep=""),format="fasta"))
  }else if(method == "clustal"){
    return(seqinr::read.alignment(paste(dir,"clustalw.fa",sep=""),format="fasta"))
  }else if(method == "muscle"){
    return(seqinr::read.alignment(paste(dir,"muscle.fa",sep=""),format="fasta"))
  }else if(method == "multialin"){
    return(seqinr::read.alignment(paste(dir,"multialin.fa",sep=""),format="fasta"))
  }else if(method == "tcoffee"){
    return(seqinr::read.alignment(paste(dir,"tcoffee.fa",sep=""),format="fasta"))
  }
  
  
}

loadDataAsMatrixRespectingOrder <- function(file){
  f <- read.table(file,sep="\t")
  rnames <- unique(as.character(f$V1))
  cnames <- unique(as.character(f$V2))
  mat <- matrix(data=NA,nrow=length(rnames),ncol=length(cnames),dimnames=list(rnames,cnames))
  l <- nrow(f)
  for(k in 1:l){
    r <- as.character(f$V1[k])
    c <- as.character(f$V2[k])
    v <- as.numeric(as.character(f$V3[k]))
    mat[r,][c] <- v
  }
  return(mat)
}


loadGO <- function(what){
  if(what == "localisation"){
    return(read.table("/Users/deghou/remote/projects/liposomes/newdata/data.on.the.04.02.2013/enrichment/prot.localisation",sep="\t"))
  }else if(what == "function"){
    return(read.table("/Users/deghou/remote/projects/liposomes/newdata/data.on.the.04.02.2013/enrichment/prot.function",sep="\t"))
  }else if(what == "process"){
    return(read.table("/Users/deghou/remote/projects/liposomes/newdata/data.on.the.04.02.2013/enrichment/prot.process",sep="\t"))
  }
}

minusLipids <- function(mat,lip){
  lip.hom <- loadHomogeneousLipids()
  lip.con <- loadControlsLipids()
  lip.het <- loadHeterogeneousLipids()
  if(lip == "hom"){
    return(reduceMatrix(mat,lip.hom$V1,col=TRUE))
  }else if(lip == "het"){
    return(reduceMatrix(mat,lip.het$V1,col=TRUE))
  }else if(lip == "con"){
    return(reduceMatrix(mat,lip.con$V1,col=TRUE))
  }
  return(mat.x)
}

writeAlignmentsFiles <- function(){
  setwd("/Users/deghou/remote/projects/liposomes/newdata/data.on.the.04.02.2013/motif.search/")
  files <- list.files()
  files <- subset(files,grepl(pattern="proteins.sequence*",files))
  for(file in files){
    print(file)
    seqaln(aln=read.fasta(file),exepath="/Users/deghou/utils/",file=paste("/Users/deghou/remote/projects/liposomes/newdata/data.on.the.04.02.2013/motif.search/alignment.",file,sep=""))
  }
}



doMultipleAlignment <- function(){
  setwd("/Users/deghou/remote/projects/liposomes/newdata/data.on.the.04.02.2013/motif.search/")
  files <- list.files()
  for(file in files){
    print(file)
    proteins <- as.character(read.table(file)$V1)
    sequences <- sapply(proteins,getSequence)
#    sequences <- sequences[complete.cases(sequences)]
    print(paste(length(proteins),length(sequences)))
  }
}



loadProteinSequences <- function(){
  protein.sequences <- read.fasta("/Users/deghou/remote/projects/liposomes/data/proteins_restored/sequence/FULL.PROT.UNIQUE.SEQUENCE",as.string=TRUE)
  return(protein.sequences)
}

stringenceToRevisit <- function(){
  file <- read.table("/Users/deghou/remote/projects/liposomes/to.revisit.unconsistent.together.tolerance.2")
  proteins <- unique(as.character(file$V4))
  lipids <- unique(as.character(file$V5))
  for(protein in proteins){
    for(lipid in lipids){
      pl <- subset(file,subset=(file$V4 == protein & file$V5 == lipid))
      if(max(pl$V6) < 0.01){
        write.table(pl,file="/Users/deghou/remote/projects/liposomes/to.revisit.unconsistent.together.tolerance.2.ALL.UNDER.0.01",append=TRUE)
      } 
    }
  }
}

doMatrixOutOfScreening <- function(just.three=TRUE,including.controls=FALSE,values="sr",method="mean"){
  
#  data <- loadNewData(including.controls=including.controls,just.three=just.three)
#  data = data.dec
#  data = data.dec.1
  data = dd
  proteins <- unique(as.character(data$protein))
  lipids <- unique(as.character(data$lipid))
  mat <- matrix(data=NA,nrow=length(proteins),ncol=length(lipids),dimnames=list(proteins,lipids))
  for(pp in proteins){
    sub.p <- subset(data,data$protein == pp)
    for(ll in lipids){
      sub.pl <- subset(sub.p,sub.p$lipid == ll)
      if(nrow(sub.pl) > 0){
        v <- c()
        if(values == "sr"){
          v <- sub.pl$signal.ratio
          if(method == "mean"){
            mat[pp,][ll] <- mean(v) 
          }
        }else{
          v <- sub.pl$man.annotation
          if(length(which(v == -1)) > length(which(v == 1))){
            mat[pp,][ll] <- -1
          }else if(length(which(v == 1)) > length(which(v == -1))){
            mat[pp,][ll] <- 1
          }else{
            mat[pp,][ll] <- 0
          }
        } 
      }
    }
  }
  return(mat)
}


loadNewData <- function(including.controls = FALSE,just.three = TRUE){
  if(including.controls == TRUE){
    if(just.three == TRUE){
      return(read.table("/Users/deghou/remote/projects/liposomes/newdata/data.on.the.04.02.2013/all.without.cc.without.bq.wihtout.d.q.withoutNA.withoutSyntaxine.without.45678"))
    }else{
      return(read.table("/Users/deghou/remote/projects/liposomes/newdata/data.on.the.04.02.2013/all.without.cc.without.bq.wihtout.d.q.withoutNA.withoutSyntaxine"))      
    }
  }else{
    if(just.three){
  #    return(read.table("/Users/deghou/remote/projects/liposomes/newdata/data.on.the.04.02.2013/CACACACA",header=TRUE))
#      return(read.table("/Users/deghou/remote/projects/liposomes/newdata/data.on.the.04.02.2013/all.without.cc.without.bq.wihtout.d.q.withoutNA.withoutSyntaxine.withoutControls.without.45678.withoutPrecipitate",header=TRUE))
      return(read.table("/Users/deghou/remote/projects/liposomes/newdata/xx",header=TRUE,sep="\t"))
    }else{
      return(read.table("/Users/deghou/remote/projects/liposomes/newdata/data.on.the.04.02.2013/all.without.cc.without.bq.wihtout.d.q.withoutNA.withoutSyntaxine.withoutControls"))    
    }
  }
}

loadGoodRuns <- function(){
  return(as.character(read.table("/Users/deghou/remote/projects/liposomes/to.consider")$V1))
}


loadManualBinding <- function(run){
  annotations.dir <- "/Users/deghou/remote/projects/liposomes/annotations/manual-binding"
  if(file_test("-f",paste(annotations.dir,paste("manual-binding-",run,".txt",sep=""),sep="/"))){
    return(read.table(paste(annotations.dir,paste("manual-binding-",run,".txt",sep=""),sep="/")))
  }else{
    return("NULL")
  }
  return(tolower(c$V1))
}

loadProteinAdjacentDomain <- function(){
  return(as.matrix(read.delim("../../deghou/remote_copy/projects/liposomes/data/runs/final_files/central_file_good/files/lower_case/signal_ratio/domains/domain.composition")))
}


loadInfoAllAvailabaleReplicatesRatioAndDecision <- function(){

return(read.table("/Users/deghou/remote_copy/projects/liposomes/data/runs/final_files/central_file_good/files/lower_case/signal_ratio/lipid_interaction/data",sep="\t"))

}

loadCrossContaminated <- function(){
  file <-"/Users/deghou/remote/projects/liposomes/cleanup/crosscontamination_Sat_Dec_29_14:59:44_2012"
  return(read.table(file))
}

loadBadQuality <- function(){
  file <-"/Users/deghou/remote/projects/liposomes/cleanup/badquality"
  return(read.table(file))
}

loadDubiousQuality <- function(){
  file <-"/Users/deghou/remote/projects/liposomes/cleanup/dubiousquality"
  return(read.table(file))
}

loadManualColoringLocaton <- function(){
        l <- read.table("/Users/deghou/remote_copy/projects/liposomes/data/proteins_restored/level1/col",sep="\t")
        u <- as.character(as.vector(l$V2))
        v <- as.character(as.vector(l$V1))
        ll <- list()
        for(k in 1:length(u)){
                print(v[k])
                ll[[tolower(u[k])]] <- v[k]

        }
#       print(ll)
        return(ll)
}

loadManualColoring <- function(){
	l <- read.table("/Users/deghou/remote/projects/liposomes/data/lipids/ind.col",sep="\t")
	u <- as.character(as.vector(l$V2))
	v <- as.character(as.vector(l$V3))
	ll <- list()
	for(k in 1:length(u)){
#		print(v[k])
		ll[[tolower(u[k])]] <- v[k]
	
	}
#	print(ll)
	return(ll)
}

loadHomogeneousLipids <- function(){
	return(read.table("/Users/deghou/remote_copy/projects/liposomes/data/lipids/lipids.individualities",sep=" "))
}


loadHeterogeneousLipids <- function(){
	return(read.table("/Users/deghou/remote_copy/projects/liposomes/data/lipids/lipids.mixture",sep=" "))
}

loadControlsLipids <- function(){
        return(read.table("/Users/deghou/remote_copy/projects/liposomes/data/lipids/controls",sep="\t"))
}


loadDataScreening <- function(){
	return(read.table("/Users/deghou/remote_copy/projects/liposomes/data/runs/final_files/central_file_good/files/lower_case/signal_ratio/prot_lip_round1_round2_round3_TEST_WITH_DEC_AVG_NEW_ROUND_AVG_LESS_PROTEINS",sep="\t"))
}

loadProteinInformation <- function(){
	prot_info_file <- "/Users/deghou/remote/projects/liposomes/data/proteins_restored/big_protein_information_file.3"
	d <- read.table(prot_info_file,sep="\t",header=TRUE)
	return(d);
}

loadLipidInformation <- function(){
	lip_info_file <- "/Users/deghou/remote_copy/projects/liposomes/data/lipids/big_lipid_info.txt"
	d <- read.table(lip_info_file,sep="\t",header=TRUE)
	return(d);
}



createMappingLipName <- function(){
  lip.info <- loadLipidInformation()
  l <- list()
  for(k in 1:length(lip.info$lipid)){
  	id <- tolower(as.character(lip.info$lipid[k])) 	
  	name <- tolower(as.character(lip.info$name[k]))
  #	print(paste(id,name))
  	l[[id]] <- name
  }
  return(l)
}


replaceLipidByNames <- function(m){
	l <- createMappingLipName()
	n <- c()
	for(k in 1:length(rownames(m))){
	 id <- as.character(rownames(m)[k])
	 name <- as.character(l[id])
#   print(name)
#	print(strsplit(name,"_")[[1]][1])
# n[k] <- strsplit(name,"_")[[1]][1]
   n[k] <- name
	}
	rownames(m) <- n
	return(m)
}




setMapRunProtein <- function(sens="runprot") {
  file <- "/Users/deghou/remote_copy/projects/liposomes/data/runs/final_files/central_file_good/files/lower_case/signal_ratio/cleaning/mapping.run.prot"
  f <- read.table(file)
  map <- list()
  colkey = f$V1
  colval = f$V2
  if(!sens == "runprot"){
    colkey = f$V2
    colval = f$V1
  }
  run <- ""
  for(k in 1:length(colkey)){
    r <- as.character(colkey[k])
    if(run != r){
       map[[r]] <- as.character(colval[k]) 
    }
    run = k;
  }
  return(map)
}

loadDataAsMatrix <- function(values="sr",file=NULL){
  if(is.null(file)){
    d <- read.table("/Users/deghou/remote_copy/projects/liposomes/data/runs/final_files/central_file_good/files/lower_case/signal_ratio/prot_lip_round1_round2_round3_TEST_WITH_DEC_AVG_NEW_ROUND_AVG_LESS_PROTEINS",sep="\t")
    if(values == "sr"){
      return(as.matrix(xtabs(d$V3 ~ d$V1 + d$V2)))
    }else if(values == "dec"){
      return(as.matrix(xtabs(d$V4 ~ d$V1 + d$V2)))
    } 
  }else{
    d <- read.table(file,sep="\t")
    return(as.matrix(xtabs(d$V3 ~ d$V1 + d$V2)))
    
  }

}





