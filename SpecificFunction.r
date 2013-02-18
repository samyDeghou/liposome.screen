#!/usr/bin/r


source("/Users/deghou/PHdomain/CleanData.R")
source("/Users/deghou/PHdomain/RobustnessScreening.R")





searchPositionAcrossAlignmentWithNoIntersectionBetweenGroups <- function(){
  setwd("/Users/deghou/remote/projects/liposomes/newdata/alignments/")
  
  methods = list.dirs()
  methods = methods[-1]
  for(method in methods){
    setwd(method)
    lipids = list.files()
    results = ""
    print(method)
    for(lipid in lipids){
      #print(lipid)
      alignment = seqinr::read.alignment(lipid,format="fasta")
      res = performSearch(alignment)
      if(!is.null(res)){
#        results = paste(results,res)
        res = paste(sapply(res,FUN=function(x){paste(lipid,x,sep="\t")},USE.NAMES=FALSE),collapse='\n')
 #       print(results)
        write.table(x=res,file="results",append=TRUE)
      }
    }
    #write.table(x=results,file="results")
    setwd("../")
  }
}


setGroupsDecisionsOnNultipleAlignments <- function(){
  methods = c("multialin","muscle","tcoffee","clustal","mafft")       
  hom.lip = as.character(loadHomogeneousLipids()$V2)
  for(method in methods){
    print(method)
    for(lip in hom.lip){
#      print(paste(method,hom.lip))
      setAlignmentGroupForLipid(lip,method)
    }
  }
}

#@param an alignment in fasta format
#@return position in the multiple alignment displaying no intersecntion between group 1 and 2 (binders and non binders)
performSearch <- function(alignment){
  positions = c()
  for(k in 1:nchar(alignment$seq[[1]])){
    seq1 = sapply(alignment$seq[grep("GROUP:1",alignment$nam)],FUN=function(x){return(substr(x,k,k))})
    seq2 = sapply(alignment$seq[grep("GROUP:2",alignment$nam)],FUN=function(x){return(substr(x,k,k))})
    seq3 = sapply(alignment$seq[grep("GROUP:3",alignment$nam)],FUN=function(x){return(substr(x,k,k))})
    #      print(append(k,intersect(seq1,seq2)))
    #     print(length(intersect(seq1,seq2)))
    if(length(intersect(seq1,seq2)) == 0){
      positions = append(positions,k)
      #        print(paste(file,k))
      #        print(seq2)
      #       print(seq3)
    }
  }
  if(length(positions) > 0){
    return(positions)
  }else{
    return(NULL)
  }
}

detectManually <- function(){
#  setwd("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/grouping.auto.generated/")
  setwd("/Users/deghou/phd/ressources/data/liposomes/data/motif.search/")
  files = list.files(pattern='global.fa')
  for(file in files){
    print(file)
    alignment = read.alignment(file=file,format="fasta")
    res = performSearch(alignment)
    if(!is.null(res)){
      print(res)
    }
  }
}


whoHaspattern <- function(pattern=""){

  pattern = '[m|l|v|f|y|i].k.[g|a|s|p].*[k|r]*.*[r|k].[r].[f|l]'
  prot <- c()
  have = c()
  
  for(k in 1:84){
     name <- as.character(protein.informations$reserved_name[k])
     seq <- as.character(protein.informations$sequence[k])
#     name = "pdk1"
#     seq = "gsnieqyihdldsnsfeldlqfsedekrlllekqaggnpwhqfvennlilkmgpvdkrkglfarrrqllltegphlyyvdpvnkvlkgeipwsqelrpeaknfktffvhtpnrtyylmdpsgnahkwcrkiqevwrqryqshpdaavq"
     res = sub(pattern,x=seq,perl=TRUE,replacement="ZZZZZ")
     res1 = strsplit(res,"")[[1]]
     int <- length(grep(c("Z","Z","Z","Z","Z"),res1))
     if(int == 5){
#       print("FOUND 1 !!")
       have = append(have,name)
#       print(name)
       prot <- append(prot,name)
#       print(seq)
#       print(append(res)
     }
  }
  return(have)  
}



quickLook <- function(prot,lip){
  return(subset(data,subset=data$protein %in% prot & data$lipid.composition %in% lip))
}


predictSecondaryStructures <- function(protein){

  seq <- subset(protein.informations$sequence,subset = protein.informations$reserved_name == protein)
  name <- subset(protein.informations$sequence,subset = protein.informations$fasta_name == protein)
  write.fasta(sequences=seq,names=name,file.out="/Users/deghou/remote/utils/tools/psipred/tmp.seq.to.be.deleted.fa")
  setwd("/Users/deghou/remote/utils/tools/psipred/")
  launch.psi.pred <- "./runpsipred_single tmp.seq.to.be.deleted.fa"
  delete.seq.file <- "rm tmp.seq.to.be.deleted.fa"
  file <- paste("/Users/deghou/remote/projects/liposomes/data/proteins_restored/sequence/predicted.secondary.structure/",protein,".secondary.structure.psi.pred",sep="")
  mv.to.dir <- paste("mv tmp.seq.to.be.deleted.horiz",file)
  rm.examples <- "rm tmp.seq.to.be.deleted.*"
  system(launch.psi.pred)
  system(delete.seq.file)
  system(mv.to.dir)
  system(rm.examples)
  setwd("/Users/deghou")
}

#@data the data set to work with
#@protein.to.work.with sometimes we want to exclude some proteins .. worth to have it then
setAlignmentGroupForLipid <- function(lipid,method){
  id <- getLipidId(lipid)
  protein.informations <- loadProteinInformation()
  #binders and non binders results
  res <- getLipidBindersNonBindersAndUnclear(id,data,protein.to.work.with)
  res.binders.fasta <- as.character(subset(protein.informations$fasta_name,subset = protein.informations$reserved_name %in% as.character(res$binders)))
  res.unclear.fasta <- as.character(subset(protein.informations$fasta_name,subset = protein.informations$reserved_name %in% as.character(res$unclear)))
  #make sure we can work with these proteins
  res.binders.fasta <- res.binders.fasta[which(res.binders.fasta %in% protein.to.work.with)]
  res.unclear.fasta <- res.unclear.fasta[which(res.unclear.fasta %in% protein.to.work.with)]
  #the global alignment
  alignment <- loadMultipleAlignment(method=method)
  all <- alignment$nam
  
  proteins <- tolower(alignment$nam)
  ar <- proteins
  if(length(res.binders.fasta) != 0){
    ar <- addGroupsToMultipleAlignment(res.binders.fasta,"GROUP:1|",ar) 
  }
  if(length(res.unclear.fasta) != 0){
    ar <- addGroupsToMultipleAlignment(res.unclear.fasta,"GROUP:3|",ar) 
  }
  non.binders <- subset(ar,subset=!substr(ar,0,5) == "GROUP")
  ar <- addGroupsToMultipleAlignment(non.binders,"GROUP:2|",ar)
  alignment$nam <- ar
  print(paste("writing alignment",method," about lipid",lipid))
  
  for(k in 1:alignment$nb){
#    if(lipid == "DHS_10" && as.character(alignment$nam[[k]]) == "GROUP:3|avo1|1056-117"){
#      real.indice <- which(tolower(all) == substr(alignment$nam[[k]],9,nchar(alignment$nam[[k]])))
#    }
 #   print(alignment$nam[[k]])
  #  print(substr(alignment$nam[[k]],9,nchar(alignment$nam[[k]])))
    real.indice <- which(tolower(all) == substr(alignment$nam[[k]],9,nchar(alignment$nam[[k]])))
    dir = "/Users/deghou/remote/projects/liposomes/newdata/alignments/"
    
    write(x=paste(paste(">",alignment$nam[[k]],sep=""),alignment$seq[[real.indice]],sep="\n"),file=paste(dir,method,"/",lipid,sep=""),append=TRUE)
#    write(x=paste(paste(">",alignment$nam[[k]],sep=""),alignment$seq[[real.indice]],sep="\n"),file=paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/grouping.auto.generated/",lipid,sep=""),append=TRUE)
    
  }
#  seqinr::write.fasta(sequences=alignment$seq,names=alignment$nam,file.out=paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/grouping.auto.generated.2/",lipid,sep=""))
#  return(alignment)
}

addGroupsToMultipleAlignment <- function(to.be.tagged,tag,all.proteins){
  f <- function(x){
    return(paste(tag,x,sep=""))
  }
  untagged <- setdiff(all.proteins,to.be.tagged)
  tagged <- sapply(to.be.tagged,f,USE.NAMES=FALSE)
  return(append(tagged,untagged))
}


preapreBreaks <- function(mat.vec,vec1,vec2){
  #sort the vector
  mat.vec <- mat.vec[order(mat.vec)]
  min <- mat.vec[1]
  max <- mat.vec[length(mat.vec)]
  for(k in 1:length(breaks)){
    actual.break <- breaks[k]
    below.break <- NULL;
    if(k == 1){
      below.break <- 0
    }else{
      below.break <- breaks[k-1]
    }
    #sub vector in between the two breaks
    vec.in.between <- mat.vec[which(mat.vec < actual.break & mat.vec > below.break)]
    mat.vec.u <- l
  }  
}






plotOccurenceAAinGroup <- function(position,lipid,group,alignment=NULL){
  fasta = NULL
  if(is.null(alignment)){
    fasta <- "null"    
  }else{
    fasta <- alignment
  }
  df <- getPositionSummary(position,lipid,group,fasta)
  plot.new()
  x = 0.1
  y = 0
  inc = 1/nrow(df)
  for(k in 1:nrow(df)){
    value <- as.character(df$occurence[k]/sum(df$occurence))
    value <- as.numeric(value) * 100
    value <- paste(substr(value,0,4),"%")
    text(x,y,label=as.character(df$s[k]),cex=as.numeric(as.character(df$occurence[k]))/3,col="blue")
    text(0.5,y,label=value,cex=1)
    y = y + inc
  }
}

getPositionSummary <- function(position,lipid,group,alignment=NULL){
  fasta = NULL
  if(is.null(alignment)){
    fasta <- loadMultipleAlignment()    
  }else{
    fasta <- alignment
  }
  seq <- getPosAA(fasta,position,group)
  seq <- seq[-which(seq == "")]
  df <- characterOccurence(seq)
  return(df)
}

getPosAA <- function(alignment,pos,group){
  sequence <- ""
  for(s in 1:alignment$nb){
    name <- alignment$nam[s]
    if(substring(alignment$nam[s],0,7) == paste("GROUP:",group,sep="")){
      sequence = append(sequence,substr(alignment$seq[[s]],pos,pos))
    }
  }
  return(sequence)
}



generateBindersNonBindersPerLipid <- function(){
  prot.to.work.with <- as.character(read.table("/Users/deghou/remote/projects/liposomes/newdata/motif.search/to.work.with")$V1)
  prot.to.wirk.with.names <- as.character(subset(protein.informations$reserved_name,subset = protein.informations$fasta_name %in% prot.to.work.with))
  lip.hom <- loadHomogeneousLipids()
  data <- loadNewData()
  for(k in 1:nrow(lip.hom)){
    lip.id <- tolower(as.character(lip.hom$V1[k]))
    lip.name <- lip.hom$V2[k]
 #   print(lip.id)
    df <- getLipidBindersNonBindersAndUnclear(lip.id,data,prot.to.work.with)  
   #   print(paste(lip.id,length(df$binders),length(df$nonbinders)))
      binders <- df$binders
      unclear <- df$unclear
      nonbinders <- subset(prot.to.wirk.with.names,subset = !(prot.to.wirk.with.names %in% df$binders) & !(prot.to.wirk.with.names %in% df$unclear))
      for(protein in binders){
        sequence <- as.character(subset(protein.informations$sequence,subset = protein.informations$reserved_name == protein))
        fasta.name <- as.character(subset(protein.informations$fasta_name,subset = protein.informations$reserved_name == protein))
        write(x=paste(">GROUP:1","|",substr(fasta.name,start=2,stop=nchar(fasta.name)),"\n",sequence,sep=""),file=paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/",lip.name,sep=""),append=TRUE)
      }
      for(protein in nonbinders){
        sequence <- as.character(subset(protein.informations$sequence,subset = protein.informations$reserved_name == protein))
        fasta.name <- as.character(subset(protein.informations$fasta_name,subset = protein.informations$reserved_name == protein))
        write(x=paste(">GROUP:2","|",substr(fasta.name,start=2,stop=nchar(fasta.name)),"\n",sequence,sep=""),file=paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/",lip.name,sep=""),append=TRUE)
      }
      for(protein in unclear){
        sequence <- as.character(subset(protein.informations$sequence,subset = protein.informations$reserved_name == protein))
        fasta.name <- as.character(subset(protein.informations$fasta_name,subset = protein.informations$reserved_name == protein))
        write(x=paste(">GROUP:3","|",substr(fasta.name,start=2,stop=nchar(fasta.name)),"\n",sequence,sep=""),file=paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/",lip.name,sep=""),append=TRUE)
      }    
  }
}
idd <- function(lipid){
  return(which(ll == lipid))
}
#for this function load :
# loadNewData()
# proteinInformation()
# proteinToWorkWith("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/to.work.with")
getLipidBindersNonBindersAndUnclear <- function(lip.id,data,protein.to.work.with){
  data.lip <- subset(data,subset= data$lipid == lip.id)
  proteins <- unique(as.character(data.lip$protein))
  protein.with.at.least.one.positive.annotation <- unique(as.character(subset(data.lip$protein,data.lip$man.annotation == 1)))
  binders <- c()
  nonbinders <- c()
  unclear <- c()
  for(protein in protein.with.at.least.one.positive.annotation){
  # print(protein)
    fasta.name <- as.character(subset(protein.information$fasta_name,subset = protein.information$reserved_name == protein))
    #print(fasta.name)
    if(protein != "hsv2" && fasta.name %in% protein.to.work.with){
      annotations <- as.character(subset(data.lip$man.annotation, subset = data.lip$protein == protein))
      ratios <- as.numeric(as.character(subset(data.lip$signal.ratio, subset = data.lip$protein == protein)))
      dec <- makeDecision(annotations,ratios)
#      print(dec)
      if(dec == "unclear"){
        unclear = append(unclear,protein)
      }else if(dec == "bind"){
        binders = append(binders, protein)
      }else if(dec == "does not bind"){
        nonbinders = append(nonbinders, protein)
      }
    }
  }
  return(list("binders"=binders,"non binders"=nonbinders,"unclear"=unclear))
}


plotBestMotifCandidates <- function(d){
#  ggplot(d, aes(x=SH_score, y = mR_score, colour = W_score,size = Z_score)) +
#    geom_point() + geom_vline() + geom_hline() + opts(title = "Identification of funciton residues")
#  + scale_color_gradient(low="white",hihg="black")
  print("kpk")
  qplot(
    xlab = "SH score",
    ylab = "mR score",
    SH_score,mR_score,data=d,colour=W_score,size=Z_score) + geom_hline(aes(yintercept=0.35)) + 
    scale_color_gradient(low="yellow",high="blue")
    labs(title = "Identification of functional residues")
}

getBestMotifCandidates <- function() {
  lipids <- "cardiolipin_10\nceramide_10\nceramide-1p_10\ndag_5\ndhs_10\ndhs-1p_7\ndhs-1p_10\ndi-hydro-ceramide_10\nceramide_10\ndopa_10\ndopg_10\ndopi_10\ndopi3p_10\ndopi5p_10\ndopi34p2_10\ndopi35p2_10\ndopi45p2_7\ndopi45p2_10\ndopi345p3_10\ndops_10\ndppi4p_10\nphs_10\nphs-1p_10\nphyto-ceramide_10\npopc_10\npope_10\ns-1p_10\nsphingosine_10"
  lipids <- toupper(lipids)
  lipids <- strsplit(lipids,"\n")[[1]]
  d <- data.frame()
  couleurs <- rainbow(length(lipids))
  ind <- 0
  for(lipid in lipids){
#    print(lipid)
    ind = ind + 1
    if(file.exists(paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/output.results/",lipid,"_SH",sep=""))){
      results <- parseResultsForLipid(lipid)
      results.bests <- selectBestCandidates(results)
      results.bests$col <- rep(couleurs[ind],dim(results.bests)[1])
      results.bests$lip <- rep(lipid,dim(results.bests)[1])
      #    print(results.bests)
      d <- rbind(d,results.bests) 
    }else{
      print(paste(file," does not exist"))
    }
  }
  print(names(d))
  d$motif.1 <- as.character(d$motif.1)
  d$motif.2 <- as.character(d$motif.2)
  d$motif.3 <- as.character(d$motif.3)
  dd <- data.frame()
  X <- c(-3.0,0.98)
  Y <- c(1.53,1.26)
  motif.1 <- c("AMS","KR")
  motif.2 <- c("QRP","KR")
  motif.3 <- c("KVE","RK")
  SH_score <- c(0.99,0.01)
  mR_score <- c(0.5,-.19)
  col <- c("caca","caca")
  lip <- c("positive control","negative control")
  dd <- data.frame(X,Y,SH_score,mR_score,motif.1,motif.2,motif.3,col,lip)
  d <- rbind(d,dd)
  names(d) <- c("Z_score","W_score","SH_score","mR_score","motif.1","motif.2","motif.3","col","lipids")  
  return(d)
}

selectBestCandidates <- function(d){
  return(head(d[with(d, order(X,-Y)), ],n=10))
}


parseResultsForLipid <- function(lipid){
  file.sh <- paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/output.results/",lipid,"_SH",sep="")
  file.mr <- paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/output.results/",lipid,"_MR",sep="")
  con <- file(file.sh, open = "r")
  save = FALSE
  d <- data.frame()
  X <- c()
  Y <- c()
  motif.1 <- c(); motif.2 <- c(); motif.3 <- c()
  names <- c()
  SH_score <- c()
  mR_score <- c()
  a <- 0;
  more = FALSE
  while(length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if(nchar(line) > 17 && substr(line,0,18) == "Alignment Position"){
      save = TRUE
      print(strsplit(line,"\t")[[1]])
      if(length(strsplit(line,"\t")[[1]]) == 10){
        more = TRUE
      }
      line <- readLines(con, n = 1, warn = FALSE)
    }
    if(save){
      a = a + 1
      t <- strsplit(line,"\t")
      t <- t[[1]]
      X[length(X) + 1] <- t[3]
      motif.1[length(motif.1) + 1] <- t[6]
      motif.2[length(motif.2) + 1] <- t[8]
      if(more){
        motif.3[length(motif.3) + 1] <- t[10]        
      }else{
        motif.3[length(motif.3) + 1] <- "NO THIRD GROUP"
      }
      # CAREFUL !! 1 - SH_score
      SH_score[length(SH_score) +1] <- 1 - as.numeric(t[2])
    }
  }
  save = FALSE
  a = 0
  con <- file(file.mr, open = "r")
  while(length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if(nchar(line) > 17 && substr(line,0,18) == "Alignment Position"){
      save = TRUE
      line <- readLines(con, n = 1, warn = FALSE)
    }
    if(save){
      a = a + 1
      t <- strsplit(line,"\t")
      t <- t[[1]]
      Y[length(Y) + 1] <- t[5]
      mR_score[length(mR_score) + 1] <- t[2]
    }
  }
#  print(paste(length(X),length(Y)))
  # print(a)
  X <- as.numeric(as.character(X))
  Y <- as.numeric(as.character(Y))
  mR_score <- as.numeric(as.character(mR_score))
  return(data.frame(X,Y,SH_score,mR_score,motif.1,motif.2,motif.3))
}


downloadHarmonyResults <- function(){
  data <- read.table("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/output.results/tmp",sep="\t")
  urls <- as.character(data$V2)
  lipids <- as.character(data$V1)
  for(k in 1:length(lipids)){
    lipid <- lipids[k]
    url <- urls[k]
    if(url != "same sequence"){
      print(paste(lipid,url))
      sh.results <- paste(url,"SH.out",sep="")
      mr.results <- paste(url,"MR.out",sep="")
      write(x=RCurl::getURL(sh.results),file=paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/output.results/",lipid,"_SH",sep=""))
      write(x=RCurl::getURL(mr.results),file=paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/subset/output.results/",lipid,"_MR",sep="")) 
    }
  }
}



plotZWscoresFromAlignments <- function(){
  X <- c(-5,-4,-2,0,1)
  Y <- c(-1.3,-3.4,-0.89,3.4,1.9)
  col <- c("red","blue","purple","orange","grey")
  label <- c("D_R_T","B_H_H","WCECW","CCW","grey")
  d <- data.frame(X,Y,col,label)
  
}


doThing <- function(){
  col <- unique(as.vector(which(coop.mat == 1,arr.ind=TRUE)[,"col"]))
  col1 <- col[1:13]
  col2 <- col[14:25]
  col3 <- col[26:37]
  col4 <- col[38:49]
  col5 <- col[50:54]
  lm <- createMappingLipName()
  par(mfrow=c(3,4))
  for(c in col5){
    lipid.id <- colnames(coop.mat)[c]
    lipid.name <- as.character(lm[(which(lm == lipid.id))])
    proteins.pos <- names(which(coop.mat[,c] == 1))
    proteins.neg <- names(which(coop.mat[,c] == -1))
    print(lipid.name)
    if(length(proteins.neg) > 0){
      print("STOOOOOPP")
    }
    checkCooperativityOfLipids(lipid.name,proteins.pos,proteins.neg)
  }
}

checkCooperativityOfLipids <- function(mixture.c,proteins.to.label = c(),proteins.to.label.neg = c()){
  mixture <- strsplit(mixture.c,":")
  mixture <- mixture[[1]]
  concentrations <- mixture[which(nchar(mixture) == 2)]
  if(length(concentrations) == 0){
    concentrations <- mixture[which(nchar(mixture) == 1)]    
  }
  lipid.names <- mixture[which(nchar(mixture) > 2)]
  first <- paste(mixture[1],mixture[length(mixture)],sep="_")
  if(first %in% colnames(mat.27) & mixture[2] %in% colnames(mat.27)){
    c1 <- mat.27[,first]
    c2 <- mat.27[,mixture[2]]
#    c1 <- mat.27[,mixture[2]]
#    c2 <- mat.27[,first]
    
    c12 <- mat.110[,mixture.c]
    c3 <- c12 - c1
    reprep <- rep(1,times=length(c1))
    nnn <- names(c1)
    d <- data.frame(x=as.vector(c1),y=as.vector(c12),names=as.character(as.vector(nnn)),z=as.character(as.vector(reprep)))
    sss <- rep("black",times=length(c1))
    test <- cor.test(c1,c12)
    fit <- lm(c12~c1)
    for(kk in 1:length(d$names)){
      prot <- d$names[kk]
      if(length(proteins.to.label.neg) == 0){
        if(!(prot %in% proteins.to.label)){
          nnn[kk] <- ""
        }else{
          reprep[kk] <- 2
          sss[kk] <- "red"
        } 
      }else{
        if(prot %in% proteins.to.label){
          reprep[kk] <- 2
          sss[kk] <- "red"
        }else if(prot %in% proteins.to.label.neg){
          reprep[kk] <- 3
          sss[kk] <- "green"          
        }else{
          nnn[kk] <- ""
        }
      }
    }
    d <- data.frame(x=as.vector(c1),y=as.vector(c12),names=as.character(as.vector(nnn)),z=as.character(as.vector(reprep)),size=sss)
    plot(c1,c12,main=paste(
      mixture.c,paste("R^2= ",
      substr((test$estimate * test$estimat),0,4)),sep="\n"),
      xlab=mixture[2],ylab=mixture.c,cex=0.1,cex.main=1)    
  #  plot(c1,c12,main=paste(mixture.c,paste("R^2= ",substr((test$estimate * test$estimat),0,4)),sep="\n"),xlab=mixture[1],ylab=mixture.c,cex=0.1,cex.main=1)    
    abline(fit,col="orange")
    if(length(proteins.to.label) > 1){
      legend("topright",replaceCharacterInString(toString(proteins.to.label),",","\n"),cex=0.6)
    }else{
      legend("bottomright",proteins.to.label)
    }
#    legend("topleft",,cex=0.6)
    dd <- subset(d,d$z==2)
    ddd <- subset(d,d$z==3)
    lll <- nrow(dd)
    llll <- nrow(ddd)
    for(item in 1:lll){
    #  print(paste(dd$x[item],dd$y[item],dcol="red",labs=dd$names[item]))
      points(dd$x[item],dd$y[item],col="darkgreen")
      calibrate::textxy(X=dd$x[item],Y=dd$y[item],dcol="darkgreen",labs=dd$names[item],cx=0.7)
    }
    if(length(proteins.to.label.neg) > 0){
      for(item in 1:llll){
        #  print(paste(dd$x[item],dd$y[item],dcol="red",labs=dd$names[item]))
        points(ddd$x[item],ddd$y[item],col="purple")
        calibrate::textxy(X=ddd$x[item],Y=ddd$y[item],dcol="purple",labs=ddd$names[item],cx=0.7)
      }
    }
  }else{
    print(paste("One of the two members of " , mixture.c, " was not present in the individual lipids"))
  }
  return(d)
}

plotProtCoopProfile <- function(proteins,color=NULL){
  coop.mat <- doCooperativness()
  par(las=2)
  par(cex.axis=0.7)
#  print("cooperative matrix done")
 # lip.het <- loadHeterogeneousLipids()
  if(is.null(color)){
    color = sample(colors(),size=length(proteins))
  }
  if(length(proteins) == 1){
    pp <- proteins[1]
    profile <- coop.mat[pp,][which(!is.na(coop.mat[pp,]))]
    max <- max(profile)
    min <- min(profile)
    plot(profile,main=pp,xlab="",ylab="Observed - Expected",ylim=c(min,max))
    lines(profile,col=color[1])
#    axis(1,83,colnames(coop.mat))
  }else{
    mmax <- 0
    mmin <- 67
    for(c in proteins){
      pp <- c
      profile <- coop.mat[pp,][which(!is.na(coop.mat[pp,]))]
      if(max(profile) > mmax){
        mmax <- max(profile)
      }
      if(min(profile) < mmin){
        mmin <- min(profile)
      }
    }
    print(paste(mmin,mmax))
    for(k in 1:length(proteins)){
      pp <- proteins[k]
    #  print(pp)
      if(k == 1){
        profile <- coop.mat[pp,][which(!is.na(coop.mat[pp,]))]
        plot(profile,main="Lipid cooperativiy",xlab="",ylab="Observed(GFP/Cy5) - Expected(GFP/Cy5)",ylim=c(mmin,mmax))
        lines(profile,col=color[1])
        axis(side=1,at=1:length(colnames(coop.mat)),colnames(coop.mat))
      }
      profile <- coop.mat[pp,][which(!is.na(coop.mat[pp,]))]      
      lines(profile,col=color[k])
    }
    mat <- coop.mat
    print(min)
    print(max)
    plotrix::color.legend(85,mmin - 2,86,mmax + 0.1,legend=proteins,rect.col=color,gradient="y",cex=0.6,align="rb")
  }
}


makeDecisionHeatmap <- function(){
#   data <- loadNewData()
#   mat.110 <- doMatrixOutOfScreening()
   lip.het <- loadHeterogeneousLipids()
   mat.27 <- reduceMatrix(mat.110,lip.het$V1,col=TRUE)
   mat.27 <- t(replaceLipidByNames(m=t(mat.27)))
   mat.110 <- t(replaceLipidByNames(m=t(mat.110)))
   proteins <- unique(as.character(data$protein))
   hets <- tolower(as.character(lip.het$V2))
   homs <-tolower(as.character(lip.hom$V2))
   coop.mat <- matrix(data=-1,nrow=length(proteins),ncol=length(hets),dimnames=list(proteins,hets))
  # print(length(which(is.na(coop.mat))))
   for(pp in proteins){
     for(ll in hets){
#      print(paste(pp,ll))
      sr <- subset(data$signal.ratio,data$protein == pp & data$lipid.composition == ll)
      man <- subset(data$man.annotation,data$protein == pp & data$lipid.composition == ll)
    #  print(pp)
    #  print(sr)
    #  print(man)
      if(pp == "boi2-ph" & ll == "dopi34p2:di-hydro-ceramide_10:10" ){
        print(sr)
        print(man)
        print(dec)
      }
      if(length(sr) > 0 & length(man) > 0){
        dec = makeDecision(man,sr)
        if(pp == "boi2-ph" & ll == "dopi34p2:di-hydro-ceramide_10:10" ){
          print(sr)
          print(man)
          print(dec)
        }
        if(dec == "bind"){
          coop.mat[pp,][ll] <- 1
        }else if(dec == "unclear"){
          coop.mat[pp,][ll] <- 0
        } 
      }else{
        coop.mat[pp,][ll] <- 0        
      }
    }
  }
   return(coop.mat)
}



doCooperativness <- function(){
#  data <- loadNewData()
#  mat.110 <- doMatrixOutOfScreening()
#  lip.het <- loadHeterogeneousLipids()
#  mat.27 <- reduceMatrix(mat.110,lip.het$V1,col=TRUE)
#  mat.27 <- t(replaceLipidByNames(m=t(mat.27)))
#  mat.110 <- t(replaceLipidByNames(m=t(mat.110)))
#  proteins <- unique(as.character(data$protein))
#  hets <- tolower(as.character(lip.het$V2))
#  coop.mat <- matrix(data=NA,nrow=length(proteins),ncol=length(hets),dimnames=list(proteins,hets))
 # print(length(which(is.na(coop.mat))))
  for(protein in proteins){
    print(protein)
    for(het in hets){
#      value <- getResultsFromLipidComposition(protein,het)
#      value <- getSpecials(protein,het)
      value <- getSpecialsWithLessStringence(protein,het)
 #     print(paste(protein,het,value))
      coop.mat[protein,][het] <- value
    }
  }
#  print(length(which(is.na(coop.mat))))
  coop.mat[which(is.na(coop.mat))] = 0
  heatmap.2(coop.mat,Rowv=NULL,Colv=NULL,trace="none",na.color="blue",col="redgreen",margins=c(10,4),key=FALSE)
  return(coop.mat)
}


getSpecialsWithLessStringence <- function(protein,test){
  test.c <- strsplit(test,":")
  concentrations <- test.c[[1]][which(nchar(test.c[[1]]) == 2)]
  if(length(concentrations) == 0){
    concentrations <- test.c[[1]][which(nchar(test.c[[1]]) == 1)]    
  }
  lipid.names <- test.c[[1]][which(nchar(test.c[[1]]) > 2)]
  inds <- c()
  for(k in 1:(length(concentrations) )){
    inds <- append(inds,paste(as.character(lipid.names[k]),as.character(concentrations[k]),sep="_"))
  }
  inds <- append(inds,lipid.names[length(lipid.names)])
  observed.value <- mat.110[protein,][test]
 # print(paste(protein,test,observed.value,m1.83[protein,][test]))
  if(is.na(observed.value) && is.na(m1.83[protein,][test])){
    return(NA)
  }
  if(observed.value > th && m1.83[protein,][test] == 1){
    ppp = protein
    for(ind.lip in inds){
      value <- as.numeric(mat.27[protein,][ind.lip])
      ann <- m1.27[protein,][ind.lip]
      if(!is.na(value) && !(is.na(ann))){
        mmaa = subset(data$man.annotation,subset=data$protein==ppp & data$lipid.composition==ind.lip)
        ssrr = subset(data$man.annotation,subset=data$protein==ppp & data$lipid.composition==ind.lip)
        dec = makeDecision(mmaa,ssrr)
        if(dec != "does not bind"){
          return(0)
        } 
      }else{
        return(NA)
      }
    }
    return(1)
  }else if (observed.value < th && m1.83[protein,][test] == -1){
    ppp = protein
    for(ind.lip in inds){
      value <- as.numeric(mat.27[protein,][ind.lip])
      ann <- m1.27[protein,][ind.lip]
      if(!is.na(value) && !is.na(ann)){
        mmaa = subset(data$man.annotation,subset=data$protein==ppp & data$lipid.composition==ind.lip)
        ssrr = subset(data$man.annotation,subset=data$protein==ppp & data$lipid.composition==ind.lip)
        dec = makeDecision(mmaa,ssrr)
        if(dec == "does not bind"){
          return(0)
        } 
      }else{
        return(NA)
      }
    }
    return(-1)
  }else{
    return(NA)
  }
}



getSpecials <- function(protein,test){
  test.c <- strsplit(test,":")
  concentrations <- test.c[[1]][which(nchar(test.c[[1]]) == 2)]
  if(length(concentrations) == 0){
    concentrations <- test.c[[1]][which(nchar(test.c[[1]]) == 1)]    
  }
  lipid.names <- test.c[[1]][which(nchar(test.c[[1]]) > 2)]
  inds <- c()
  for(k in 1:(length(concentrations) )){
    inds <- append(inds,paste(as.character(lipid.names[k]),as.character(concentrations[k]),sep="_"))
  }
  inds <- append(inds,lipid.names[length(lipid.names)])
  observed.value <- mat.110[protein,][test]
 # print(paste(protein,test,observed.value,m1.83[protein,][test]))
  if(is.na(observed.value) && is.na(m1.83[protein,][test])){
    return(NA)
  }
  if(observed.value > th && m1.83[protein,][test] == 1){
    for(ind.lip in inds){
      value <- as.numeric(mat.27[protein,][ind.lip])
      ann <- m1.27[protein,][ind.lip]
      if(!is.na(value) && !(is.na(ann))){
        if(value > th || ann == 1){
          return(0)
        }
      }else{
        return(NA)
      }
    }
    return(1)
  }else if (observed.value < th && m1.83[protein,][test] == -1){
    for(ind.lip in inds){
      value <- as.numeric(mat.27[protein,][ind.lip])
      ann <- m1.27[protein,][ind.lip]
      if(!is.na(value) && !(is.na(ann))){
        if(value < th || ann == -1){
          return(0)
        }
      }else{
        return(NA)
      }
    }
    return(-1)
  }else{
    return(NA)
  }
}

getResultsFromLipidComposition <- function(protein,test){
  test.c <- strsplit(test,":")
  concentrations <- test.c[[1]][which(nchar(test.c[[1]]) == 2)]
  if(length(concentrations) == 0){
    concentrations <- test.c[[1]][which(nchar(test.c[[1]]) == 1)]    
  }
  lipid.names <- test.c[[1]][which(nchar(test.c[[1]]) > 2)]
  inds <- c()
  for(k in 1:(length(concentrations) )){
    inds <- append(inds,paste(as.character(lipid.names[k]),as.character(concentrations[k]),sep="_"))
  }
  inds <- append(inds,lipid.names[length(lipid.names)])
  observed.value <- mat.110[protein,][test]
  expected.value <- NA
  for(ind.lip in inds){
    value <- as.numeric(mat.27[protein,][ind.lip])
    if(!is.na(value)){
      if(is.na(expected.value)){
        expected.value = 0
      }
      expected.value <- expected.value + value   
    }
  }
  if(is.na(expected.value)){
  #  print(paste(protein,test))
    return(NA)
  }else{
    diff <- observed.value  - expected.value
    return(as.numeric(diff))
  }  
}


interessant <- function(hc){
  l <- loadManualColoring()
  df <- extractMapToVectors(l)
  plot(as.dendrogram(hc))
  plotrix::color.legend(28,0,29,0.5,df$b,df$a,cex=0.6,gradient="y",align="rb")
  ccc <- doDoDo(hc)
}

doDoDo <- function(hc){
  par(cex=0.1)
  plot(as.dendrogram(hc))
  sorted.labels <- hc$labels[as.vector(hc$order)]
  l <- loadManualColoring()
  x.pos <- 0
  ccc <- c()
  for(lab in sorted.labels){
    x.pos <- x.pos + 1
    cc <- strsplit(lab,":")[[1]]
    vec.col <- c()
  #  print(cc)
    for(xx in cc[1:length(cc) - 1]){
     # print(xx)
      if(!grepl("_[5|7|10]",xx)){
        xx = paste(xx,"_",cc[length(cc)],sep="")
      }
      col <- l[[xx]]
      if(!xx %in% names(l)){
        col <- "white"
        print(paste(xx,"white"))
      }
      print(xx)
      print(col)
      ccc <- append(ccc,col)
#      print(paste(xx,col,sep= "    "))
      vec.col <- append(vec.col,col)
  #    print(vec.col)
  #    print(cc)
      plotrix::color.legend(x.pos - 1,.9,x.pos + 0.45,1.1,c(),rev(vec.col),gradient="y")
#      print(vec.col)
    }
  }
#  return(ccc)
}




plotProtProfile <- function(proteins,color=NULL,fileOut,lipids="all",x="ids",size.axis.labels = 0.7,dim.width=300,dim.height=300){
 # mat <- doMatrixOutOfScreening()
  th <- 0.03720974
  if(is.null(color)){
    color = sample(colors(),size=length(proteins))
  }
  print(length(color))
  if(lipids == "het"){
    lip.hom <- loadHomogeneousLipids()
    mat <- reduceMatrix(mat,lip.hom$V1,col=TRUE)
    if(x == "names"){
      mat <- t(replaceLipidByNames(t(mat)))      
    }
  }else if (lipids == "hom"){
    lip.het <- loadHeterogeneousLipids()
    mat <- reduceMatrix(mat,lip.het$V1,col=TRUE)
    if(x == "names"){
      mat <- t(replaceLipidByNames(t(mat)))      
    }
  }else if(lipids == "all" && x == "names"){
    mat <- t(replaceLipidByNames(t(mat)))    
  }
  print(dimnames(mat))
  fit <- "/Users/deghou/remote/projects/liposomes/plots/for.ivana/"
  f <- paste(fit,fileOut,sep="")
  par(cex.axis = size.axis.labels)
  par(cex.main=0.5)
  par(las=2)
  par(mar=c(7,4,1,6))
  cc <- c()
#  mat <- doMatrixOutOfScreening()
  if(length(proteins) == 1){
    plot(mat[proteins[1],],main=proteins,xlab="",ylab="GFP/Cy5")
    lines(mat[proteins[1],],col=color[1])
    axis(1,1:110,colnames(mat))
  }else{
    for(k in 1:length(proteins)){
      cc <- append(cc,mat[proteins[k],])
    }
    cc <- cc[which(complete.cases(cc))]
    mmax <- max(cc)    
    max <- 0.2 * (max(cc)/10)
    min <- 10 * (max(cc)/10)
    print(paste(mmax,max,min))
    first <- proteins[1]
    print(f)
  postscript(file=f,width=dim.width,height=dim.height,horizontal=FALSE)
    plot(mat[proteins[1],],main="lipid binding profil",xlab="",ylab="GFP/Cy5",ylim=c(0,mmax))
    lines(mat[proteins[1],],col=color[1])
    axis(1,1:(length(colnames(mat))),colnames(mat))
    abline(h=th,col="red")
    if(lipids == "hom"){
      plotrix::color.legend(length(colnames(mat)) + 5,max,length(colnames(mat)) + 6,min,proteins,color,gradient="y",cex=size.axis.labels)      
    }else if(lipids == "het"){
      plotrix::color.legend(length(colnames(mat)) + 15,max,length(colnames(mat)) + 19,min,proteins,color,gradient="y",cex=size.axis.labels)            
    }else if (lipids == "all"){
      plotrix::color.legend(length(colnames(mat)) + 23,max,length(colnames(mat)) + 24,min,proteins,color,gradient="y",cex=size.axis.labels)
    }
    for(k in 2:length(proteins)){
      lines(mat[proteins[k],],col=color[k])
    }
    dev.off()
  }
  
  
}

isConsitent <- function(protein,lipid,sr.vector,man.vector){
  ma <- (length(which(sr.vector > th))) > (length(which(sr.vector < th)))
  mu <- (length(which(sr.vector < th))) > (length(which(sr.vector > th)))
  mp <- (length(which(man.vector == 1))) > (length(which(man.vector == -1)))
  mn <- (length(which(man.vector == 1))) < (length(which(man.vector == -1)))
  if( ma && mp){
    return(1)
  }else if(ma && mn){
    return(-1)
  }else{
    write(file="/Users/deghou/remote/projects/liposomes/newdata/motif.search/UNDECIDABLE.CASES",paste(protein,lipid,sr.vector,man.vector,sep="\t"),append=TRUE)
    return(0)
  }  
}

getBinders <- function(binders=TRUE,criteria="sr",ll){
  proteins <- as.character(unique(data$protein))
  c <- c()
  for(pp in proteins){
    sub <- subset(data,subset = (data$lipid == ll) & (data$protein == pp))
    if(nrow(sub) > 0){
      if(criteria == "sr"){
        mean <- mean(sub$signal.ratio)
    #    print(paste(ll,pp))
        if(binders && mean > th){
          c <- append(c,pp)
        }else if((!binders) && (mean < th)){
          c <- append(c,pp)
        }
      }else{
        sr <- sub$signal.ratio
        man <- sub$man.annotation
        if(length(unique(man)) == 1 && unique(man) == 1 && binders){
          c <- append(c,pp)
        }else if(length(unique(man)) == 1 && unique(man) == -1 && (!binders)){
          c <- append(c,pp)
        }else if(length(unique(man)) == 2){
          a <- isConsitent(pp,ll,sr,man)
          if(a == 1 && binders){
            c <- append(c,pp)
          }else if(a == -1 && (!binders)){
            c <- append(c,pp)
          }
        }
      }
    }
  }
  return(c)
}

splitDataIndividualLipidsBinders <- function(){
  data <- loadNewData()
  proteins <- unique(data$protein)
  lipids <- unique(data$lipid)
  lip.hom <- loadHomogeneousLipids()
  ids <- lip.hom$V1
  for(id in ids){
    print(id)
    name <- paste("proteins_binding_",as.character(lip.hom$V2[which(lip.hom$V1 == id)]),sep="")
    name.s <- paste("proteins.sequence_binding_",as.character(lip.hom$V2[which(lip.hom$V1 == id)]),sep="")
    name.d <- paste("proteins.domain.sequence_binding_",as.character(lip.hom$V2[which(lip.hom$V1 == id)]),sep="")
    
    id.lower <- tolower(id)
    proteins <- getBinders(binders=FALSE,ll=id.lower,criteria="man")
#    proteins.above.sr <- as.character(unique(subset(data$protein,subset=(data$lipid == id.lower & data$signal.ratio > th))))
#    proteins.below.sr <- as.character(unique(subset(data$protein,subset=(data$lipid == id.lower & data$signal.ratio < th))))
#    proteins.bind <- as.character(unique(subset(data$protein,subset=(data$lipid == id.lower & data$man.annotation == "1"))))
#    proteins.dont <- as.character(unique(subset(data$protein,subset=(data$lipid == id.lower & data$man.annotation == "-1"))))
    
#    proteins <- proteins.above.sr
    #WRITE THE PROTEIN NAMES
    write(x=proteins,file=paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/man.dont/names/",name,sep=""))
    #WRITE THE PROTEIN SEQUENCE  
    sequences <- sapply(proteins,getSequence)
    write(x=sequences,file=paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/man.dont/sequence.full/",name.s,sep=""))    
    #WRITE THE PROTEIN DOMAIN
    domains <- sapply(proteins,getDomainSequence)
    write(x=domains,file=paste("/Users/deghou/remote/projects/liposomes/newdata/motif.search/man.dont/sequence.domain/",name.d,sep=""))    
  }
}

getDomainSequence <- function(reserved.name){
  fasta.header <- prot.info$fasta_name[which(prot.info$reserved_name == reserved.name)]
  fasta.sequence <- prot.info$sequence[which(prot.info$reserved_name == reserved.name)]
  return(paste(fasta.header,fasta.sequence,sep="\n"))
}

getSequence <- function(reserved.name){
  uniprot <- as.character(prot.info$uniprot[which(prot.info$reserved_name == reserved.name)])
  seq.fa <- ""
  if(uniprot != "null"){
    # print(paste(reserved.name,uniprot))
    fasta.name <- subset(names(prot.sequences), grepl(tolower(uniprot),names(prot.sequences)))   
    # print(fasta.name)
    seq <- prot.sequences[[fasta.name]]
    seq.fa <- paste(paste(">",attr(seq,"name"),sep=""),as.character(seq[1]),sep="\n")     
    return(seq.fa)
  }else{
    return(NA)
  }
}


doCol <- function(){
  col.v <- c()
  a = 255
  for(k in 1:20){
    if(k == 1){
      a <- a - 15
    }else{
      a <- a - 10
    }
    col <- rgb(red=a,green=a,blue=255,max=255)
    print(col)
    col.v[k] <- rgb(col)
  }
  return(col.v)
}

addPredictionsToFile <- function(including.controls=FALSE,just.three=TRUE,including.syntaxine=FALSE){
  data <- loadNewData(including.controls=including.controls,just.three=just.three)
  threshold <- 0
  if(including.controls == TRUE){
    threshold = 0.02896941
  }else{
    threshold = 0.03720975
  }
  len <- length(data$V1)
  for(k in 1:len){
    n <- data[k,]
    if(data$V5[k] > threshold){
      n$V11 <- "PREDICTED.YES"
    }else{
      n$V11 <- "PREDICTED.NO"
    }
    write.table(file="remote/projects/liposomes/newdata/all.without.cc.without.bq.wihtout.d.q.withoutNA.withoutSyntaxine.withoutControls.without.45678.with.predictions",n,append=TRUE)
  }
}




