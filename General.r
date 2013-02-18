  #!/usr/bin/r


colLabSam  =  function(n) {
  if(is.leaf(n)) {
    a <- attributes(n)
    # clusMember - a vector designating leaf grouping
    # labelColors - a vector of colors for the above grouping
    labCol <- as.character(labelColors[which(names(labelColors) == a$label)])
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
}

#@dec.map = makeDecisionHeatMap()
printProtProfileBelowDendrogramHet <- function(hc){
  labelColors = loadManualColoring()
#  plot(dendrapply(as.dendrogram(hc),colLabSamHet),main = "Homogeneous lipid cluster by family")
  doDoDo(hc)
  order = hc$labels[hc$order]
  all.vaues = as.numeric(as.vector(mat.all))
  #  all.values = all.values[order(all.values)]
  col = topo.colors(length(all.values))
  x = -0.3
  
  for(lip in order){
    
    print(lip)
    x = x + 1
    #    profile <- getProteinProfileOfLipid(lip)
    #    profile = profile[order(profile)]
    profile <- dec.map[,lip]
    #    index.in.all.values.that.have.the.same.value = which(all.values %in% profile)[1:103]
    #    colors = col[index.in.all.values.that.have.the.same.value]
    if(x == 0.7){
      plotrix::color.legend(x,-2,x+ 0.7,-1.3,rownames(coop.col),rev(profile),cex=0.2,gradient="y",align="lt")       
    }else{
      plotrix::color.legend(x,-2,x+ 0.7,-1.3,c(""),rev(profile),cex=0.6,gradient="y",align="rb")  
    }
    text(x,0,length(which(profile == "indianred")),cex=0.5,col="blue")
  }
}



#@call
#hc = hclust(d = amap::Dist(t(mat.all.27), "pearson"))
printProtProfileBelowDendrogram <- function(hc){
  labelColors = loadManualColoring()
  plot(dendrapply(as.dendrogram(hc),colLabSam),main = "Homogeneous lipid cluster by family")
  order = hc$labels[hc$order]
  all.vaues = as.numeric(as.vector(mat.all))
#  all.values = all.values[order(all.values)]
  col = topo.colors(length(all.values))
  x = -0.3
  
  for(lip in order){
     
    print(lip)
    x = x + 1
#    profile <- getProteinProfileOfLipid(lip)
#    profile = profile[order(profile)]
    profile <- coop.col[,lip]
#    index.in.all.values.that.have.the.same.value = which(all.values %in% profile)[1:103]
#    colors = col[index.in.all.values.that.have.the.same.value]
    if(x == 0.7){
      plotrix::color.legend(x,-1.5,x+ 0.7,-.8,rownames(coop.col),rev(profile),cex=0.2,gradient="y",align="lt")       
    }else{
      plotrix::color.legend(x,-1.5,x+ 0.7,-.8,c(""),rev(profile),cex=0.6,gradient="y",align="rb")  
    }
    text(x,0.5,length(which(profile == "indianred")))
    }
  }


getProteinProfileOfLipid <- function(lipid,sr=TRUE,man=FALSE){
  if(sr){
#    mat.all = doMatrixOutOfScreening()
    return(as.numeric(mat.all[,which(colnames(mat.all) == lipid,useNames=FALSE)]))
  }else{
    
  }
}

getHomMatWithName <- function(m){
  m1 <- reduceMatrix(m,loadHeterogeneousLipids()$V1,col=TRUE)
  m2 <- t(replaceLipidByNames(t(m1)))
  return(m2)
}



addProteinConcentrationToMainFile <- function(){
 data <- loadNewData()
data["protein.concentration"] <- seq(1,24797)
 prots = unique(data$protein)
 runs = unique(data$run)
 for(run in runs){
   conc <- getProteinConcentration(run)
   print(paste(run,conc))
   data$protein.concentration[which(data$run == run)] = conc
 }
 return(data)
}

makeDecForProtAndLip <- function(prot,lip){
  sr = subset(data$signal.ratio,subset = data$protein == prot & data$lipid.composition == lip)
  man = subset(data$man.annotation,subset = data$protein == prot & data$lipid.composition == lip)

  if(length(sr) == 0 || length(man) == 0){
    return("unclear")
  }else{
    dec = makeDecision(man,sr)
    return(dec) 
  }
}

docacaca <- function(){
  data <- loadNewData()
  protein <- unique(data$protein)
  lipids <- unique(data$lipid.composition)
#  hom <- tolower(as.character(loadHomogeneousLipids()$V2))
  hom = hc$labels[hc$order]
  mat = matrix(data=-1,nrow=length(protein),ncol=length(hom),dimnames=list(protein,hom))
  for(p in protein){
    print(p)
    for(l in hom){
      dec = makeDecForProtAndLip(p,l)
    
      if(dec == "bind"){
        print(p)
        mat[p,][l] <- 1
      }else if(dec == "unclear"){
        mat[p,][l] <- 0
      }
    }
  }
  return(mat)
}


makeDecision <- function(annotations,ratios){
  th = 0.03720974
  prop.pos = length(which(annotations == 1)) / length(annotations);prop.neg = 1 - prop.pos;prop.above.th = length(which(ratios > th)) / length(ratios);prop.below.th = 1 - prop.above.th
  maj.above = (prop.pos > 0.5); maj.below = (prop.neg > 0.5); maj.th.above = (prop.above.th > 0.5); maj.th.below = (prop.below.th > 0.5) 
  
  if((maj.above & maj.th.above) || (maj.above & prop.above.th == 0.5)){
    return("bind")
  }else if((maj.below & maj.th.below) || (maj.above & prop.below.th == 0.5)){
    return("does not bind")
  }else if(prop.pos == 0.5 & maj.th.above){
    return("bind")
  }else if(prop.neg == 0.5 & maj.th.below){
    return("does not bind")
  }else{
    return("unclear")
  }
  
}


addPHsuffix <- function(x){
  f <- function(x){
    
    return(paste(x,"-ph",sep=""))
  }
  return(sapply(x,FUN=f,USE.NAMES=FALSE))
}

getRidOfPHsuffix <- function(x){
  f <- function(x){
    return(substr(x,0,match("-",strsplit(x,"")[[1]]) - 1))
  }
  return(sapply(x,FUN=f,USE.NAMES=FALSE))
}


getSimpleProteinName <- function(vector){
  f <- function(x){
    return(strsplit(x,split="\\|")[[1]][2])
  }
  c <- sapply(vector,FUN=f,USE.NAMES=FALSE)
  return(c)
}

fromLipidGet <- function(lipid,what){
  f <- function(x){
    if(startsWith(x,what)){
      return(x)
    }
  }
  alignment <- loadMultipleAlignment(lipid)
  c <- sapply(alignment$nam,f,USE.NAMES=FALSE)
  c <- as.character(c)
  c <- c[-grep("NULL",c)]
  return(c)
}


startsWith <- function(str,s){
  return(substr(str,0,nchar(s) ) == s)
}


getLipidId <- function(lipid.name){
  lip.het <- loadHeterogeneousLipids()
  lip.hom <- loadHomogeneousLipids()
  all <- append(lip.het,lip.hom)
  df <- data.frame()
  df <- rbind(df,lip.het)
  df <- rbind(df,lip.hom)
  id <- tolower(as.character(df$V1[which(df$V2 == lipid.name)]))
  return(id)
}

getLipidName <- function(lipid.name){
  lip.het <- loadHeterogeneousLipids()
  lip.hom <- loadHomogeneousLipids()
  all <- append(lip.het,lip.hom)
  df <- data.frame()
  df <- rbind(df,lip.het)
  df <- rbind(df,lip.hom)
  id <- tolower(as.character(df$V2[which(df$V1 == lipid.name)]))
  return(id)
}



testColM = function(k){
  colfunc <- colorRampPalette(c("gray67", "white"))
  plot(rep(1,k),col=colfunc(k),pch=19,cex=3)
}
#imprimer  heatmap avec s = 0.5 et k = 6



#pour la deuxiemen matrice : garder 0.03 et 50 (enfin on va voir les valeurs des deux autre sproteines ... )
generateHeatmaps <- function(s,k,z,y){
#  mat <- m  
  mat = mat.12
#  c = c("white")
  c <- colorpanel(y + 2,"cornflowerblue","goldenrod")
#  c <- append(c,"orange")
#  int <- interval(0.037,0.1,n=50)
#  for(s in int){
#print(s)
 #   dir1 = paste("/Users/deghou/Desktop/manu/",s,sep="")
#    dir.heatmaps <- paste(dir1,"/heatmaps",sep="")
#    dir.legends <- paste     (dir1,"/legends",sep="")
#    cmd1 <- paste("mkdir",dir1)
#    cmd2 <- paste("mkdir",dir.heatmaps)
#    cmd3 <- paste("mkdir",dir.legends)
#    system(cmd1)
  ##    system(cmd2)
#    system(cmd3)
 #   for(k in 15:100){
 #     print(k)
  mat.vec = as.vector(mat)
#      pdf(file=paste("/Users/deghou/Desktop/manu/",s,"/heatmaps/heatmap",k,".pdf",sep=""),width=1height=20)
#      B <- c(interval(min(mat.vec),s,k) , interval( s,max(mat.vec),(210 - k)))
    B <- c(interval(min(mat.vec),s,k) , interval( s,z,y))
#  B <- c()
#      B <- B[-(k+1)]
#          B <- B[-(k+1)]
      heatmap.2(mat,breaks=B,col=c,Colv=FALSE,Rowv=FALSE,trace="none",margins=c(10,30),sepwidth=c(2,2),sepcolor="black")
#      dev.off()
#      pdf(file=paste("/Users/deghou/Desktop/manu/",s,"/legends/heatmap",k,".legend.pdf",sep=""),width=10,height=20)
 #     plot.new()
#      B <- signif(x=B,digits=2)
#      plotrix::color.legend(0,-0.1,0.1,1.1,legend=B,rect.col=c,gradient="y",cex=0.3)
#      dev.off()
#    }
#  }
}

characterOccurence <- function(seq){
  s <- unique(seq)
  occurence <- c()
  df <- data.frame("s","occurence")
  ind <- 0
  a <- c()
  b <- c()
  for(k in s){
    g <- length(grep(pattern=k,seq))
    occurence <- append(occurence,g,after=TRUE)
    a <- append(a,k)
    b <- append(b,g)

  }
  df <- data.frame(a,b)
  names(df) = c("s","occurence")
  return(df)
}


replaceCharacterInString <- function(string,character,value){
  ind <- grep(character,strsplit(string,"")[[1]])
  for(k in ind){
    substr(string,k,k) <- value 
  }  
  return(string)
}

adaptBreaksColEasy <- function(mat,minValue,maxValue,low,mid,high){
  b <- sort(unique(mat))
  col <- colorpanel(length(b),low,mid,high)
}


adaptBreaksCol <- function(mat,minValue,maxValue,low,mid,high){
  l <- length(subset(manu,subset=manu > minValue & manu < maxValue))
  b <- c()
  b[1] <- 0
  b[2] <- minValue
  b <- append(b,sort((subset(manu,subset=manu > minValue & manu < maxValue))))
  b[length(b)] <- maxValue  
  col <- c()
  col[1] <- low
  cc <- colorpanel(l + 10 , low,mid,high)[10:l]
  col <- append(col,cc)
  col[length(col)] <- high
}

  colLab <- function(n) {
    if(is.leaf(n)) {
      a <- attributes(n)
    #  print(a$label)
      # clusMember - a vector designating leaf grouping
      # labelColors - a vector of colors for the above grouping
      if(a$label == "spo71c-ph" || a$label == "spo71n-ph"){
        labCol <- "blue";
      }else if(a$label == "osh3-ph" || a$label == "osh3ct-ph"){
        labCol <- "green";
      }else if(a$label == "bud4ct_fused-ph" || a$label == "bud4ct_1-ph"){
        labCol <- "pink";
      }else if(a$label == "fyve_eea1-ph" || a$label == "fyve_eea1_short-ph"){
        labCol <- "red";
      }else if(a$label == "rtt106-ph" || a$label == "rtt106ct-ph"){
        labCol <- "brown";
      }else if(a$label == "bem3-ph" || a$label == "bem3ct-ph"){
        labCol <- "cyan";
      }else if(a$label == "spo14-ph" || a$label == "spo14ct-ph"){
        labCol <- "purple";
      }else if(a$label == "caf120_2-ph" || a$label == "caf120_fused-ph"){
        labCol <- "orange";
      }else{
        labCol = "black"
      }
      print(labCol)
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
  }




tmpFun <- function(proteins){
  prot.info <- loadProteinInformation()
  unips <- sort(unique(prot.info$uniprot))
  col <- c()
  for(prot in proteins){
    unip <- as.character(prot.info$uniprot[which(prot.info$reserved_name == prot)])
    col <- append(col,which(unips == unip))
  }
  return(col)
}
getLeavesFromCluster <- function(hc){
  return(hc$labels[as.vector(hc$order)])  
}

showpanel <- function(col){
  image(z=matrix(1:100, ncol=1), col=col, xaxt="n", yaxt="n" )
}

makeAssociativeArry <- function(ind_key,ind_value,df){
  l = list()
  for(k in 1:length(df[[ind_key]])){
    key = as.character(df[[ind_key]][k])
    value = as.character(df[[ind_value]][k])
    l[[key]] <- value
  }
  return(l)
}

extractMapToVectors <- function(map){
	a <- c()
	b <- c()
	c1 <- paste(map,names(map),sep="___")
	for(k in c1){
		#print(strsplit(k,"___")[[1]][1])
		a <- append(a,strsplit(k,"___")[[1]][1])
		b <- append(b,strsplit(k,"___")[[1]][2])
	}
	l <- list()
	print(length(a))
	print(length(b))
	l$a <- a
	l$b <- b
	return(l)
	
}


buildIndFeatColPalette <- function(){
	ind.feat <- read.table("/Users/deghou/remote_copy/projects/liposomes/data/proteins_restored/sub_cellular_localisation",sep="\t")$V1
	ind.feat.col.mapping <- list()
	c <- 1
	for(k in ind.feat){
		c <- c + 1;
		ind.feat.col.mapping[[tolower(k)]] <- c
	}
	return(ind.feat.col.mapping)
}



buildUnique <- function(vec,sep){
	lipid_id_color_vec <- list()
	col_vec <- c()
	for(k in 1:length(vec)){
	lipid_mixture <- strsplit(as.character(vec[k]),sep)[[1]]
	lipid_mixture <- unique(lipid_mixture)
	for(z in 1:length(lipid_mixture)){
	    li <- as.character(lipid_mixture[z])
	    
	    if(nchar(li) > 2){
#		print(col_vec)
	     	col_vec <- append(col_vec,li)
	    }
	   }
#	   print(length(col_vec))
	 }
	return(unique(col_vec))
	
}

reduceMatrix <- function(mat,v,col=FALSE,row=FALSE,discard=TRUE){
	m <- mat;
  w <- c()
  if(col){
    w <- colnames(mat)
  }else{
    w <- rownames(mat)
  }
  if(discard){
    print("will discard everything that is found in the vector")
    for(k in 1:length(v)){
      if(row){
        i <- which(tolower(rownames(m)) == tolower(v[k]));
    #    print(i)
        m <- m[-i,];
      }else if(col){
        i <- which(tolower(colnames(m)) == tolower(v[k]));
        m <- m[,-i]
      }
    } 
  }else{
    print("will discard everything that is not found in the vector")
    for(k in 1:length(w)){
      if(!w[k] %in% v){
        i <- which(rownames(m) == w[k])
        print(paste(w[k],"has to be removed"))
        if(row){
          m <- m[-i,]
        }else{
          m <- m[,-i]
        }
      }
    }
  }
	return(m)
}


writeMatriceToPairwiseFile <- function(m, f){

 for(k in 1:length(rownames(m))){
	a <- as.character(rownames(m)[k])
	aa <- strsplit(a,"\\|")[[1]][1]
  print(aa)
	for(l in k:length(colnames(m))){
		b <- as.character(colnames(m)[l])
		bb <- strsplit(b,"\\|")[[1]][1]		
		v <- m[a,][b]
		s <- paste(aa,bb,v,sep="\t")
		
		write(s,file=f,append=TRUE)
	}
 } 
}

getCurves <- function(df,cutoff){
  pred <- df[1]
  lab <- df[2]
  tp <- c()
  fp <- c()
  fn <- c()
  tn <- c()
  for(k in 1:length(cutoff)){
    ctf <- cutoff[k]
    predicted.positive <- length(which(pred > ctf))
    true.positive <- length(subset(length(which(pred > ctf)),lab == 1))
    false.positive <- length(subset(length(which(pred > ctf)),lab == -1))
    false.negative <- length(subset(length(which(pred < ctf)),lab == 1))
    true.negative <- length(subset(length(which(pred < ctf)),lab == -1))
    tp[k] <- true.positive
    fp[k] <- false.positive
    fn[k] <- false.negative
    tn[k] <- true.negative
  }
  plot(cutoff,tp)
  lines(cutoff,tp,col="red")
  lines(cutoff,fp,col="green")
  lines(cutoff,tn,col="gray")
  lines(cutoff,fn,col="blue")
}

getTp <- function(df,threshold){
  pred <- df[1]
  lab <- df[2]
  return(length(which(subset(lab,pred > threshold) == 1)))
}

getFn <- function(df,threshold){
  pred <- df[1]
  lab <- df[2]
  return(length(which(subset(lab,pred < threshold) == 1))) 
}


interval <- function(min,max,n){
  inc <- (max - min)/n
  b <- c()
  b <- append(b,min)
  a <- min 
  for(k in 1:n){
  	a <- a + inc
  	b <- append(b,a)
  }
  return(b)
}


fillDiagonalWith <- function(m,value){
	c <- rownames(m)
	for(k in c){
		m[k,][k] <- value
	}
	return(m)
}

getProteinsInTreeCut <- function(hc,k=NULL,h=NULL){
  groups = cutree(hc, k=k,h=h)
  groups = sort(groups)
  res <- list()
  proteins <- names(groups)
  for(protein in proteins){
  #  print(res)
    value <- as.numeric(groups[protein])
    l <- length(res)
 #   print(paste("length =",l,"value =",value)) 
    if(value > l){
      res[[value]] <- c(protein);
    }else{
      c <- res[[value]]
      c <- append(c,protein)
      res[[value]] <- c
    }
  }
  return(res)  
}

drawPieFromList <- function(l){
  c1 <- c()
  for(k in l){
    c1 <- append(c1,k)
  }
  return(c1)
}

  drawPieEnrichment <- function(proteins,what,ncol){
    pie(table(setPieEnrichment(proteins,what)$en),col=rainbow(ncol),cex=1)
  }
  
  setPieEnrichment <- function(proteins,what){
    prot.info <- loadProteinInformation() 
    enrichment <- NULL
    if(what == "localisation"){
      enrichment <- read.table("/Users/deghou/remote/projects/liposomes/newdata/enrichment/prot.localisation",sep="\t")
    }else if(what == "function"){
      enrichment <- read.table("/Users/deghou/remote/projects/liposomes/newdata/enrichment/prot.function",sep="\t")
    }else if(what == "process"){
      enrichment <- read.table("/Users/deghou/remote/projects/liposomes/newdata/enrichment/prot.process",sep="\t")
    }
    en <- c()
    for(protein in proteins){
    #  print(protein)
      unip <- as.character(prot.info$uniprot[which(prot.info$reserved_name == protein)])
      if(unip != "null"){
        print(unip)
        enn <- as.character(enrichment$V2[which(enrichment$V1 == unip)])
      #  print(enn)
        en <- append(en,enn)
      }
    }
    u <- unique(en)
    uu <- c()
    for(k in en){
      uu <- append(uu,which(u == k))
    }

    return(data.frame(en,uu))
  }

  
setPieForGroupOfProtein <- function(proteins,col=1){
  l <- loadProteinInformation()
  
  ratios <- list()
  #for each protein
  for(protein in proteins){
    index = which(l$reserved_name == protein)
    if(col == 1){
      adjacent_domains = as.character(l$adjacent_domain_profile[index]) 
    }else if(col == 2){
      adjacent_domains = as.character(l$ffunction[index])       
    }else if(col == 3){
      adjacent_domains = as.character(l$process[index]) 
    }else if(col == 4){
      adjacent_domains = as.character(l$localisation[index]) 
    }
    individuals = strsplit(adjacent_domains,";")
    #get the adjacent domain
    for(ind in individuals[[1]]){
      if(!is.numeric(ratios[[ind]])){
        ratios[[ind]] <- 1
      }else{
        ratios[[ind]] <- ratios[[ind]] + 1
      }
    }
  }
  
  c <- c()
  for(k in names(ratios)){
    print(k)
   c <- append(c,ratios[[k]])
  }
  
  return(ratios)
}


transformMatrixToTabFile <- function(m,output_file){
if(!is.matrix(m)){
 m <- as.matrix(m);
}
vars <- dimnames(m)[[2]];
#dimnames(m) <- list(vars,vars);
print(vars)
for(i in 1:length(vars)){
        print(paste(i , " : " , vars[i]));
        var1 <- vars[i];
	k <- i + 1;
        for(j in k:length(vars)){
                var2 <- vars[j];	
		print(var1)
                value <- m[var1,][var2];
                s <- paste(var1,var2,value,sep="\t");
                write(s, file = output_file, append = TRUE, sep = "\n");
        }
  }
}

getClusterLeavesXaxis <- function(hc, level=length(hc$height), init=TRUE){   

  if(init){
	.count <<- 0
    	.topAbsis <<- NULL
    	.heights <<- NULL


  }

  if(level<0) {
    .count <<- .count + 1
    return(.count)
  }

  node <- hc$merge[level,]
  le <- abs(hc, node[1], init=FALSE)
  ri <- abs(hc, node[2], init=FALSE)
 

  mid <- (le+ri)/2  

  .topAbsis <<- c(.topAbsis, mid)
  .heights <<- c(.heights, hc$height[level])  

  invisible(mid)
  return(.topAbsis)
} 



reverseList <- function(l){
  ll <- list()
  c <- names(l)
  for(k in c){
    if(nchar(l[""]))
    v <- c[k]
    
  }
}

