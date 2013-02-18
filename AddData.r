#!/usr/bin/r

new.data.dir <- "/Users/deghou/remote/projects/liposomes/newdata/per.runs/"
annotations.dir <- "/Users/deghou/remote/projects/liposomes/annotations/"
setwd(new.data.dir)
map.r.p <- setMapRunProtein()
list.of.new.runs <- list.files();
list.of.new.runs <- grep("run*",list.of.new.runs,value=TRUE)

#this file contains all the protein/lipids interacitons. It needs to be cleaned up afterwards !
final.file.with.120.lipids <- paste(new.data.dir,"final.file.with.120.lipids.NOT.CLEANED.UP.ALL.RUNS.2",sep="")
#write(file=final.file.with.120.lipids,"caca",append=TRUE)
print(final.file.with.120.lipids)
prot.round <- list()
good.runs <- loadGoodRuns()
#list.of.new.runs <- list.of.new.runs[1:10]
for(run in list.of.new.runs){
  print(run)
  if(run %in% good.runs){
    prot <- map.r.p[[run]]
    round <- "round"
    if(is.null(prot)){
      write(file="/Users/deghou/remote/projects/liposomes/newdata/WARNINGS",append=TRUE,paste("RUN",run,"IS NOT IN THE HASH NOR IN THE MAPPING FILE NOR IN IVANA'S EXCEL FILE. HENCE NOT TO BE TAKEN INTO CONSIDERATION"))
      paste("RUN",run,"IS NOT IN THE HASH NOR IN THE MAPPING FILE NOR IN IVANA'S EXCEL FILE. HENCE NOT TO BE TAKEN INTO CONSIDERATION")
    }else{
      prot <- as.character(prot)
      if(prot %in% names(prot.round)){
        prot.round[[prot]] <- as.numeric(as.character(prot.round[[prot]])) + 1
        round <- paste(round,prot.round[[prot]],sep="")
      }else{
        prot.round[[prot]] <- 1
        round <- paste(round,1,sep="")
      }
      #  if(prot == "NULL"){
      #    print(pate("CAREFULL ! ",run," has no protein name"))
      #  }
#      print(paste(prot,round))
      sum.data.file <- read.table(paste(run,grep("ratiolip*",list.files(run),value=TRUE),sep="/"),header=TRUE)
      run.manual.binding <- loadManualBinding(run)
      run.manual.binding <- tolower(run.manual.binding$V1)
      if(run.manual.binding != "NULL"){
        len <- length(sum.data.file$lipid)
        to.out <- ""
        for(k in 1:len){
          lipid <- tolower(as.character(sum.data.file$lipid[k]))
          mean <- tolower(as.character(sum.data.file$mean[k]))
          sd <- tolower(as.character(sum.data.file$sd[k]))
          n <- tolower(as.character(sum.data.file$n[k]))
          name <- tolower(as.character(sum.data.file$name[k]))
          nb <- tolower(as.character(sum.data.file$liposomes[k]))
          man <- 0
       #   print(lipid)
        #  print(run.manual.binding)
          if(lipid %in% run.manual.binding){
         #   print("CACACACACACA")
            man <- 1
          }else{
            man <- -1
          }
          to.out <- paste(to.out,paste(round,run,prot,lipid,mean,sd,n,name,nb,man,sep="\t"),sep="\n")
          #    write(file=final.file.with.120.lipids,paste(round,run,prot,lipid,mean,sd,n,name,nb,man,sep="\t"),append=TRUE)
        }
        #  print(file)
        write(file=final.file.with.120.lipids,to.out,append=TRUE)  
      }else{
        write(file="/Users/deghou/remote/projects/liposomes/newdata/WARNINGS",append=TRUE,paste("NO MANUAL BINDING ANNOTATIONS FOR",run))
        print(paste("NO MANUAL BINDING ANNOTATIONS FOR",run))
      }   
    } 
  }

}
