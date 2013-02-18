#/usr/bin/r


createProtFeatureCompositionColorPalette <- function(df,prot_col_index=NULL,feature_col_index=NULL,sep=";",directly=FALSE,manual=FALSE){
  feature.ind.color.palette <- list()
  prot.feature.composition.color.palette <- list()
	prot.info <- df
  if(manual == FALSE && directly == FALSE){
    ll <- loadProteinInformation()
    print(dim(ll))
    feature.ind.col.palette <- createIndFeatures(ll,8,1,";",FALSE,FALSE)
    print((paste(length(feature.ind.col.palette)," features have been loaded")))
  } 
	if(manual == TRUE){
		feature.ind.color.palette <- loadManualColoringLocaton()
		print((paste(length(feature.ind.color.palette)," features have been loaded")))
	#	feature.ind.color.palette <- createIndFeatures(df,prot_col_index,feature_col_index,sep,directly,manual)
	}

	if(directly == TRUE){
		ii <- length(df)
	}else{
		ii <- length(prot.info[[prot_col_index]])
	}
	for(k in 1:ii){
		if(directly == TRUE){
			protein <- tolower(as.character(df[k]))
			features <- strsplit(as.character(df[k]),sep)
			features[[1]] <- unique(features[[1]])
		}else{
			protein <- tolower(as.character(prot.info[[prot_col_index]][k]))
			features <- strsplit(as.character(prot.info[[feature_col_index]][k]),sep)
			features[[1]] <- unique(features[[1]])
		}
#		print(features)
		c <- c()
	#	print(paste(protein, " --->", features))
		for(i in features[[1]]){
			i <- tolower(i)	
			if(nchar(i) > 2){
				col <- feature.ind.col.palette[[i]]
#				print(paste(protein,i,col))
#				print(paste(protein," ----> ",features[[1]], "(", col , ")"))
				if(is.null(col)){
					print(paste(i , " does not have a color ", "( for ", protein , " )"));
				}else{
					c <- append(c,col)
				}
			}
		}
  #  print(paste(protein , "--->", c))
		prot.feature.composition.color.palette[[protein]] <- c		
#    print(prot.feature.composition.color.palette[[protein]])
    print(length(prot.feature.composition.color.palette))
	}
	#print(feature.ind.color.palette)
	return(prot.feature.composition.color.palette)
}




createIndFeatures <- function(df,prot_col_index=NULL,feature_col_index=NULL,sep=";",directly=FALSE,manual=FALSE){
  print(paste(prot_col_index,feature_col_index,sep,directly,manual))
  feature.ind.color.palette <- list()
        prot.feature.composition.color.palette <- list()
        prot.info <- df

        if(directly == TRUE ){
                if(manual == TRUE){
			print("will load pre existing")
                        feature.ind.color.palette <- loadManualColoringLocaton()
		
                }else{
                        u <- buildUnique(df,sep)
                        for(k in 1:length(u)){
                        #print(u[k])
                                feature.ind.color.palette[[as.character(u[k])]] <- k
                        }
                }
        }else{

                if(sep != ";"){
                        u <- prot.info[feature_col_index]
                        u <- unique(u)

                        for(k in 1:length(u[[1]])){
                #print(u[k])
                                feature.ind.color.palette[[as.character(u[[1]][k])]] <- k
                        }
                }else{
                        u <- buildUnique(prot.info[[feature_col_index]],sep)

                        for(k in 1:length(u)){
                #print(u[k])
                                feature.ind.color.palette[[as.character(u[k])]] <- k
                        }

                }
        }
	return(feature.ind.color.palette)
}





plotLabelColorsToClusterPlot <- function(cluster,lab_mixture,title,width,height=1){
	d <- as.dendrogram(cluster);
	y <- max(cluster$height) - (max(cluster$height)/(3/height))
	height <- max(cluster$height) / 5
	dev.new(width=18,height=9.3)
	par(cex=0.8)
	plot(cluster,main=title,cex.lab=3)
	dd <- "num1ct-ph"
	print(lab_mixture[[dd]])
	labels <- labels(d)
	i <- 0;
	col.labels <- c()
	for(lab in labels){
		i <- i + 1
		cv <- lab_mixture[[as.character(lab)]]
		print(cv)
		if(is.null(cv)){
			print(paste("problem with ",lab))
		}else{
#			print(paste(lab, " ----> ", cv))
#			plotrix::color.legend(i,y,i+width,y+height,col.labels,cv,gradient="y")
		}
		
	}
}

