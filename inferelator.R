###################################################
# Inferelator version B, fewer false negatives
###################################################

###########################################
# preclust.and.inferelate
###########################################

preclust.and.inferelate <- function( clusterStack, predictors, data, predictor.mat=NULL, preclust=30, filter.by.aic=5, col.map=NULL, tau=10,  funcs="min", cutoff=0.01, cv.k=10, cv.choose=c("min","min+1se"), logit.scale=0.25, quiet=F, kmeans.iter=200, kmeans.nstart=25, scale.coeffs=T, plot=F){

# clean up Cytoscape files to make sure each run is fresh
unlink("*.noa")
unlink("*.sif")
unlink("*.eda")

# make sure all the TFs given exist in the ratios given
predictors <- predictors[predictors %in% rownames(data)]
predictors.genetic <- grep("VNG", predictors, value=T)
predictors.env <- setdiff(predictors, predictors.genetic)
# normalize the ratios of the genes to mean 0 and sd 1
data <- mean.variance.normalize(data)
# slice out the ratios of the supplied predictors
tf.matrix <- data[ predictors, ]

	if ( ! is.na( preclust ) && preclust != 0 ) {
	# "precluster" ONLY THE GENETIC predictors into TFGROUPs (don't want to mix environmental factors)
	# but append the environmental ones so we can have env+TFGROUP combos later	
	predictor.mat <- preclust(data = data, tf.matrix = data[predictors.genetic,] , clust.count = preclust, kmeans.iter = kmeans.iter, kmeans.nstart = kmeans.nstart)
		if ( length(predictors.env) > 0){
		genetic.names <- rownames(predictor.mat)	
		# add the real -environmental- predictors to predictor matrix...
		predictor.mat <- rbind(predictor.mat, data[predictors.env,])
		rownames(predictor.mat) <- c(genetic.names,predictors.env)		
		}
	}

	if ( ! is.na(funcs) && length(funcs) > 0 ){
	# make combinatory predictors from TFGROUPs 
	predictor.mat.ands <- combine.with.ands(data=predictor.mat, predictors = rownames(predictor.mat) , funcs = funcs)
	}

	results<-list()

	for ( i in 1:length(clusterStack) ) {
		k <- clusterStack[[ i ]]$k
		predicted <- clusterStack[[i]]$rows
		predictor.mat <- tf.matrix
		cat(paste("\n Working on bicluster #", k, "using", length(predictors), "predictors\n\n"))
		
		if ( ! is.na( filter.by.aic ) && filter.by.aic != 0 ){
		# filter single TFGROUPs and combinatory TFGROUPs -separately- by AIC, keeping the best in both cases
		best.singletons <- filter.by.aic(data = data, predicted = predicted, predictor.matrix = predictor.mat, top.aic.to.keep = filter.by.aic)
		best.combos <- filter.by.aic(data = data, predicted = predicted, predictor.matrix = predictor.mat.ands, top.aic.to.keep = filter.by.aic)
		}
		
		predictor.mat <- rbind(best.singletons, best.combos)
		# add the -environmental- predictor's profiles to use as starting points for hidden variable
#		if ( length(predictors.env) > 0){
#		hidden.mat <- data[predictors.env,]
#		rownames(hidden.mat) <- paste(predictors.env, ".hidden", sep="")
#		predictor.mat <- rbind(predictor.mat, hidden.mat)			
#		}
		
	
		# SET PREDICTOR MATRIX TO A GLOBAL VARIABLE (FOR DEBUGGING PURPOSES ONLY! ALSO, WILL OVERWRITE SO ONLY LAST CLUSTER'S IS INVOLVED)
		predictor.mat.debug <<- predictor.mat
		
		# this line actually runs Inferelator on the single biclust
		singleresult <- inferelator(predicted=predicted, predictors=predictors, data=data, predictor.mat=predictor.mat, col.map=col.map, tau=tau, funcs=funcs, cutoff=cutoff, cv.k=cv.k, cv.choose=cv.choose, logit.scale=logit.scale, scale.coeffs=scale.coeffs, quiet=quiet, plot=plot)
		
		# if there are non-zero coefficients, add result for single cluster to the bigger resultset
		if (length(singleresult) > 0){
			names( singleresult ) <- paste( k, names( singleresult ), sep="." )
			results[[ k ]] <- singleresult
			#print(paste("got non-zero results, adding to", k))
			mean.profile <- apply(data[predicted, ], 2, mean)
			if( plot ) {
				matplot(t(rbind( mean.profile , predictor.mat[sub("[0-9]+.","",names(singleresult)),])), type="l", col=c("red",rep("grey",length(singleresult))), main=paste("Bicluster", k, " ", clusterStack[[i]]$nrows, "genes"), ylab="Normalized expression", xlab="Conditions")
				legend("bottomright", paste(c( "mean biclust profile", sub("[0-9]+.","",names(singleresult))), " ", c("",substr(singleresult,1,5)) ), col=c("red", rep("grey",length(singleresult))), lty=c(1,1:length(singleresult)))
				#legend("bottomleft", clusterStack[[i]]$rows)
			}
		} else {
		# if there are null coefficients (null = didn't meet cut-off), still return a cluster number so we know it was processed	
			results[[ k ]] <- k	
			#print(paste("no results ! k is", k))
		}
		
}
return(results)		

} # end of preclust.and.inferelate function


###########################################
# preclust
###########################################
preclust <- function(data, tf.matrix, clust.count, kmeans.iter, kmeans.nstart) {
## cluster tf's w/ kmeans before feeding into LARS
## otherwise too many predictors and LARS will fail 
## manually set seed so will get same clusters each time 
   		
   	set.seed(31337)
   	
   	data.c <- kmeans( tf.matrix, clust.count, iter.max=kmeans.iter, nstart=kmeans.nstart )	
   	result <- data.c$centers
   	 		
   	## name the kmean-clustered centers
   	rownames( result ) <- paste( "TFGROUP", 1:clust.count, sep="" )
   	
   	# if there are files already, don't overwrite them, since each preclust() run will result
   	# in the same preclusts (see above), no need to write if not necessary
   	if ( ! file.exists("tfGroups.noa") && ! file.exists("types.noa")) {
   		write("tfGroups", "tfGroups.noa")
   		write("type", "types.noa")
   		for ( i in 1:length(data.c$size)){
    		write( paste(paste("TFGROUP", i, " = ", "(", sep="") ,paste(names(which(data.c$cluster==i)), collapse="::"), ")", sep="") , "tfGroups.noa", append=T)
    		write(paste("TFGROUP", i, " = (regulator)", sep=""), "types.noa", append=T)	
   		}
   	} 
	
cat( "Preclusted with k-means, predictor matrix is", nrow( result ), "x", ncol( result ), "\n" )	
return(result)
} #end of preclust function


filter.by.aic <- function(data, predicted, predictor.matrix, top.aic.to.keep){
#this function takes a matrix of ratios and using AIC, keeps only the top <specified number>
# i.e.  ones with lowest AIC's , throwing out the rest	

# calculate mean profile of the cluster
	mean.profile <- apply(data[predicted, ], 2, mean)	
	best.preclusts<-numeric()
		
	for (j in 1:nrow(predictor.matrix)){
		paic <- numeric()
 		paic[j] <- AIC(lm(mean.profile ~ as.numeric(predictor.matrix[j,])))
 		names(paic)[j]<- rownames(predictor.matrix)[j]
 		#cat(paste("working on predictor", rownames(predictor.matrix)[j], "which has an AIC of", paic[j], "\n"))
 		best.preclusts<-c(best.preclusts, paic[j])
	}
	
	# keep only the best <specified number> of predictors

	best.preclusts <- names(sort(unlist(best.preclusts)))[1:top.aic.to.keep]
	result<- predictor.matrix[best.preclusts,]

	cat( "Filtered by AIC, predictor matrix is now ", nrow( result ), "x", ncol( result ), "\n" )	
	return(result)	
} #end of filter.by.aic function


###########################################
# combine.with.ands
###########################################
combine.with.ands <- function(data, predictors, funcs){
	## Make array of min()s and max()s:
    if ( ! is.null( funcs ) && ! is.na( funcs ) ) {
      for ( func in funcs ) {
        tmp <- sapply( predictors, function(i) t( sapply( predictors, function(j) apply( data[ c( i, j ), ], 2, func ) ) ), simplify=F )
        #names( tmp ) <- predictors
        tmp.2 <- sapply( predictors, function(i) rownames( tmp[[ i ]] ) <<- paste( i, predictors, func, sep="." ) )
        tmp.2 <- tmp[[ 1 ]]; for ( i in 2:length( tmp ) ) tmp.2 <- rbind( tmp.2, tmp[[ i ]] )
        tmp <- tmp.2[ ! 1:nrow( tmp.2 ) %in% sapply( paste( predictors, predictors, func, sep="." ),
                                                    grep, rownames( tmp.2 ) ), ]
        ## Combine full predictor mat (individual tfb's and their min/max pairs)
        ## Remove e.g. B.A.min if there is an A.B.min
        result <- unique(tmp)       
      }
    }
	
cat( "Combined with ands, predictor matrix is now ", nrow( result ), "x", ncol( result ), "\n" )
return(result)
} #end of combine.with.ands function


###########################################
# inferelator
###########################################
inferelator <- function( predicted, predictors, data, predictor.mat=NULL, col.map=NULL, tau=10, tau.best=NULL , funcs="min", cutoff=0.01, cv.k=10, cv.choose=c("min","min+1se"), logit.scale=0.25, quiet=T, scale.coeffs=T, plot=F, ... ){
# this function actually runs lars and output coefficients	
	
	predicted <- predicted[ predicted %in% rownames( data ) ]
  	predictors <- predictors[ predictors %in% rownames( data ) ]
  	predicted <- predicted[ ! predicted %in% predictors ]
  
	out.tmp <- apply( data[ predicted, ,drop=F ], 2, mean )
  	in.tmp <- predictor.mat

  cat(paste("Tau is currently", tau, "\n"))
  # Convert predicted variable to the difference equation, i.e. from y to tau/dt * (y_i - y_(i-1)) + y_(i-1)
  if ( ! is.null( col.map ) ) {
  	cat("Time series data supplied, converting predictors to difference equation... \n\n")
    conds <- colnames( data )
    good.i <- ( ( col.map[ conds, "isTs" ] == TRUE ) &
                ( col.map[ conds, "is1stLast" ] == "m" | col.map[ conds, "is1stLast" ] == "l" ) ) |
                ( col.map[ conds, "isTs" ] == FALSE & col.map[ conds, "is1stLast" ] == "e" )
    prevs <- col.map[ conds, "prevCol" ]
    del.ts <- col.map[ conds, "delta.t" ]
    del.ts[ del.ts < 1 ] <- 1
    out.tmp <- ( ( tau / del.ts ) * ( out.tmp - out.tmp[ prevs ] ) ) + out.tmp[ prevs ]
    out.tmp[ out.tmp > 3.0 ] <- 3.0
    out.tmp[ out.tmp < -3.0 ] <- -3.0
    in.tmp <- in.tmp[ ,prevs ]

  } else { cat("NO time series data supplied, assuming STEADY-STATE... \n\n") }
  
  ##predictor.mat <- rbind( data[ predicted, ], predictor.mat )
  ##colnames( in.tmp ) <- colnames( out.tmp )
  predictor.mat <- rbind( out.tmp, in.tmp )
  rownames( predictor.mat )[ 1 ] <- predicted[ 1 ]
   ########## FIND PARAMS USING LARS

  ##if ( ! quiet ) cat( predicted, ":\n" )
  
  df.tmp <- t( in.tmp ) ##predictor.mat[ predictors, ] )
  output <- as.numeric( out.tmp ) ##predictor.mat[ predicted[ 1 ], ] ) ## Linear output
  ##output <- plogis( predictor.mat[ predicted, ], scale=logit.scale ) ## Logistic output (how to determine offset?)
  if ( plot ) {
    par( mfrow=c( 2, 2 ) )
    ##plot( output, plogis( output ), scale=logit.scale, main=paste( predicted, collapse=" ", sep=" " ) )
  }

  ## Run lars on predicted ~ predictors
  require( lars )
  lars.tmp <- try( lars( df.tmp, output, type="lasso", ... ), silent=quiet )
  if ( class( lars.tmp ) == "try-error" ) {
    tries <- 1; while( tries <= 20 && class( lars.tmp ) == "try-error" ) {
      lars.tmp <- try( lars( df.tmp, output, type="lasso", ... ), silent=quiet )
      tries <- tries + 1
    } }
  if ( class( lars.tmp ) == "try-error" ) return( numeric() )
  if ( plot ) plot( lars.tmp, main=paste( predicted, collapse=" ", sep=" " ) )

  ## Do lars CV and ...
  cv.lars.tmp <- try( cv.lars( df.tmp, output, K=cv.k, type="lasso", plot.it=plot ), silent=quiet )
  if ( class( cv.lars.tmp ) == "try-error" ) {
    tries <- 1; while( tries <= 20 && class( cv.lars.tmp ) == "try-error" ) {
    cv.lars.tmp <- try( cv.lars( df.tmp, output, K=cv.k, type="lasso", plot.it=plot ), silent=quiet )
    tries <- tries + 1
  } }
  if ( class( cv.lars.tmp ) == "try-error" ) return( numeric() )

  ## ... choose min CV (either absolute min if cv.choose="min" or min+1se if cv.choose="min+1se")
  min.i <- which.min( cv.lars.tmp$cv )
  min.err <- cv.lars.tmp$cv.error[ min.i ]
  ##if ( cv.min ) thresh.cv <- min( cv.lars.tmp$cv ) else
  thresh.cv <- min( cv.lars.tmp$cv ) + min.err

  if ( cv.choose[ 1 ] == "min+1se" ) best.s <- min( which( cv.lars.tmp$cv <= thresh.cv ) )
  else if ( cv.choose[ 1 ] == "min" ) best.s <- which.min( cv.lars.tmp$cv )
  if ( plot ) lines( rep( cv.lars.tmp$fraction[ best.s ], 2 ), range( cv.lars.tmp$cv ), col=2, lty=2, lwd=3 )
  if ( ! quiet ) cat( "cv: ", cv.choose, min.i, min.err, best.s, cv.lars.tmp$cv[ best.s ],
                     cv.lars.tmp$fraction[ best.s ], "\n" )
  coeffs <- coef.lars( lars.tmp, s=cv.lars.tmp$fraction[ best.s ], mode="fraction" )
  cv.err <- cv.lars.tmp$cv[ best.s ]
  coeffs <- coeffs[abs( coeffs ) >= cutoff ]
  if ( scale.coeffs ) {
   # plug coeffs from LARS back into LM, coeffs will be larger. Not necessarily kosher according to some people.
   if(length(coeffs) > 0){
        coeffs.s <- list();
            for (i in 1:length(coeffs)){
				if (length(coeffs) == 1 ){ins <- predictor.mat[names(coeffs),]
                } else {ins <- t(predictor.mat[names(coeffs),])}
				coeffs.s <- coef(lm(output ~ ins))
				coeffs.s <-  coeffs.s[2:length(coeffs.s)]
				names(coeffs.s) <- names(coeffs)
				#coeffs.debug <<- coeffs.s
            }
        coeffs <- unlist(coeffs.s)
    }
 }

#

#cat(del.ts) 
# df.tmp is a matrix of predictors
#cat(colnames(df.tmp))
#cat(rownames(predictor.mat))
# output is a vector of the biclust mean profile
# e.star.of.t <- predictor.mat.hidden * exp( unlist(t(apply(as.matrix(output),2,FUN=function(x)x-x[1]))  / tau)
# y <- output
# z.betas <- 
# <- unlist(t(apply(as.matrix(output),2,FUN=function(x)x-x[1])) 
# to.be.optimized <- function(tau, y, z.betas, e.star.of.t) {
#	return((y - z.betas - e.star.of.t)^2 )
# }
# optimize(to.be.optimized, tau=tau, y=output, interval=c(-5,5), )
 
 
 if ( ! quiet ) { cat( "NONZERO COEFFS:\n" ); print( coeffs ) }
 return(coeffs)


} # end of inferelate function


###########################################
# predictelator
###########################################
predictelator <- function(coeffs, clusterStack, oldratios, newratios, col.map, env.map=NULL, tfgroups.file, tau=10){
	# read in the coefficients named list
	# expect list to be of format <clusternumber>.predictorA.predictorB.func
	# !!!! VERY IMPORTANT: IF FORMAT OF COEFFS OUTPUT CHANGES, MUST CHANGE THIS FUNCTION ACCORDINGLY !!!!
	set.seed(31337)
	newratios <- mean.variance.normalize(newratios)
   	oldratios <- mean.variance.normalize(oldratios)
   	
   	#pdf("inferelator_diag.pdf", width=10, height=6)
   	
	tmp3 <- build.profiles.from.weights(coeffs=coeffs, ratios=newratios, tfgroups.file=tfgroups.file)
	mean.old.profiles <- matrix(ncol=ncol(oldratios), nrow=0)
	for (i in 1:length(clusterStack)){
		mean.old.profiles <- rbind(mean.old.profiles,apply(oldratios[clusterStack[[i]]$rows, ], 2, mean))		
   		rownames(mean.old.profiles)[i] <- paste("bicluster"	, i, sep="")
   	}
	tmp3.steady <- tmp3
	rmsd.steady <- NULL
	for (i in 1:nrow(tmp3.steady)){
		biclustnum<-rownames(tmp3.steady)[i]
		rmsd.s <- rmsd(mean.old.profiles[biclustnum,], tmp3.steady[biclustnum,])
		rmsd.steady <- c(rmsd.steady, rmsd.s)
	}
	rmsd.steady <<- rmsd.steady
	if(! is.null( col.map )){ par(mfcol=c(1,2)) }
	hist(rmsd.steady, main="RMSD's (steady-state)")
	
	#tmp3.steady <<- tmp3.steady
	mean.old.profiles <<- mean.old.profiles
	
	if ( ! is.null( col.map ) ){
		# if have col.map (i.e. we have time-series data),  put in the time component
		conds <- colnames( newratios )
		prevs <- col.map[ conds, "prevCol" ]
		
		nexts <- list()
		for (i in 1:length(conds)){nexts[i]= conds[i+1]}
		nexts[length(conds)] = conds[length(conds)]
		nexts <- unlist(nexts)
		
		del.ts <- col.map[ conds, "delta.t" ]
		del.ts[ del.ts < 1 ] <- 1
		out.tmp <- tmp3[,conds ] - tau * ( tmp3[ ,nexts ] - tmp3[ ,conds ]) / del.ts
		out.tmp[ out.tmp > 3.0 ] <- 3.0
		out.tmp[ out.tmp < -3.0 ] <- -3.0
		#out.tmp <<- out.tmp
		tmp3 <- out.tmp
		#print(tau * (tmp3[ good.i ] - tmp3[ prevs ]) / del.ts)
		#par(mfrow=c(3,1), ask=T)
	
		# calculate root-mean-variance-squared to measure how good our fit is
		rmsd.timeseries <- NULL
		for (i in 1:nrow(tmp3)){
			biclustnum<-rownames(tmp3)[i]
			rmsd.s <- rmsd(mean.old.profiles[biclustnum,], tmp3[biclustnum,])
			rmsd.timeseries <- c(rmsd.timeseries, rmsd.s)
			 #The section below plots out each bicluster mean profile INDIVIDUALLY, compared to steady-state predictions and timeseries
			#plot(mean.old.profiles[biclustnum,], type="l", main=paste(biclustnum, "mean profile"), ylim=c(-3,3), col="black")	
			#plot(tmp3.steady[biclustnum,], type="l", main=paste("steady-state, RMSD =", rmsd(mean.old.profiles[biclustnum,], tmp3.steady[biclustnum,])), ylim=c(-3,3), col="blue")
			#plot(tmp3[biclustnum,], type="l", main=paste("time-series, RMSD =", rmsd(mean.old.profiles[biclustnum,], tmp3[biclustnum,])), ylim=c(-3,3), col="red")
		}
		rmsd.timeseries <<- rmsd.timeseries
		hist(rmsd.timeseries, main="RMSD's (timeseries)")	
	}		
	return(tmp3)
	
} # end of predictelator function

###########################################
# build.profiles.from.weights
###########################################
# !!!! VERY IMPORTANT: IF FORMAT OF COEFFS OUTPUT CHANGES, MUST CHANGE THIS FUNCTION ACCORDINGLY !!!!
build.profiles.from.weights <- function(coeffs, ratios, tfgroups.file=NULL){
# this function takes a nested list of coeffs from LARS and construct new gene expression profiles
# if a list of tfgroups are provided, builds profiles for them too
	if (! is.null(tfgroups.file) ){
	   	tfgroups <- list()
	   	# retrieve the TFGROUPs from the tfgroup.noa file
	   	# and calculate their profiles for the old conditions
	   	tfgroup.ratios <- matrix(nrow=0, ncol=ncol(ratios))
	   	file.in <- as.matrix(read.csv2(tfgroups.file))
	    for ( k in 1:length(file.in)){
	    	rowname <- strsplit(file.in[k], " = ")[[1]][1]
	    	tfgroups <- c(tfgroups, strsplit(sub("\\)","", sub("\\(", "",as.character(strsplit(file.in[k], " = ")[[1]][2]))), "::"))
	    	names(tfgroups)[k] <- rowname
	    	tmp <- apply(as.matrix(ratios[tfgroups[[k]], ]), 2, mean)
	    	tfgroup.ratios <- rbind(tfgroup.ratios, tmp)
	    	rownames(tfgroup.ratios)[k] <- rowname    
	    }
	    # tfgroup.ratios is a matrix of the ratios for the TFGROUPs
	    ratios <- rbind(ratios, tfgroup.ratios)
	} 
 	tmp2 <- 0
	tmp3 <- matrix(nrow=0, ncol=ncol(ratios))
	for ( i in 1:length(coeffs) ) {
	cat(paste("Working on biclust", i, "\n"))	
		if( length(coeffs[[i]]) > 1 ) {	
			for ( j in 1:length(coeffs[[i]])){
				nodes <- unlist(strsplit(names(coeffs[[i]][j]), "\\."))
				weight <- coeffs[[i]][j]
				 	print(nodes)
				# there is a single predictor, multiply the predictor profile by the weight
				if (length(nodes) == 2) {
					#print(paste("Single predictor", nodes[2], "regulates", nodes[1], "with a weight of", weight) )
			  		tmp <- ratios[nodes[2], ] * weight	   	
			   	} else if (length(nodes) == 4){			    	
			   		# there are 2 predictors, get their weighted ratios as well as the combining function (min or max)
					#print(paste(nodes[2], nodes[4], "with", nodes[3], "regulates", nodes[1], "with a weight of", weight))
					func <- nodes[4]
					tmp <- apply(rbind(ratios[nodes[2], ], ratios[nodes[3], ]), MARGIN=2, FUN=func) * weight
					}
					tmp2 <- tmp2 + tmp
				}
				tmp3 <- rbind(tmp3, tmp2)	
				tmp2 <- 0
				#name the results
				rownames(tmp3)[dim(tmp3)[1]] <- paste("bicluster",i,sep="")
				colnames(tmp3) <- colnames(ratios)
		}
	}
	return(tmp3)
} # end of build.profiles.from.weights function



###########################################
# mean.variance.normalize
###########################################
mean.variance.normalize <- function(input.matrix){
	means <- apply( input.matrix, 1, function(x) mean(x, na.rm=T))
	sds <- apply( input.matrix, 1, function(x) sd(x, na.rm=T) )
	input.matrix <- apply( input.matrix, 2, "-", means )
	input.matrix <- apply( input.matrix, 2, "/", sds )
	return(input.matrix)
} # end of mean.variance.normalize function


###########################################
# rmsd
###########################################
rmsd <- function(predicted, observed){
# calculates the root-mean-square deviance between two vectors of numbers
	out<- sqrt(mean(abs(predicted - observed)^2))
	return(out)
} # end of rmsd function


###########################################
# make.network.files
###########################################
make.network.files <- function(inf.result, clusterStack, sif.filename){
# this function takes an Inferelator result and a cM clusterStack
# and outputs a network file and the associated edge and node
# attribute files for visualizing the network in Cytoscape

# IMPORTANT: this function depends on the output returned by inferelator()
# change the below accordingly if the output there changes !

	out<- unlist(inf.result)
	# write the weights of the TFGROUPs to an edge-attribute file 
	# at the same time create the actual network .sif
	write("weight (java.lang.Double)", "weights.eda")
	
	gatecount = 1
	for (j in 1:length(out)){
	nodes<- strsplit(names(out)[j], "\\.")[[1]]
	#print(out[j])
	# make a .sif file for the actual network
	# .sif format : [node]<tab>[relationship]<tab>[node]

	#print(gatecount)
	# TODO : load this into a data frame and write.table() or something similar 
	# if there is only 1 TFGROUP, create two nodes
	if (length(nodes) == 2) {
	    write(paste(nodes[2], "activates", nodes[1]), sif.filename, append=T)
	    
	    write(paste(nodes[2], "(activates)", nodes[1],  "=", out[j]), "weights.eda", append=T)
	    
		} else if (length(nodes) == 4){
		# there are 2 TFGROUP's, create an AND gate (a Y-shaped segment with 4 nodes and 3 edges)
		
		write(paste(nodes[2], "combines", paste("AND-", gatecount, sep="")), sif.filename, append=T)
		write(paste(nodes[3], "combines", paste("AND-", gatecount, sep="")), sif.filename, append=T)
		write(paste(paste("AND-", gatecount, sep=""), "activates", nodes[1]), sif.filename, append=T)
		
		write(paste(paste("AND-", gatecount, sep=""),"(activates)", nodes[1],  "=", out[j]), "weights.eda", append=T)
		
		write(paste(paste("AND-",gatecount,sep=""), "=", "(logicGate)"), "types.noa", append=T)
		gatecount = gatecount + 1
		
		# there is only a cluster (no significant coeffs from inferelator), just write the cluster as a node to file
		} else if (length(nodes) == 1){
		write(nodes[1], sif.filename, append=T)
		}
	
	}
	# get attributes of the biclusts from clusterStack (e.g. no. of genes, conds, p-vals)
	# and write to appropriate node-attribute files
	write("clusterGenes", "clusterGenes.noa")
	write("clusterConditions", "clusterConditions.noa")	
	write("clusterGeneCount", "clusterGeneCount.noa")
	write("clusterConditionCount", "clusterConditionCount.noa")
	write("clusterMotifPValues", "clusterMotifPValues.noa")
	write("clusterResiduals", "clusterResiduals.noa")
	
	for (i in 1:length(clusterStack)){
	 	write( paste(paste(i, " = ", "(", sep="") ,paste(clusterStack[[i]]$rows, collapse="::"), ")", sep="") , "clusterGenes.noa", append=T)	
	 	write( paste(paste(i, " = ", "(", sep="") ,paste(clusterStack[[i]]$cols, collapse="::"), ")", sep="") , "clusterConditions.noa", append=T)
	 	write( paste(i, " = ", clusterStack[[i]]$nrows , sep="") , "clusterGeneCount.noa", append=T)
	 	write( paste(i, " = ", clusterStack[[i]]$ncols, sep="") , "clusterConditionCount.noa", append=T)
	 	write( paste(i, " = ", clusterStack[[i]]$e.val, sep="") , "clusterMotifPValues.noa", append=T)
	 	write( paste(i, " = ", clusterStack[[i]]$resid, sep="") , "clusterResiduals.noa", append=T)
	 	write( paste(i, " = (cluster)",sep=""), "types.noa", append=T )
	
	}
	
} # end of make.network.files function