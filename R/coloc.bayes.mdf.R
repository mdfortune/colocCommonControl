
#'Colocation for two traits with a common control
#'Generates bayes factors for each plausible one SNP model
#' 
#'
#'@title bayesian colocalisation; two traits
#'@export
#'@param df1 A dataframe, containing response and potential explanatory variables for the dataset.
#'@param snps The SNPs to consider as potential explanatory variables
#'@param response The name of the response variable in \code{df1}
#'@param priors A list of priors over the hypotheses 
#'@param pp.thr posterior probability threshold used to trim SNP list.  Only SNPs with a marginal posterior probability of inclusion greater than this with one or other trait will be included in the full BMA analysis
#'@param r2.trim for pairs SNPs with r2 greater than \code{r2.trim}, only one SNP will be retained.  This avoids numerical instability problems caused by including two highly correlated SNPs in the model.
#'@param quiet suppress messages about how the model spaced is trimmed for BMA
#'@return a list of posterior probabilities that each SNP is causitive to both traits, and the corresponding SNPs
#'@author Mary Fortune
coloc.bayes <- function(df1,snps=setdiff(colnames(df1),response),response="Y",priors=list(c(1,1,1,1,1)),r2.trim=0.99,pp.thr=0.005,quiet=TRUE) {
    #we consider all models which contain at most 1 snp for each trait
    snps <- unique(snps)
    n.orig <- length(snps)
    if(n.orig<2)
        return(1)
    prep <- prepare.df(df1, snps, r2.trim=r2.trim, dataset=1, quiet=quiet)
    df1 <- prep$df
    snps <- prep$snps
    
    if(!quiet)
        cat("Dropped",n.orig - length(snps),"of",n.orig,"SNPs due to LD: r2 >",r2.trim,"\n",length(snps),"SNPs remain.\n")
    
    ## remove any completely predictive SNPs
    f1 <- as.formula(paste("Y ~", paste(snps,collapse="+")))
    capture.output(lm1 <- multinom(f1, data=df1,maxit=1000))
    while(any(is.na(coefficients(lm1)))) {
        drop <- which(is.na(coefficients(lm1))[-1])
        if(!quiet)
            cat("Dropping",length(drop),"inestimable SNPs (most likely due to colinearity):\n",drop,"\n")
        snps <- snps[-drop]
        f1 <- as.formula(paste("Y ~", paste(snps,collapse="+")))
        capture.output(lm1 <- multinom(f1, data=df1,maxit=1000))
    }
    n.clean <- length(snps)
    #remove SNPs with low posterior probabilities in the individual models.
	#extract just those samples relating to trait1
	df.trait1<-df1[which(df1[,1]!=2),c("Y",snps)]
	#remove columns if necessary to ensure a full rank X
	intX1<-cbind(Intercept=1,df.trait1[,-1])
	intrank<-qr(intX1)$rank
	if ( intrank < ncol(intX1)){
		#we must remove at k columns
		k<-ncol(intX1)-intrank
		drop<-c()
		for(ii in 1:n.clean) { 
			if (length(drop)<k){ 
				xsub <- intX1[,setdiff(colnames(intX1),snps[ii])]
				droprank<- qr(xsub)$rank
				if (droprank==intrank){
					drop<-c(drop,snps[ii])
					intX1<-xsub
				}
			}
		}
		snps<-setdiff(snps,drop)
		df.trait1<-df.trait1[,- which(colnames(df.trait1) %in% drop)]
		df1<-df1[,-which(colnames(df1) %in% drop)]
		n.clean <- length(snps)
	}
	#run the model
	f1 <- as.formula(paste("Y ~", paste(snps,collapse="+")))
	modelsep<-diag(n.clean)
	mod.trait1<-glib(x=df.trait1[,-1], y=df.trait1[,1], error="binomial", link="logit",models=modelsep)
	pp.trait1<-mod.trait1$bf$postprob[,2]
	whichsnps1<-snps[which(pp.trait1>pp.thr)]
	#extract just those samples relating to trait2
	df.trait2<-df1[which(df1[,1]!=1),c("Y",snps)]
	#remove columns if necessary to ensure a full rank X
	intX2<-cbind(Intercept=1,df.trait2[,-1])
	intrank<-qr(intX2)$rank
	if ( intrank < ncol(intX2)){
		#we must remove at k columns
		k<-ncol(intX2)-intrank
		drop<-c()
		for(ii in 1:n.clean) { 
			if (length(drop)<k){ 
				xsub <- intX2[,setdiff(colnames(intX2),snps[ii])]
				droprank<- qr(xsub)$rank
				if (droprank==intrank){
					drop<-c(drop,snps[ii])
					intX2<-xsub
				}
			}
		}
		snps<-setdiff(snps,drop)
		df.trait2<-df.trait2[,- which(colnames(df.trait2) %in% drop)]
		df2<-df2[,-which(colnames(df2) %in% drop)]
		n.clean <- length(snps)
	}
	#run the model
	f2 <- as.formula(paste("Y ~", paste(snps,collapse="+")))
	modelsep<-diag(n.clean)
	mod.trait2<-glib(x=df.trait2[,-1], y=df.trait2[,1]/2, error="binomial", link="logit",models=modelsep)
	pp.trait2<-mod.trait2$bf$postprob[,2]
	whichsnps2<-snps[which(pp.trait2>pp.thr)]
	whichsnps<-union(whichsnps1,whichsnps2)
	snps<-intersect(whichsnps,snps)
    n.clean <- length(snps)
	cat("We consider ",n.clean, " snps in the final analysis. \n")
    #covert to a binomial model so we can run glib
    f1 <- as.formula(paste("Y ~ 1 | ", paste(snps,collapse="+")))
    binmod<-mlogit2logit(f1,data=df1,choices=0:2,base.choice = 1)$data
    index<-grep("z",colnames(binmod))
    binX<-binmod[,index]
    #remove z_1 (we do not want it in our model matrix)
    binX<-binX[,-1]
    #extract the new reponse
    binY<-binmod[,"Y.star"]
	## remove SNPs to ensure that binX is full rank when given an intercept
	intbinX<-cbind(Intercept=1,binX)
	intrank<-qr(intbinX)$rank
	if ( intrank < ncol(intbinX)){
		#we must remove at k columns
		k<-ncol(intbinX)-intrank
		drop<-c()
		for(ii in 1:n.clean) {
			if (intrank < ncol(intbinX)){ 
				xsub <- intbinX[,setdiff(colnames(intbinX),paste0(c("z_1.","z_2."),snps[ii]))]
				droprank<- qr(xsub)$rank
				if (droprank>intrank-2){
					drop<-c(drop,snps[ii])
					intbinX<-xsub
					intrank<-qr(intbinX)$rank
				}
			}
		}
		bindrop<-c(paste0("z_1.",drop),paste0("z_2.",drop))
		snps<-setdiff(snps,drop)
		binX<-binX[,- which(colnames(binX) %in% bindrop)]
		n.clean <- length(snps)
	}
	#which models are we testing?
    models<-makebinmod(n.clean,1)
    category<-apply(models,1,whichcat)
    #run glib
    mods1 <- glib(x=binX, y=binY, error="binomial", link="logit",models=models)
    #-2 log the bf with a flat prior
    twologB10=mods1$bf$twologB10[,2]
    logB10=0.5*twologB10
    logbf<-rep(0,5)
    logbf[1]<-wlogsum(0.5*twologB10[which(category==0)])
    logbf[2]<-wlogsum(0.5*twologB10[which(category==1)])
    logbf[3]<-wlogsum(0.5*twologB10[which(category==2)])
    logbf[4]<-wlogsum(0.5*twologB10[which(category==3)])
    logbf[5]<-wlogsum(0.5*twologB10[which(category==4)])
    postprobs<-vector("list",length(priors))
    for (ii in 1:length(priors)){
        prior<-priors[[ii]]
        postlogbf<-logbf+log(prior)
        postprob<-(exp(postlogbf))
        postprob<-postprob/sum(postprob)
        cat("for prior: ", prior/sum(prior), "\n")
        cat("we have posterior: ",postprob, "\n")
        cat("--- \n")
        postprobs[[ii]]=postprob
    }
    tmp <- 0.5*twologB10[which(category==4)]
    pp<-(exp( tmp - wlogsum(tmp) ))
    return(list("postprob"=pp,"snps"=snps,"models"=models,"bf"=logB10))
}

#'Colocation for two traits with a common control
#'Merges SNPs with high r2 into tags prior to analysis
#'Generates bayes factors for each plausible one tag model
#'
#'
#'@title bayesian colocalisation; two traits; with tagging
#'@export
#'@param df1 A dataframe, containing response and potential explanatory variables for the dataset.
#'@param snps The SNPs to consider as potential explanatory variables
#'@param response The name of the response variable in \code{df1}
#'@param priors A list of priors over the hypotheses 
#'@param pp.thr posterior probability threshold used to trim SNP list.  Only SNPs with a marginal posterior probability of inclusion greater than this with one or other trait will be included in the full BMA analysis
#'@param r2.trim If a pairs of SNPs has r2 greater than \code{r2.trim}, they are put in the same tag
#'@param quiet suppress messages about how the model spaced is trimmed for BMA
#'@return a list of posterior probabilities that each tag is causitive to both traits, the tag names, and the corresponding SNPs
#'@author Mary Fortune
coloc.bayes.tag <- function(df1,snps=setdiff(colnames(df1),response),response="Y",priors=list(c(1,1,1,1,1)),r2.trim=0.99,pp.thr=0.005,quiet=TRUE) {
    #we consider all models which contain at most 1 snps for each trait, using tagging
    snps <- unique(snps)
    n.orig <- length(snps)
    
    #generate tags using r2.trim
    invisible(capture.output(X<-new("SnpMatrix",as.matrix(df1[,-1]))))
    tagkey <- tag(X,tag.threshold=r2.trim)
    tags<-unique(tagkey)
    n.clean=length(tags)
    tagsize<-rep(0,n.clean)
    for (ii in 1:n.clean){
        tagsize[ii]<-length(which(tagkey==tags[ii]))
    }
    
 	#make binomial model
    f1 <- as.formula(paste("Y ~ 1 | ", paste(tags,collapse="+")))
    x1 <- df1[,tags]
    n.clean <- length(tags)
    binmod<-mlogit2logit(f1,data=df1,choices=0:2,base.choice = 1)$data
    index<-grep("z",colnames(binmod))
    binX<-binmod[,index]
    #remove z_1
    binX<-binX[,-1]
	binY<-as.matrix(binmod[,"Y.star"])
	colnames(binY)="Y"


   ## remove tags to ensure that binX is full rank when given an intercept
	intbinX<-cbind(Intercept=1,binX)
	intrank<-qr(intbinX)$rank
	if ( intrank < ncol(intbinX)){
		#we must remove at k columns
		k<-ncol(intbinX)-intrank
		drop<-c()
		for(ii in 1:n.clean) {
			if (intrank < ncol(intbinX)){ 
				xsub <- intbinX[,setdiff(colnames(intbinX),paste0(c("z_1.","z_2."),tags[ii]))]
				droprank<- qr(xsub)$rank
				if (droprank>intrank-2){
					drop<-c(drop,tags[ii])
					intbinX<-xsub
					intrank<-qr(intbinX)$rank
				}
			}
		}
		bindrop<-c(paste0("z_1.",drop),paste0("z_2.",drop))
		tags<-setdiff(tags,drop)
		binX<-binX[,- which(colnames(binX) %in% bindrop)]
		n.clean <- length(tags)
	}

	# remove any completely predictive tags
	 f <- as.formula(paste("Y ~", paste(colnames(binX),collapse="+")))
    capture.output(lm1 <- glm(f,data=as.data.frame(cbind(binY,binX)),family="binomial"))
    while(any(is.na(coefficients(lm1)))) {
        z.drop <- which(is.na(coefficients(lm1))[-c(1,ncol(binX)+1)])
		drop<-unique(gsub("z_2.","",gsub("z_1.","",z.drop)))
        if(!quiet)
            cat("Dropping",length(drop),"inestimable tags (most likely due to colinearity):\n",drop,"\n")
        bindrop<-c(paste0("z_1.",drop),paste0("z_2.",drop))
		tags<-setdiff(tags,drop)
		binX<-binX[,- which(colnames(binX) %in% bindrop)]
        f <- as.formula(paste("Y ~", paste(colnames(binX),collapse="+")))
        capture.output(lm1 <- glm(f,data=as.data.frame(cbind(binY,binX)),family="binomial"))
    }

    n.clean<-length(tags)
	binY<-binmod[,"Y.star"]
    
    #remove tags with low posterior probabilities in the individual models.
	#extract just those samples relating to trait1
	df.trait1<-df1[which(df1[,1]!=2),c("Y",tags)]
	#remove columns if necessary to ensure a full rank X
	intX1<-cbind(Intercept=1,df.trait1[,-1])
	intrank<-qr(intX1)$rank
	if ( intrank < ncol(intX1)){
		#we must remove at k columns
		k<-ncol(intX1)-intrank
		drop<-c()
		for(ii in 1:n.clean) { 
			if (length(drop)<k){ 
				xsub <- intX1[,setdiff(colnames(intX1),tags[ii])]
				droprank<- qr(xsub)$rank
				if (droprank==intrank){
					drop<-c(drop,tags[ii])
					intX1<-xsub
				}
			}
		}
		tags<-setdiff(tags,drop)
		df.trait1<-df.trait1[,- which(colnames(df.trait1) %in% drop)]
		df1<-df1[,-which(colnames(df1) %in% drop)]
		n.clean <- length(tags)
	}
	#run the model
	f1 <- as.formula(paste("Y ~", paste(tags,collapse="+")))
	modelsep<-diag(n.clean)
	mod.trait1<-glib(x=df.trait1[,-1], y=df.trait1[,1], error="binomial", link="logit",models=modelsep)
	pp.trait1<-mod.trait1$bf$postprob[,2]
	whichtags1<-tags[which(pp.trait1>pp.thr)]
	#extract just those samples relating to trait2
	df.trait2<-df1[which(df1[,1]!=1),c("Y",tags)]
	#remove columns if necessary to ensure a full rank X
	intX2<-cbind(Intercept=1,df.trait2[,-1])
	intrank<-qr(intX2)$rank
	if ( intrank < ncol(intX2)){
		#we must remove at k columns
		k<-ncol(intX2)-intrank
		drop<-c()
		for(ii in 1:n.clean) { 
			if (length(drop)<k){ 
				xsub <- intX2[,setdiff(colnames(intX2),tags[ii])]
				droprank<- qr(xsub)$rank
				if (droprank==intrank){
					drop<-c(drop,tags[ii])
					intX2<-xsub
				}
			}
		}
		tags<-setdiff(tags,drop)
		df.trait2<-df.trait2[,- which(colnames(df.trait2) %in% drop)]
		df1<-df1[,-which(colnames(df1) %in% drop)]
		n.clean <- length(tags)
	}
	#run the model
	f2 <- as.formula(paste("Y ~", paste(tags,collapse="+")))
	modelsep<-diag(n.clean)
	mod.trait2<-glib(x=df.trait2[,-1], y=df.trait2[,1]/2, error="binomial", link="logit",models=modelsep)
	pp.trait2<-mod.trait2$bf$postprob[,2]
	whichtags2<-tags[which(pp.trait2>pp.thr)]
	whichtags<-union(whichtags1,whichtags2)
	tags<-intersect(whichtags,tags)
    	n.clean <- length(tags)
	cat("We consider ",n.clean, " tags in the final analysis. \n")

    tagsize<-tagsize[which(unique(tagkey) %in% tags)]

    #covert to a binomial model so we can run glib
    f1 <- as.formula(paste("Y ~ 1 | ", paste(tags,collapse="+")))
    binmod<-mlogit2logit(f1,data=df1,choices=0:2,base.choice = 1)$data
    index<-grep("z",colnames(binmod))
    binX<-binmod[,index]
    #remove z_1 (we do not want it in our model matrix)
    binX<-binX[,-1]
    #extract the new reponse
    binY<-binmod[,"Y.star"]
	## remove tags to ensure that binX is full rank when given an intercept
	intbinX<-cbind(Intercept=1,binX)
	intrank<-qr(intbinX)$rank
	if ( intrank < ncol(intbinX)){
		#we must remove at k columns
		k<-ncol(intbinX)-intrank
		drop<-c()
		for(ii in 1:n.clean) {
			if (intrank < ncol(intbinX)){ 
				xsub <- intbinX[,setdiff(colnames(intbinX),paste0(c("z_1.","z_2."),tags[ii]))]
				droprank<- qr(xsub)$rank
				if (droprank>intrank-2){
					drop<-c(drop,tags[ii])
					intbinX<-xsub
					intrank<-qr(intbinX)$rank
				}
			}
		}
		bindrop<-c(paste0("z_1.",drop),paste0("z_2.",drop))
		tags<-setdiff(tags,drop)
		binX<-binX[,- which(colnames(binX) %in% bindrop)]
		n.clean <- length(tags)
	}
	#which models are we testing?
    models<-makebinmod(n.clean,1)
    category<-apply(models,1,whichcat)
    #run glib
    mods1 <- glib(x=binX, y=binY, error="binomial", link="logit",models=models)
    #-2 log the bf with a flat prior
    logB10=0.5*mods1$bf$twologB10[,2]
    logbf <- numeric(5)
    #find which SNPs are present in each model
    if (n.clean>1){
        n1 <- as.vector( models[,1:n.clean] %*% tagsize )
        n2 <- as.vector( models[,(n.clean+1):(2*n.clean)] %*% tagsize )
    }else{
        n1<-rep(tagsize,length(category))
        n2<-rep(tagsize,length(category))
    }
    logbf[1]<-wlogsum(logB10[which(category==0)])
    wh1 <- which(category==1)
    logbf[2]<-wlogsum(logB10[wh1], n1[wh1])
    wh2 <- which(category==2)
    logbf[3]<-wlogsum(logB10[wh2], n2[wh2])
    wh3 <- which(category==3)
    wh4 <- which(category==4)
    logbf[4]<-wlogsum(c(logB10[wh3],logB10[wh4]), c(n1[wh3] * n2[wh3],n1[wh4]*(n1[wh4]-1)))    
    logbf[5]<-wlogsum(logB10[wh4], n1[wh4]) ## assumes n1==n2 for cat 4, can't see why this wouldn't be true, but untested
    postprobs<-vector("list",length(priors))
    for (ii in 1:length(priors)){
        prior<-priors[[ii]]
        postlogbf<-logbf-max(logbf)+log(prior)
        postprob<-(exp(postlogbf))
        postprob<-postprob/sum(postprob)
        cat("for prior: ", prior/sum(prior), "\n")
        cat("we have posterior: ",postprob, "\n")
        cat("--- \n")
        postprobs[[ii]]=postprob
    }
    tmp <- logB10[which(category==4)]
    pp<-(exp( tmp - wlogsum(tmp) ))
    snplist<-vector("list",length(pp))
    for (ii in 1:length(pp)){
        snplist[[ii]]<-names(which(tagkey==tags[ii]))
    }
    return(list("postprob"=pp,"tags"=tags,"snps"=snplist,"models"=models,"bf"=logB10))
}

##' Internal function, makebinmod
##'
##' This function takes in a line from the model matrix
##' and computes which of the categories it corresponds to
##' @title makebinmod
##' @param p the number of SNPs
##' @param m the maximum number of causitive snps for each trait
##' @return a numeric matrix giving the models 
##' @author Mary Fortune
makebinmod<-function(p,m){
    #returns a model matrix for the binomial equivalent model
    #p=number of snps present
    #m=max number of snps in each model
    if (m>p) {m=p}
    snplist<-vector("list",2*m)
    for (ii in 1:(2*m)){snplist[[ii]]=0:p}
    whichsnps<-expand.grid(snplist)
    if (m >1){
        whichsnps<-unique(t(rbind(apply(whichsnps[,1:m],1,sort),apply(whichsnps[,(m+1):(2*m)],1,sort))))
    }
    whichsnps.t1<-as.matrix(whichsnps[,1:m])
    whichsnps.t2<-as.matrix(whichsnps[,(1+m):(2*m)])
    nummod<-nrow(whichsnps)
    models.t1<-matrix(0,nummod,p)
    models.t2<-matrix(0,nummod,p)
    for (ii in 1:nummod){
        for (jj in 1:m){
            if (whichsnps.t1[ii,jj]>0){
                models.t1[ii,whichsnps.t1[ii,jj]]=1
            }
            if (whichsnps.t2[ii,jj]>0){
                models.t2[ii,whichsnps.t2[ii,jj]]=1
            }
        }
    }
    models<-cbind(models.t1,models.t2,rep(1,nummod))
    if (m>1){
        models<-unique(models)
    }
    return(models)
}

##' Internal function, whichcat
##'
##' This function takes in a line from the model matrix
##' and computes which of the categories it corresponds to
##' @title whichcat
##' @param line numeric vector
##' @return a number corresponding to the category 
##' @author Mary Fortune
whichcat<-function(line){
    #puts the model given in line into one of the five categories
    #assumes m=1
    p<-(length(line)-1)/2
    t1<-line[1:p]
    t2<-line[(p+1):(2*p)]
    if (sum(t1)==0 & sum(t2)==0){
        return(0)
    }else if (sum(t2)==0){
        return(1)
    }else if (sum(t1)==0){
        return(2)
    }else if (sum(abs(t1-t2))==0){
        return(4)
    }else{
        return(3)
    }
}


##' Internal function, wlogsum
##'
##' This function calculates the log of the sum of the exponentiated
##' logs taking out the max, i.e. insuring that the sum is not Inf
##' This sum is weighted by some constants w
##' @title wlogsum
##' @param x numeric vector
##' @param w numeric vector
##' @return my.max + log(sum(exp(x - my.max )*w))
##' @author Chris Wallace
wlogsum <- function(x, w=NULL) {
    if (length(x)==1){
        if(is.null(w)) {
            return(x)
        }else{
            return(x*w)
        }
    }
    my.max <- max(x) ##take out the maximum value in log form
    if (my.max == -Inf){return(-Inf)}
    if(is.null(w)) {
        my.max + log(sum(exp(x - my.max )))
    } else {
        my.max + log(sum(exp(x - my.max )*w))
    }
}


##'Returns the r2 values between each pair of SNPs
##'@title find r2
##'@param X a SnpMatrix
##'@return a matrix of the r2 values
##'@author Chris Wallace
myr2 <- function(X) {
  r2 <- ld(X,
           depth=ncol(X)-1,
           symmetric=TRUE,
           stat="R.squared")
if(any(is.na(r2))) {
    r2.na <- as(is.na(r2),"matrix")
    use <- rowSums(r2.na)>0
## work around for r2=NA bug.
    r2.cor <- as(cor(as(X[,use,drop=FALSE],"numeric"), use="pairwise.complete.obs")^2,"Matrix")
    r2[ which(r2.na) ] <- r2.cor[ which(r2.na[use,use]) ]
}
  diag(r2) <- 1
return(r2)
}


##' Derive tag SNPs for a SnpMatrix object using heirarchical clustering
##'
##' Uses complete linkage and the \code{\link{hclust}} function to define clusters, then cuts the tree at 1-tag.threshold
##' @title tag
##' @param X a SnpMatrix
##' @param snps colnames of the SnpMatrix object to be used
##' @param tag.threshold threshold to cut tree, default=0.99
##' @param samples optional, subset of samples to use
##' @return character vector, names are \code{snps}, values are the tag for each SNP
##' @author Chris Wallace
##' @export
tag <- function(X,tag.threshold=0.99, snps=NULL, samples=NULL) {
if(!is.null(snps) || !is.null(samples))
    X <- X[samples,snps]
  r2 <- myr2(X)
   D <- as.dist(1-r2)
   hc <- hclust(D, method="complete")
   clusters <- cutree(hc, h=1-tag.threshold)
   snps.use <- names(clusters)[!duplicated(clusters)]
   r2.use <- r2[snps.use, colnames(X), drop=FALSE]
   tags <- rownames(r2.use)[apply(r2.use,2,which.max)]
   names(tags) <- colnames(r2.use)
return(tags)
}
