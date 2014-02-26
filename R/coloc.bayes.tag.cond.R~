#'Colocation for two traits with a common control, conditioning on SNPs
#'Merges SNPs with high r2 into tags prior to analysis
#'Generates bayes factors for each plausible one tag model
#'
#'
#'@title bayesian colocalisation; two traits; with tagging; conditional
#'@export
#'@param df1 A dataframe, containing response and potential explanatory variables for the dataset.
#'@param snps The SNPs to consider as potential explanatory variables
#'@param response The name of the response variable in \code{df1}
#'@param priors A list of priors over the hypotheses 
#'@param pp.thr posterior probability threshold used to trim SNP list.  Only SNPs with a marginal posterior probability of inclusion greater than this with one or other trait will be included in the full BMA analysis
#'@param r2.trim If a pairs of SNPs has r2 greater than \code{r2.trim}, they are put in the same tag
#'@param quiet suppress messages about how the model spaced is trimmed for BMA
#'@param cond1 A list of the SNPs we are to condition upon for trait1
#'@param cond2 A list of the SNPs we are to condition upon for trait2
#'@return a list of posterior probabilities that each tag is causitive to both traits, the tag names, and the corresponding SNPs
#'@author Mary Fortune
coloc.bayes.tag.cond <- function(df1,snps=setdiff(colnames(df1),response),response="Y",priors=list(c(1,1,1,1,1)),r2.trim=0.99,pp.thr=0.005,quiet=TRUE,cond1=NULL,cond2=NULL) {
    cond<-union(cond1,cond2)
    #we consider all models which contain at most 1 snps for each trait, using tagging
	#we condition on the SNPs present in cond
	if (length(cond)==0){
		return(coloc.bayes.tag(df1,snps=snps,response=response,priors=priors,r2.trim=r2.trim,pp.thr=pp.thr,quiet=quiet))
	}
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
	orig.tags<-tags
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
    while(any(is.na(coefficients(lm1)[-c(1,ncol(binX)+1)]))) {
        z.drop <- colnames(binX)[which(is.na(coefficients(lm1))[-c(1,ncol(binX)+1)])]
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

    #seperate the tags and the conditioning SNPs
	tagdrop<-c()
	for (tt in 1:n.clean){
		firsttag<-tags[tt]
		fulltag<-names(tagkey[which(tagkey==firsttag)])
		numcond<-length(which(fulltag %in% cond))
		orig.place<-which(orig.tags==firsttag)
		if (numcond>0){
			#there is at least one conditional SNP, which we must remove
			if (numcond==length(fulltag)){
				#the tag is entirely from cond - remove the tag
				tagdrop<-c(tagdrop,tt)
				tagsize[orig.place]=tagsize[orig.place]-numcond
			}
			else{
				#remove numcond SNPs from the tag and rename if necessary
				tagsize[orig.place]=tagsize[orig.place]-numcond
				if (firsttag %in% cond){
					tags[tt]<-setdiff(fulltag,cond[which(cond %in% fulltag)])[1]		
					orig.tags[orig.place]<-tags[tt]
					tagkey[which(tagkey==firsttag)]<-tags[tt]
				}
			}
		}
	}
	if (length(tagdrop)>0){
		tags<-tags[-tagdrop]
	}

	n.clean<-length(tags)
    #remove tags with low posterior probabilities in the individual models.
	#extract just those samples relating to trait1
	df.trait1<-df1[which(df1[,1]!=2),c("Y",tags,cond1)]
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
		n.clean <- length(tags)
	}
	if ( intrank < ncol(intX1)){
		#we need to remove one of cond1
		#we must remove at k columns
		k<-ncol(intX1)-intrank
		drop<-c()
		for(ii in 1:length(cond1)) { 
			if (length(drop)<k){ 
				xsub <- intX1[,setdiff(colnames(intX1),cond1[ii])]
				droprank<- qr(xsub)$rank
				if (droprank==intrank){
					drop<-c(drop,cond1[ii])
					intX1<-xsub
				}
			}
		}
		cond1<-setdiff(cond1,drop)
		df.trait1<-df.trait1[,- which(colnames(df.trait1) %in% drop)]
		cond<-union(cond1,cond2)
	}
	#run the model
	modelsep<-cbind(diag(n.clean),matrix(1,ncol=length(cond1),nrow=n.clean))
	mod.trait1<-glib(x=df.trait1[,-1], y=df.trait1[,1], error="binomial", link="logit",models=modelsep)
	pp.trait1<-mod.trait1$bf$postprob[,2]
	whichtags1<-tags[which(pp.trait1>pp.thr)]
	#extract just those samples relating to trait2
	df.trait2<-df1[which(df1[,1]!=1),c("Y",tags,cond2)]
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
		n.clean <- length(tags)
	}
	if ( intrank < ncol(intX2)){
		#we need to remove one of cond2
		#we must remove at k columns
		k<-ncol(intX2)-intrank
		drop<-c()
		for(ii in 1:n.clean) { 
			if (length(drop)<k){ 
				xsub <- intX2[,setdiff(colnames(intX2),cond2[ii])]
				droprank<- qr(xsub)$rank
				if (droprank==intrank){
					drop<-c(drop,cond2[ii])
					intX2<-xsub
				}
			}
		}
		cond2<-setdiff(cond2,drop)
		df.trait2<-df.trait2[,- which(colnames(df.trait2) %in% drop)]
		cond<-union(cond1,cond2)
	}
	#run the model
	modelsep<-cbind(diag(n.clean),matrix(1,ncol=length(cond2),nrow=n.clean))
	mod.trait2<-glib(x=df.trait2[,-1], y=df.trait2[,1]/2, error="binomial", link="logit",models=modelsep)
	pp.trait2<-mod.trait2$bf$postprob[,2]
	whichtags2<-tags[which(pp.trait2>pp.thr)]
	whichtags<-union(whichtags1,whichtags2)
	tags<-intersect(whichtags,tags)
    n.clean <- length(tags)


	cat("We consider ",n.clean, " tags in the final analysis, and condition upon ",length(cond), " SNPs \n")
	if (n.clean ==0){
		print("for prior:  0.99979 1e-04 1e-04 1.00021e-08 1e-05")
		print(" we have posterior:  0 0 0 0 0") 
		print("---")
		return (1)	
	}
    

    #covert to a binomial model so we can run glib
	#note that the variables corresponding to the conditioning SNPs are included here
    f1 <- as.formula(paste("Y ~ 1 | ", paste(c(tags,cond),collapse="+")))
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
	n.clean <- length(tags)
	tagsize<-tagsize[unlist(lapply(tags,function(x){which(orig.tags %in% x)}))]
	#which models are we testing?
	#the matrix of tag models only
    tagmodels<-makebinmod(n.clean,1)
    category<-apply(tagmodels,1,whichcat)
	#add the conditionals to this
	cond1row<-rep(0,length(cond))
	cond1row[which(cond %in% cond1)]=1
	cond1models<-t(matrix(rep(cond1row,times=nrow(tagmodels)),nrow=length(cond)))
	cond2row<-rep(0,length(cond))
	cond2row[which(cond %in% cond2)]=1
	cond2models<-t(matrix(rep(cond2row,times=nrow(tagmodels)),nrow=length(cond)))
	models<-cbind(tagmodels[,1:n.clean],cond1models,tagmodels[,(n.clean+1):(2*n.clean)],cond2models,rep(1,times=nrow(tagmodels)))
    #run glib
    mods1 <- glib(x=binX, y=binY, error="binomial", link="logit",models=models)
    #-2 log the bf with a flat prior
    logB10=0.5*mods1$bf$twologB10[,2]
    logbf <- numeric(5)
    #find which SNPs are present in each model
    if (n.clean>1){
        n1 <- as.vector( models[,1:n.clean] %*% tagsize )
        n2 <- as.vector( models[,(n.clean+length(cond)+1):(2*n.clean+length(cond))] %*% tagsize )
    }else{
        n1<-rep(tagsize,length(category))
        n2<-rep(tagsize,length(category))
    }
    logbf[1]<-logsum(logB10[which(category==0)])
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
    pp<-(exp( tmp - logsum(tmp) ))
    snplist<-vector("list",length(pp))
    for (ii in 1:length(pp)){
        snplist[[ii]]<-setdiff(names(which(tagkey==tags[ii])),cond)
    }
    return(list("postprob"=pp,"tags"=tags,"snps"=snplist,"models"=models,"bf"=logB10))
}
