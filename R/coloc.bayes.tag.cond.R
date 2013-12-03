coloc.bayes.tag.cond <- function(df1,snps=setdiff(colnames(df1),response),response="Y",priors=list(c(1,1,1,1,1)),r2.trim=0.99,pp.thr=0.005,quiet=TRUE,cond=NULL) {
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

	#seperate the tags and the conditioning SNPs
	for (cc in 1:length(cond)){
		#which tag is conditional SNP cc in, and what is its index
		ccsnp<-cond[cc]
		cctag<-tagkey[which(names(tagkey)==ccsnp)]
		ccindex<-which(tags==cctag)
		if (length(ccindex)>0){
			if (cctag==ccsnp){
				if (tagsize[ccindex]==1){
					#remove this tag
					tags<-tags[-ccindex]
					tagsize<-tagsize[-ccindex]
					tagkey<-tagkey[-ccindex]
				}else{
					#index this tag by a different snp
					tagsize[ccindex]=tagsize[ccindex]-1
					tags[ccindex]<-setdiff(names(tagkey[which(tagkey==cctag)]),cctag)[1]
				}
			}else{
				#the tag now has size less by 1
				tagsize[ccindex]=tagsize[ccindex]-1
			}
		}
		
	}
	n.clean<-length(tags)
	cat("We consider ",n.clean, " tags in the final analysis, and condition upon ",length(cond), " SNPs \n")

    

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
	tagsize<-tagsize[which(unique(tagkey) %in% tags)]
	#which models are we testing?
	#the matrix of tag models only
    tagmodels<-makebinmod(n.clean,1)
    category<-apply(tagmodels,1,whichcat)
	#add the conditionals to this
	models<-matrix(1,nrow(tagmodels),2*n.clean+2*length(cond)+1)
	models[,1:n.clean]<-tagmodels[,1:n.clean]
	models[,(n.clean+length(cond)+1):(2*n.clean+length(cond))]<-tagmodels[,(n.clean+1):(2*n.clean)]
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
        snplist[[ii]]<-names(which(tagkey==tags[ii]))
    }
    return(list("postprob"=pp,"tags"=tags,"snps"=snplist,"models"=models,"bf"=logB10))
}
