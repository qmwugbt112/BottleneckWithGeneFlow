# Functions needed for analysis of bottlenecks with migration

##########################################################
# Log multinomial dirichlet calculation of the  
# likelihood of allelic composition 
##########################################################

# n.b. 		this log version deals well 
#			with small values of f and p
# Warning: 	p values must be non-zero

lmd<-function(	a=c(2,8),	# vector of subpopulation allelic counts
				f=0.1,		# Fst (describing variation around p)
				p=c(0.3,0.7)# allele frequencies 
				)
		{
	if (f==0) return(dmultinom(a,prob=p,log=TRUE))
	l<-1/f-1;n=sum(a);x<-l*p
	return(lgamma(l)+lgamma(n+1)-lgamma(n+l)	
			+sum(lgamma(x+a))	
			-sum(lgamma(a+1))
			-sum(lgamma(x))
			)	
	}	

# Vectorize the function so it can be evaluated over 
# a vector different f values
lmdL<-Vectorize(lmd,'f')
##########################################################


##########################################################
# Function return cutVals, which are values at which
# sum of x up to cutVals comprises cuts of the total 
# Used to return the levels used by contour-plotting function
##########################################################
##########################################################

returnCuts<-function(x,cuts,require_unique=T){
	sx<-sort(x)
	cs<-cumsum(sx)/sum(sx)
	cutVals<-rep(0,length(cuts))
	for (i in 1:length(cuts)) {
		cutVals[i]<-sx[min(which(cs>=cuts[i]))]
		}
	
	# code to break ties by subtracting a small amount
	# from the 1st value (ties are not allowed by filled.contour)
	if (require_unique) {
		for (i in length(cutVals):2)
		if (cutVals[i-1]>=cutVals[i])
		cutVals[i-1]<-cutVals[i]-abs(max(cutVals)/10^5)
		}
	return(cutVals)
}

##########################################################
# Fuction to draw number of surviving lineages for given 
#	number of lineages
# 	number of alleles
# and Ne in the previous generation
##########################################################

# test data

testCounts <- c(0,3,0,1,4,0,0,0)
prevNe 	<- 4.6
##########################################################


drawNoAncest <- function(counts=testCounts,Ne=prevNe){
	nLineages <- sum(counts)
	wholePart <- floor(Ne)
	prob <- c(rep(1,wholePart),Ne-wholePart)
	prob <- prob/sum(prob)
	return(sum(rmultinom(1,nLineages,prob)>0))
}  

##########################################################
# function to calculate a prior distribution on the mean number of migrants
# given that the population grew by g individuals in between 1st and 2nd
# census, and a gamma(a,b) prior on the number of migrants and 
# integrating over the uniform probability {0...g} that growth was due to
# immigration (not intrinsic growth) 
##########################################################
calcPrior<-function(xvals,			# values for mean migrants
					a=1,			# parameters of gamma prior 
					b=0.3,			# defaults to ~ mean=3.3, 95% =10,
					res=100,		# resolution for numerical integration
					g=4			# growth between censuses
					)
				{
					
				priorVals <- rowMeans(
										sapply(0:g*res/res,
											function(x) dgamma(xvals,a+x,b+1)
										   	)
								 )
				return(priorVals)		
					}


##########################################################
# function to return a sample of the counts of each allele
# in the previous generation, given the counts in this generation
# and the number of ancestral lineges 
# The sample depends on the mutation/error rate 
# which allows coalescence between non-alleles
##########################################################

drawCoal <- function(counts=testCounts,nAncest=4,u=10^-2){
	nLineages <- sum(counts)
	nAlleles <- length(counts)
	aNos <- c(1:nAlleles)	# used in sampling
	
	# if there are no coalsecent events return counts
	if (nLineages == nAncest) {
		return(list(uG=0,ancest=counts))
		}
	if (nLineages == 0) {
		return(list(uG=0,ancest=counts))
	}	
	
	# otherwise simulate required no coalescent events
	uG=0	
	for (i in 1:(nLineages-nAncest)) {
		nSame <- choose(counts,2)
		nDiff <- counts*(nLineages-counts)/2
		sProb <- nSame/(sum(nSame+nDiff))*(1-u)
		dProb <- nDiff/(sum(nSame+nDiff))*u
		uG <- uG + log(sum(sProb)+sum(dProb))
		choice <- sample(aNos, 1, prob=sProb+dProb)
		counts[choice] <- counts[choice] -1 	
	}
	return(list(uG=uG,ancest=counts))
}

#################################################
testDat <- cbind(
			c(0,7,0,5,9,0,0,0),
			c(rep(0.2/6,6),0.4,0.4))
dimnames(testDat)<-list(NULL,c('counts','freqs'))			
testF <- c(1:5/100)
testN <- c(12,6,3)
###################################################################
# Function to return the likelihood of an array of allele counts 
# from one locus Gdata[,1],in a subpopulation of a metapopulation.
# The likelihood is a function of the metapopulation allele 
# frequencies: Gdata[,2], and demographic parameters.
# The function works drawing a sample of possible genealogies
# (sample-size: nGenealogies).
# Genealogies are traced back until a founding generation, when ancestral
# allele frequencies are treated as multinomial dirichlet distributed.
# The distribution is specified by the metapopulation allele 
# frequencies and Fst.  Migrants are drawn from the same distribution.
# If Fst is an array, the function returns a likelihood estimate 
# for each Fst value
###################################################################
drawDistnFounders3<-function(
	Ne2N=0.5,			# ratio of Ne to census size
	meanMe=2.1,			# the average no diploid migrants per gen
	Fst=testF,			# Fst describing variance of ancestral p around metapop
	pastN=testN,		# a vector of diploid census counts (present to founding)
	gData=testDat,		# matrix 	1st col counts of each allele 
						#			2nd col allele frequency (metapopulation) 
	uRate=10^-2,		# relative prob non-alllelic coalescence (error & mutn.)
	nGenealogies=10^3,	# number of samples to draw
	convergePlot=F		# produce plot to assess number of samples needed
	){
	# Check allele array does not contain zero values
	if (any(gData[,2]==0)) stop("Allele frequencies must be non-zero")
	
		
	# Inital calculations from input variables
	nGens			<- length(pastN)			# no gens back to founding
	Ne				<- Ne2N*(pastN)*2			# convert each popn size to Ne
	pMigrant		<- meanMe/pastN				# proportion of migrants each gen
	pMigrant[nGens]	<- 1 						# all founders were immigrants	
	nAlleles		<- nrow(gData)				# num alleles (present or not)
	nLineages		<-sum(gData[,1])			# num allele copies
	
	if (any(pMigrant>1)) {
		stop(paste(
			'\n No. of migrants exceeds population size in generation ',
			(nGens:1)[pMigrant>1]
			)
		)
	}

	
			
	# Set to zero the variables to store running totals 
	# if covergence is being monitored
	if (convergePlot){ # detailed samples to monitor convergence
		dG<-matrix(0,nrow=nGenealogies,ncol=length(Fst))
		nG<-matrix(0,nrow=nGenealogies,ncol=length(Fst))
		}
		
	# numerator and denominator of the estimates 
	# (one for each Fst value)	
	numer <-rep(0,length(Fst)) 
	denom <-rep(0,length(Fst)) 

	
	
	for (iter in 1:nGenealogies){
		
		# zero the count of immigrant alleles
		migAlleles <- rep(0,nAlleles)
		
		# zero the probability count
		uP<-0	
		
		ancestAlleles <- gData[,1]
		for (gen in 1:nGens){

		
			ancestNo <- drawNoAncest( counts=ancestAlleles,
										   Ne=Ne[gen])
			
			# get allele counts after coalescence (& prob)
			draw <- drawCoal(counts=ancestAlleles,
								nAncest=ancestNo,
								u=uRate)
			
			uP <- uP + draw$uG
			ancestAlleles <- draw$ancest
			migrants <- rbinom( nAlleles,
								ancestAlleles,
								prob=pMigrant[gen] 
								)
			
			# record migrants
			migAlleles 		<- migAlleles + migrants
			# remove migrants from further consideration
			ancestAlleles 	<- ancestAlleles - migrants
			
			# for next gen start from the remaining ancestors
			print(ancestAlleles)
			} # gen
		
		
			
	
	print('\n');print(migAlleles)
	
	print(uP)
	print(lmdL(a=migAlleles,f=Fst,p=gData[,2]))	
	
	denom <-denom + 1
	numer <-numer + 
				exp(uP +
					lmdL(	a=migAlleles+ancestAlleles,	
							f=Fst,
							p=gData[,2])
					)
	
	
	# Keep running totals for convergence monitoring
	if (convergePlot){
		dG[iter,]<-denom							
		nG[iter,]<-numer
		}
	
	
	} # end each genealogy
	
	
	if (convergePlot){ # plot graphs to show convergence
		halfway<-nGenealogies%/%2
		yTop<-max(log(nG/dG)[halfway:nGenealogies,]
					,na.rm=T)*0.9
		yBot<-min(log(nG/dG)[halfway:nGenealogies,]
					,na.rm=T)*1.1
		plot(	1:nGenealogies,log(nG/dG)[,1],
				ylim=c(yBot,yTop),
				type='l',
				main="Convergence of likelihood estimate",
				ylab='Estimate',
				xlab='Sample Size (nGenealogies)')
		if ( length(Fst)>1)for (j in 2:length(Fst)) {
			lines(1:nGenealogies,log(nG/dG)[,j],col=j)
			}
			
		}
		
	return(ratio=numer/denom)
			
	} # end function

dd <- Vectorize(drawDistnFounders3,'meanMe')


