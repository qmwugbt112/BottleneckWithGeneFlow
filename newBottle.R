
source('~/Dropbox/Research by person/Jim Groombridge/bottleneckFunctions.R', chdir = TRUE)
# Read in genetic data ###

genDat<-read.csv('~/Dropbox/Research by person/jim groombridge/RNP Allele data.csv')


# Initial visualisation of the differentiation from other populations
fValues	 <- 1:600/1000	# Fst range
numPops <- 18

likeCurves <- matrix(0,numPops,length(fValues))

for (p in 1:numPops){
	c<-2+p
	# the 1 will generate Laplace estimates when normalized
	pVals<-rowSums(genDat[,-(c(1:2,p))])+1
	for (loc in levels(genDat$Locus)){
		choice <- genDat$Locus==loc
		pVec <- pVals[choice]
		pVec <- pVec/sum(pVec)
		aVec <- genDat[choice,c]
		# test for low sample size bias
		if (sum(aVec)>0){nAlleles <- sum(choice)
			subVec <- sample(1:nAlleles,
							size=4,
							prob=aVec,
							replace=T)
			aVec <- rowSums(diag(nAlleles)[,subVec])
			}					
		likeCurves[p,] <- likeCurves[p,] + lmdL(aVec,fValues,pVec)
	}			
}

popNames <- names(genDat)[-(1:2)]
nPops <- length(popNames)
colourNums <- rainbow(9)

par(mfrow=c(3,1))
for (s in c(0,6,12)){
 plot(	fValues,
		likeCurves[s+1,]-max(likeCurves[s+1,]),
		ylim=c(-6,0),
		type='l',
		ylab='Support (log likelihood)',
		xlab='Fst',
		main='Likelihood curves for Fst (each city vs rest)',
		col=colourNums[3+1])
 text(x=0.34,y=-1.8,popNames[s+1],col=colourNums[3+1])		
 for (p in 2:6){
	 lines(	fValues,
			likeCurves[p+s,]-max(likeCurves[p+s,]),
			col=colourNums[3+p])
			text(x=0.34,y=-1-p*0.8,popNames[p+s],col=colourNums[3+p])		
}		
}
#pop<-5 # Wiesbaden
#pastN<-c(1950,517,370,124,22,2)

pastN<-c(399,50,40,15,10,4)
pop<-13 # Rotterdam



# extract allele frequencies from the other populations
# the 1 will generate Laplace estimates when normalized
pVals<-rowSums(genDat[,-(c(1:2,pop))])+1	
nLoci<-nlevels(genDat$Locus)

# Parameter values to explore
mValues<-seq(0,30,,15)				# number of migrants
Fst<-seq(0.0,0.02,,7)				# Fst values

# Create list containing the genotype data
gList<-list()
for (l in levels(genDat$Locus)){
	gList[[l]]<-subset(genDat,Locus==l)[,pop]
	gList[[l]]<-cbind(gList[[l]],subset(pVals,genDat$Locus==l))
	gList[[l]][,2]<-gList[[l]][,2]/sum(gList[[l]][,2])
}

quartz()
# Test convergence
testV<-dd(
			Ne2N=0.8,			
			meanMe=10,	
			Fst=Fst,				
			pastN=pastN,				
			gData=gList[[1]],		
			uRate=10^-3,
			nGenealogies=10^4,
			convergePlot=T	
			)

# Test a range of parameter values

pt<-proc.time()
testM<- dd(		gData=gList[[1]],
					meanMe=mValues,
					Fst=Fst,
					pastN=pastN,
					uRate=10^-3,
					nGenealogies=10^2)
proc.time()-pt

################## Draw countour plot of results
filled.contour(	x=Fst,
				y=mValues,
				z=testM,
				xlab='Fst',
				ylab='Migrants per generation',
				main=paste('Likelihood surface'),	
				levels=returnCuts(testM,c(0,0.05,0.1,0.9,1)),
				key.axes=axis(
					4,
					labels=c(0,0.05,0.1,0.9,1),
					at=returnCuts(testM,c(0,0.05,0.1,0.9,1))
					)	
				)				

library(parallel)
pt<-proc.time()
# testM<-mclapply( gList,
				# ,mc.cores=length(gList)

testM <- lapply(gList[1:3],
				function(x) dd(
					gData=gList[[1]],
					meanMe=mValues,
					Fst=Fst,
					pastN=pastN,
					uRate=10^-3,
					nGenealogies=10^3)
							
				)
proc.time()-pt

Ne2Nv<-c(1:5)*5
MeV<-c(0,5,10,15,20,25)*2

calcMat2<-function(x){
	tMat<-matrix(0,nrow=length(MeV),ncol=length(Ne2Nv))
	for(n in 1:length(Ne2Nv)) for(m in 1:length(MeV)){
	tMat[m,n]<-drawDistnFounders2(
			Ne2N=Ne2Nv[n],			
			meanMe=MeV[m],	
			Fst=0.02,				
			pastN=pastN,				
			gData=x,		
			uRate=10^-3,
			sB=0.5, # tune this parameter
			nGenealogies=3*10^3,
			convergePlot=F	
			)
	}
	return(tMat)		
}

library(parallel)
pt<-proc.time()
testM<-mclapply(	gList, calcMat2, mc.cores=length(gList)
)
proc.time()-pt
save(testM,file='~/Parakeet/Roterdam2')

# sum the log likelihood matrixes for each locus
combinedM<-rep(0,nrow=nrow(testM[[1]]),ncol=ncol(testM[[1]]))
for (i in 1:length(testM)){
	combinedM<-combinedM+log(testM[[i]])-max(log(testM[[i]]))
}

# express each value as a difference from ML
combinedM<-combinedM-max(combinedM)
filled.contour(	x=MeV,
 				y=Ne2Nv,
 				z=combinedM,
 				xlab='Migrants per generation',
 				ylab='Population under-estimation',
 				main=paste('Rotterdam \n Log Likelihood surface'),	
 				#levels=returnCuts(combinedM,c(0,0.05,0.1,0.9,1)),

 					)	
 							
combinedM<-exp(combinedM)


filled.contour(	x=mValues,
 				y=Fst,
 				z=combinedM,
 				xlab='Migrants per generation',
 				ylab='Fst',
 				main=paste('Wiesbaden \n Likelihood surface'),	
 				#levels=returnCuts(combinedM,c(0,0.05,0.1,0.9,1)),
 				key.axes=axis(
 					4,
 					labels=c(	0,
 								0.05,
 								0.1,
 								0.9,
 								1),
 					at=returnCuts(combinedM,c(
 								0,
 								0.05,
 								0.1,
 								0.9,
 								1))
 					)			
 				)

normM<-lapply(testM,function(x) log(x)-max(log(x)))
resMat<-matrix(	0,
					nrow=length(mVals),
					ncol=length(Fst)
					)
for (l in 1:length(normM)) resMat<-resMatnormM[[l]]
resMat<-exp(resMat-max(resMat))


################## Draw countour plot of results
filled.contour(	x=mValues,
				y=Fst,
				z=resMat,
				xlab='Migrants per generation',
				ylab='Fst',
				main=paste('Wiesbaden \n Likelihood surface'),	
				#levels=returnCuts(resMat,c(0,0.05,0.1,0.9,1)),
				key.axes=axis(
					4,
					labels=c(0,0.05,0.1,0.9,1),
					at=returnCuts(resMat,c(0,0.05,0.1,0.9,1))
					)			
								
				)
				
				

