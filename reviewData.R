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
