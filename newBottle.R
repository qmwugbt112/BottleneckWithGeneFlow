source('https://raw.githubusercontent.com/qmwugbt112/BottleneckWithGeneFlow/master/bottleneckFunctions.R')

# print out default parameters
testN
testF

# test data
tDat <- synthDat(numLoci=20)

pt <- proc.time()
testM <- lapply(tDat,
				function(x) dd(
					gData=x,
					meanMe=seq(0,6,,10),
					Fst=testF,
					pastN=testN,
					uRate=10^-5,
					nGenealogies=10^3)				
				)
proc.time()-pt

plotLike(testM,seq(0,6,,10))


# Read in genetic data ###
genDat<-read.csv('https://raw.githubusercontent.com/qmwugbt112/BottleneckWithGeneFlow/master/RNP%20Allele%20data.csv')

#pop<-5 # Wiesbaden
#pastN<-c(1950,517,370,124,22,2)


pop<-13 # choose Rotterdam for analysis

# What are the past pop sizes ?
rotN<-c(399,50,40,15,10,4)



# extract allele frequencies from the other populations
# the 1 will generate Laplace estimates when normalized
pVals<-rowSums(genDat[,-(c(1:2,pop))])+1	
nLoci<-nlevels(genDat$Locus)

# Parameter values to explore
mValues<-seq(0,10,,11)					# number of migrants
rotFst<-seq(0.0,0.02,,7)				# Fst values

# Create list containing the genotype data
gList<-list()
for (l in levels(genDat$Locus)){
	gList[[l]]<-subset(genDat,Locus==l)[,pop]
	gList[[l]]<-cbind(gList[[l]],subset(pVals,genDat$Locus==l))
	gList[[l]][,2]<-gList[[l]][,2]/sum(gList[[l]][,2])
}

pt <- proc.time()
testM2 <- lapply(gList,
				function(x) dd(
					gData=x,
					meanMe=mValues,
					Fst=rotFst,
					pastN=rotN,
					uRate=10^-5,
					nGenealogies=10^3)				
				)
proc.time()-pt

plotLike(testM2,mValues)
