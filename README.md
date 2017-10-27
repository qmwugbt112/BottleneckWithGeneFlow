# BottleneckWithGeneFlow
R Code and example data from European Parakeets

BottleneckFunctions.R contains a series of functions that can be used in the analysis of genetic data from post-bottleneck populations, to infer the rate of gene flow that has suplemented the population since the bottleneck.

The logic is laid out in BottleneckFigure.pdf.  The analysis requires an estimate of the poulation sizes at foundation and in each generation since the bottleneck.  It then estimates whether the population grew from the founders only, or otherwise the mean per-generation immigration rate, since foundation.

An example script using the functions is in NewBottle.R

The data for the example are in RNP Allele data.csv
