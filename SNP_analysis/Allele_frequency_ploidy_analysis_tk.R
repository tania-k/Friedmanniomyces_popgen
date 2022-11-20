##this script runs allele variation analysis to determine support for diploidy for nine A.fumigatus strains with full length alignments to both MAT1 and MAT2
##Last updated 17.November.2022

##load libraries
library('vcfR')
library('pinfsc50')
library('ggpubr')
library('gridGraphics')
library('grid')


#sessionInfo() 

##set seed for reproducibility
set.seed(101)

##read in data
vcf<- read.vcfR("vcf/CCFEE_5001_v1.All.SNP.combined_selected.vcf.gz") #W/o TE's and w 260 strains
#Processed variant: 813061

knitr::kable(vcf@gt[c(1:2,11,30),1:4])

#depth each allele was sequence at 
ad <- extract.gt(vcf, element = 'AD')

#use function mansplit to extract first and second allele
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)
allele3 <- masplit(ad, record = 3)
ad1 <- allele1 / (allele1 + allele2)
ad2 <- allele2 / (allele1 + allele2)
ad12 <- allele1 / (allele1 + allele2 + allele3)
ad13 <- allele3 / (allele1 + allele2 + allele3)
ad23 <- allele2 / (allele1 + allele2 + allele3)

#plot allele frequencies in hist.

#first 5001 control:
pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5001"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5001", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5001"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5001"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5001"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5001"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5001_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5001_all
#very clean, either matches ref or does not- haploid. 

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5193"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5193", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5193"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5193"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5193"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5193"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5193_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5193_all
#one small shoulder - diploid

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5195"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5195", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5195"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5195"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5195"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5195"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5195_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5195_all
#one small shoulder - diploid

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5199"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5199", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5199"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5199"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5199"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5199"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5199_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5199_all
#two shoulders - polyploidy (tetraploid?)

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5200"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5200", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5200"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5200"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5200"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5200"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5200_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5200_all
#two shoulders - polyploidy (tetraploid?)

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5208"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5208", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5208"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5208"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5208"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5208"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5208_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5208_all
#shoulder - diploid

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_524"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_524", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_524"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_524"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_524"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_524"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_524_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_524_all
#slight bump at 1/2 

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5273"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5273", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5273"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5273"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5273"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5273"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5273_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5273_all
#slight bump at 1/2 
#trim and investigate!

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5275"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5275", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5275"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5275"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5275"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5275"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5275_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5275_all
#slight bump at 1/2 
#trim and investigate!

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5277"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5277", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5277"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5277"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5277"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5277"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5277_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5277_all
#slight shoulders

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5281"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5281", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5281"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5281"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5281"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5281"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5281_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5281_all
#very slight shoulders

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5283"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5283", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5283"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5283"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5283"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5283"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5283_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5283_all
#trim and investigate

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5307"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5307", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5307"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5307"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5307"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5307"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5307_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5307_all

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5311"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5311", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5311"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5311"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5311"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5311"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5311_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5311_all

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5311_v1"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5311_v1", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5311_v1"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5311_v1"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5311_v1"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5311_v1"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5311_v1_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5311_v1_all

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_5486"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_5486", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_5486"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_5486"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_5486"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_5486"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5486_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5486_all

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_6074"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_6074", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_6074"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_6074"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_6074"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_6074"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6074_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6074_all

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_6081"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_6081", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_6081"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_6081"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_6081"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_6081"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6081_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6081_all

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_6082"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_6082", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_6082"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_6082"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_6082"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_6082"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6082_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6082_all

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_6096"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_6096", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_6096"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_6096"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_6096"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_6096"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6096_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6096_all

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_6249"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_6249", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_6249"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_6249"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_6249"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_6249"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6249_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6249_all

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_6250"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_6250", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_6250"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_6250"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_6250"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_6250"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6250_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6250_all

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_6464"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_6464", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_6464"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_6464"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_6464"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_6464"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6464_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6464_all

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_670"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_670", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_670"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_670"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_670"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_670"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_670_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_670_all

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"CCFEE_690"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "CCFEE_690", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"CCFEE_690"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12[,"CCFEE_690"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13[,"CCFEE_690"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_690"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_690_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_690_all

#remove homozygotes that overwhelm the plot, to focus on heterozygotes.
ad2_df<- data.frame(ad2)
ad2_df[ad2_df ==0] <- NA
ad1_df<- data.frame(ad1)
ad1_df[ad1_df ==1] <- NA
ad12_df<- data.frame(ad12)
ad12_df[ad12_df ==1] <- NA
ad13_df<- data.frame(ad13)
ad13_df[ad13_df ==1] <- NA
ad23_df<- data.frame(ad23)
ad23_df[ad23_df ==1] <- NA

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5001"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5001"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5001"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5001"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5001"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5001_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5001_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5193"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5193"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5193"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5193"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5193"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5193_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5193_zoom
#one small shoulder - diploid

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5195"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5195"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5195"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5195"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5195"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5195_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5195_zoom
#one small shoulder - diploid

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5199"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5199"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5199"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5199"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5199"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5199_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5199_zoom
#two shoulders - polyploidy (tetraploid?)

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5200"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5200"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5200"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5200"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5200"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5200_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5200_zoom
#two shoulders - polyploidy (tetraploid?)

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5208"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5208"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5208"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5208"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5208"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5208_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5208_zoom
#shoulder - diploid

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_524"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_524"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_524"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_524"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_524"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_524_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_524_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5273"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5273"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5273"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5273"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5273"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5273_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5273_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5275"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5275"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5275"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5275"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5275"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5275_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5275_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5277"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5277"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5277"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5277"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5277"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5277_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5277_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5281"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5281"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5281"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5281"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5281"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5281_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5281_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5283"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5283"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5283"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5283"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5283"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5283_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5283_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5307"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5307"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5307"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5307"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5307"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5307_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5307_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5311"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5311"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5311"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5311"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5311"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5311_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5311_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5311_v1"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5311_v1"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5311_v1"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5311_v1"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5311_v1"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5311_v1_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5311_v1_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_5486"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_5486"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_5486"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_5486"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_5486"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_5486_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_5486_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_6074"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_6074"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_6074"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_6074"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_6074"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6074_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6074_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_6081"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_6081"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_6081"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_6081"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_6081"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6081_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6081_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_6082"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_6082"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_6082"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_6082"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_6082"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6082_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6082_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_6096"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_6096"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_6096"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_6096"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23[,"CCFEE_6096"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6096_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6096_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_6249"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_6249"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_6249"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_6249"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_6249"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6249_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6249_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_6250"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_6250"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_6250"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_6250"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_6250"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6250_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6250_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_6464"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_6464"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_6464"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_6464"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_6464"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_6464_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_6464_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_670"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_670"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_670"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_670"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_670"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_670_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_670_zoom

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"CCFEE_690"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab = NA, ylab = NA)
hist(ad1_df[,"CCFEE_690"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad12_df[,"CCFEE_690"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad13_df[,"CCFEE_690"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
hist(ad23_df[,"CCFEE_690"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
CCFEE_690_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
CCFEE_690_zoom

#plot all 
dev.off()
p<-ggarrange(CCFEE_5001_all, 
             CCFEE_5193_all,
             CCFEE_5195_all,
             CCFEE_5199_all,
             CCFEE_5200_all,
             CCFEE_5208_all,
             CCFEE_524_all,
             CCFEE_5273_all,
             CCFEE_5275_all,
             CCFEE_5277_all,
             CCFEE_5281_all,
             CCFEE_5283_all,
             CCFEE_5307_all,
             CCFEE_5311_all,
             CCFEE_5311_v1_all,
             CCFEE_5486_all,
             CCFEE_6074_all,
             CCFEE_6081_all,
             CCFEE_6082_all,
             CCFEE_6096_all,
             CCFEE_6249_all,
             CCFEE_6250_all,
             CCFEE_6464_all,
             CCFEE_670_all,
             CCFEE_690_all,
             ncol = 3, nrow = 4)
ggsave("ploidy_by_allele_diversity_all.pdf",p, width=28, height=32, units="in")

#plot zoom
q<-ggarrange(CCFEE_5001_zoom, 
             CCFEE_5193_zoom,
             CCFEE_5195_zoom,
             CCFEE_5199_zoom,
             CCFEE_5200_zoom,
             CCFEE_5208_zoom,
             CCFEE_524_zoom,
             CCFEE_5273_zoom,
             CCFEE_5275_zoom,
             CCFEE_5277_zoom,
             CCFEE_5281_zoom,
             CCFEE_5283_zoom,
             CCFEE_5307_zoom,
             CCFEE_5311_zoom,
             CCFEE_5311_v1_zoom,
             CCFEE_5486_zoom,
             CCFEE_6074_zoom,
             CCFEE_6081_zoom,
             CCFEE_6082_zoom,
             CCFEE_6096_zoom,
             CCFEE_6249_zoom,
             CCFEE_6250_zoom,
             CCFEE_6464_zoom,
             CCFEE_670_zoom,
             CCFEE_690_zoom,
             ncol = 3, nrow = 4)
ggsave("ploidy_by_allele_diversity_zoom.pdf",q, width=18, height=22, units="in")