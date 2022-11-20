library(gdsfmt)
library(SNPRelate)
library(dplyr)
library(wordcloud)
library(tm)
library(ggfortify)
library(ggplot2)
library(plotly)
library(phyloseq)
library(ggrepel)
library(ggbiplot)
library(RColorBrewer)


gdsfile = "plots/snps_selected.gds"
vcf.fn <- "vcf/CCFEE_5001_v1.All.SNP.combined_selected.vcf.gz"

if(!file.exists(gdsfile)){
	snpgdsVCF2GDS_R(vcf.fn, gdsfile,method="biallelic.only")
	                #option=snpgdsOption(CM002236=1,CM002237=2,CM002238=3,CM002239=4,CM002240=5,CM002241=6,CM002242=7))
}

snpgdsSummary(gdsfile)
genofile <- snpgdsOpen(gdsfile)
chroms <- read.gdsn(index.gdsn(genofile,"snp.chromosome"))
#get.attr.gdsn(index.gdsn(genofile, "snp.chromosome"))
chr <- strtoi(sub("SCAF_([0-9]+)","\\1",chroms,perl=TRUE))

pca <- snpgdsPCA(genofile,num.thread=2,autosome.only=FALSE)

pc.percent <- pca$varprop*100

id <- pca$sample.id
#id <- c("5001", "5193", "5195", "5199", "5200", "5208", "524", "5273", "5275", "5277", "5281", "5283", "5307", "5311", "5311_v1", "5486", "6074", "6081", "6082", "6096", "6249", "6250", "6464", "670", "690")

#pca$sample.id = id

head(round(pc.percent, 2))
pdf("plots/PCA_snp_plots5.pdf")
tab <- data.frame(sample.id = pca$sample.id,
                 # pop = pheno$MinimalMediaGrowth,
                  EV1=pca$eigenvect[,1], # PCA vector 1
                  EV2=pca$eigenvect[,2], # PCA vector 2
		  stringsAsFactors=FALSE)

#p1<-plot(tab$EV2, tab$EV1,
     #, col=as.integer(tab$pop),
#xlab="eigenvector 2", ylab="eigenvector 1", main="PCA SNP plot")
#text(x = pca$eigenvect[,2], y = pca$eigenvect[,1], labels = tab$sample.id, pos = 1 ,cex =0.8, offset = 0.5)

#mx <- apply(tab$EV2,5,max)
#mn <- apply(tab$EV1,5,min)
#https://blog.fellstat.com/?cat=11
p2 <-textplot(tab$EV2, tab$EV1, id, cex=1, new=TRUE, show.lines=TRUE, xlim=c(-0.15,0.65),ylim=c(-0.30,0.40))

set.seed(100)
# recode the snp.gds to support chromosomes?
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2,autosome.only=FALSE))
rv <- snpgdsCutTree(ibs.hc)
plot(rv$dendrogram, leaflab="none", main="Friedmanniomyces endolithicus Strains")

snpgdsDrawTree(rv, type="z-score", main="Friedmanniomyces endolithicus Strains")
snpgdsDrawTree(rv, main="Friedmanniomyces endolithicus Strains",
               edgePar=list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"))

table(rv$samp.group)
df = data.frame(group = rv$samp.group)
rownames(df) = pca$sample.id
write.csv(df,"plots/popset_inferred.csv")
tab <- data.frame(sample.id = pca$sample.id,
                  pop = rv$samp.group,
                  EV1=pca$eigenvect[,1], # PCA vector 1
                  EV2=pca$eigenvect[,2], # PCA vector 2
                  stringsAsFactors=FALSE)
#plot(tab$EV2, tab$EV1,
#     col=as.integer(tab$pop),
#     xlab="eigenvector 2", ylab="eigenvector 1", main="PCA SNP plot")

#CORRSNP <- snpgdsPCACorr(pca, genofile, eig.which=1:4,num.thread=2)

#savepar <- par(mfrow=c(3,1), mai=c(0.3, 0.55, 0.1, 0.25))
#for (i in 1:3)
#{
#  plot(abs(CORRSNP$snpcorr[i,]), ylim=c(0,1), xlab="", ylab=paste("PC", i),
#       col=factor(chr), pch="+")
#}

options(ggrepel.max.overlaps = Inf)

#http://rstudio-pubs-static.s3.amazonaws.com/53162_cd16ee63c24747459ccd180f69f07810.html
#https://ggrepel.slowkow.com/articles/examples.html

metadata <- read.csv("metadata.tsv", sep = "\t", row.names = 1, header = TRUE)
metadata2 <- read.csv("metadata.tsv", sep = "\t")
dataset_pcr <- metadata2[2:3]
pca_res <- prcomp(dataset_pcr, scale. = TRUE)

biplot = ggbiplot(pcobj = pca_res,
                  choices = c(1,2),
                  #label.position = "identity",
                  #label.repel = TRUE,
                  obs.scale = 1, 
                  var.scale = 1,  # Scaling of axis
                  #labels = row.names(metadata), # Add labels as rownames
                  #labels.size = 3.5,
                  varname.size = 3.5,
                  varname.abbrev = FALSE,  # Abbreviate variable names (TRUE)
                  var.axes = FALSE,      # Remove variable vectors (TRUE)
                  circle = FALSE,        # Add unit variance circle (TRUE)
                  ellipse = TRUE, 
                  groups = metadata$Ploidy) # Adding ellipses
#print(biplot)

#repel that shit https://stackoverflow.com/questions/68738673/how-to-repel-labels-in-ggplot
#turn off the legend to get the secret a's to disappear from legend https://stackoverflow.com/questions/18337653/remove-a-from-legend-when-using-aesthetics-and-geom-text
biplot2 = biplot + labs(title = "PCA of SNPs generated with 25 strains of \nFriedmanniomyces endolithicus",
                        colour = "Ploidy") + theme_bw() + expand_limits(x = c(-1, 3.5)) + 
                        geom_label_repel(mapping = aes(label = row.names(metadata), color = factor(metadata$Ploidy)), show.legend = FALSE) +
                        scale_color_manual(values = c("Diploid" = "purple", "Haploid"="orange","Triploid"="steelblue")) 
                        
#+ xlim = c(-2.0,2.65) + ylim = c(-1.5,2.50)
print(biplot2)

biplot3 = ggbiplot(pcobj = pca_res,
                   choices = c(1,2),
                   #label.position = "identity",
                   #label.repel = TRUE,
                   obs.scale = 1, 
                   var.scale = 1,  # Scaling of axis
                   #labels = row.names(metadata), # Add labels as rownames
                   #labels.size = 3.5,
                   varname.size = 3.5,
                   varname.abbrev = FALSE,  # Abbreviate variable names (TRUE)
                   var.axes = FALSE,      # Remove variable vectors (TRUE)
                   circle = FALSE,        # Add unit variance circle (TRUE)
                   ellipse = TRUE, 
                   groups = factor(metadata2$Year)) # Adding ellipses
#print(biplot3)

#repel that shit https://stackoverflow.com/questions/68738673/how-to-repel-labels-in-ggplot
biplot4 = biplot + labs(title = "PCA of SNPs generated with 25 strains of \nFriedmanniomyces endolithicus",
                        colour = "Year") + theme_bw() + expand_limits(x = c(-1, 3.5)) + 
                        geom_label_repel(mapping = aes(label = row.names(metadata), color = factor(metadata$Year)), show.legend = FALSE) +
                        scale_color_manual(values = c("Diploid" = "purple", "Haploid"="orange","Triploid"="steelblue", "1981" = "#FF5733", "1997" = "#7AFF33", "2004" = "#33FCFF", "2010" = "#FCFF33", "2011" = "#336BFF", "2016" = "#FF33AC")) 
print(biplot4)

biplot5 = ggbiplot(pcobj = pca_res,
                   choices = c(1,2),
                   #label.position = "identity",
                   #label.repel = TRUE,
                   obs.scale = 1, 
                   var.scale = 1,  # Scaling of axis
                   #labels = row.names(metadata), # Add labels as rownames
                   #labels.size = 3.5,
                   varname.size = 3.5,
                   varname.abbrev = FALSE,  # Abbreviate variable names (TRUE)
                   var.axes = FALSE,      # Remove variable vectors (TRUE)
                   circle = FALSE,        # Add unit variance circle (TRUE)
                   ellipse = TRUE, 
                   groups = factor(metadata$Elevation)) # Adding ellipses
#print(biplot)

#repel that shit https://stackoverflow.com/questions/68738673/how-to-repel-labels-in-ggplot
biplot6 = biplot + labs(title = "PCA of SNPs generated with 25 strains of \nFriedmanniomyces endolithicus",
                        colour = "Elevation") + theme_bw() + expand_limits(x = c(-1, 3.5)) + 
                        geom_label_repel(mapping = aes(label = row.names(metadata), color = factor(metadata$Elevation)), show.legend = FALSE) +
                        scale_color_manual(values = c("Diploid" = "purple", "Haploid"="orange","Triploid"="steelblue", "200" = "#FF5733" , "1000" = "#7AFF33", "2000" = "#33FCFF",  "3000" = "#FCFF33")) 
print(biplot6)
