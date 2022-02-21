#Evolution of pesticide tolerance and associated changes in the microbiome in Daphnia magna
#Lizanne Janssens, Marlies Van de Maele, Vienna Delnat, Charlotte Theys, Shinjini Mukherjee, Luc De Meester, Robby Stoks
#Journal (Year)
#R code tested on 21/02/2022


####### Install packages ####### 

#Installing Phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("decontam")

#Install package that imports QIIME artifacts into phyloseq - requires devtools first 
install.packages("devtools")
install.packages("usethis")
library(usethis)
library(devtools)
devtools::install_github("jbisanz/qiime2R")   

## Installer Code for Packages Below If You Don't Have Them 
install.packages(c("vegan","microbiome","RColorBrewer","ggplot2","gplots","cowplot"))
install.packages(c("labdsv","tidyverse","ggfittext","dplyr","stats","ggpubr"))
install.packages(c("lme4","car","effects","emmeans","MuMIn","lmerTest","lattice","afex","DCA","picante"))


####### Load packages ####### 

## Microbiome Tools
library(phyloseq)
library(decontam)
library(vegan)
library(qiime2R) 

##Plotting     
library(ggplot2)
library(RColorBrewer)
library(ggfittext)
library(ggpubr)

#Data Manipulation
library(dplyr)

#Stats Models
library(lme4)
library(car)
library(effects)
library(emmeans)
library(lmerTest)
library(afex)
library(stats)
library(picante)

#https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html
#https://microbiome.github.io/tutorials/

sessionInfo()


####### Read data files as physeq ####### 

#Session --> Set Working Directory --> To Source File Location

## Build Phyloseq Object 
physeq<-qza_to_phyloseq(features="table_16S_trim14trunc240.qza",taxonomy="taxonomy_16S_trim14trunc240.qza", 
                        metadata="sample-metadata-R.txt",tree="rooted-tree_16S_trim14trunc240.qza")
physeq

## Access each of the different elements of the phyloseq object
#Sample Data 
head(sample_data(physeq))
colnames(sample_data(physeq))
#Taxonomy
head(tax_table(physeq))
rank_names(physeq)
#ASV Matrix
head(otu_table(physeq))


####### Decontamination ####### 

#https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
#Prevalence (absence/presence) used as method 
#(not Frequency, no DNA concentration per sample available)
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="neg")
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant)
#no decomtamination found (FALSE 1455, integer(0))

#Remove negative controls from dataset
physeqDC <- subset_samples(physeq,neg==FALSE)
sample_data(physeqDC)

#ASVs within negative control samples (blankDNA, blankPCR)
#These are not removed, since not seen as contamination by decontam
physeqBlanco <- subset_samples(physeq,neg==TRUE)
sample_data(physeqBlanco)
physeqBlanco <-prune_taxa(taxa_sums(physeqBlanco)>0,physeqBlanco)
tax_table(physeqBlanco)


####### Cleaning and Filtering ####### 

#Export taxa_table
physeqDCDF=as.data.frame(tax_table(physeqDC))
# write.table(physeqDCDF, "physeqDC.txt", sep="\t")
#Show available ranks in the dataset
rank_names(physeqDC)

#Replace "uncultured*", "metagenome", "microbial*" and "unidentified*" as NA values in Species
NAlistSpecies= c("uncultured*", "metagenome", "microbial*", "unidentified*")
for (pattern in NAlistSpecies) {tax_table(physeqDC)[, "Species"] <- gsub(tax_table(physeqDC)[, "Species"], pattern = pattern, replacement =  NA)}

#Percentage of NA in Species - 87.42%
Species = tax_table(physeqDC)[,7]
apply(Species,MARGIN=2,function(x){round(mean(is.na(x))*100, digits=2)})

#Note, majority (87.42%) is NA in Species --> remove this rank
tax_table(physeqDC) <- tax_table(physeqDC)[,2:6]
#Note, NA in Phylum (no "uncharacterized", "uncultured", "unidentified") --> Remove these ASVs (2)
physeqFilter1 <- subset_taxa(physeqDC, !is.na(tax_table(physeqDC)[,1]))
physeqFilter1
#Control if there are indeed 0% NA values for the Phylum taxonomy rank, ok.
apply(tax_table(physeqFilter1)[,1],2,function(x){round(mean(is.na(x))*100, digits=2)}) 

#Note, environmental ("Lineage_IV"), plant (C), hypersaline ("vadinHA49"), 
#soil ("WD2101_soil_group", "WPS-2", "Blfdi19") or marine (CL500-29_marine_group", 
#"NS11-12_marine_group", "NS9_marine_group", "SAR324_clade(Marine_group_B)", 
#"PeM15", "mle1-27", "Pla3_lineage", "Yoonia-Loktanella", "SM2D12", "OM60(NOR5)_clade"), 
#"Chloroplast" and "Mitochondria" ASVs in Genus --> Remove these ASVs (96)
physeqFilter2 <- subset_taxa(physeqFilter1, 
                             !Genus %in% c("Mitochondria", "Chloroplast", "Lineage_IV", "LWQ8", "vadinHA49", 
                                           "WD2101_soil_group", "WPS-2", "Blfdi19", "CL500-29_marine_group", 
                                           "NS11-12_marine_group", "NS9_marine_group", 
                                           "SAR324_clade(Marine_group_B)", "PeM15","mle1-27", "Pla3_lineage", 
                                           "Yoonia-Loktanella", "SM2D12", "OM60(NOR5)_clade", "SM1A07", 
                                           "[Eubacterium]_brachy_group"))
physeqFilter2

#Note, "uncultured" in Genus, Family and Order; also "Unknown_Family" in Family; 
#NA values --> Replace these with nothing but keep ASVs!
NAlistTax=c("uncultured","Unknown_Family")
for (pattern in NAlistTax){  tax_table(physeqFilter2) <- gsub(tax_table(physeqFilter2), pattern = pattern, replacement = " Unidentified")}
tax_table(physeqFilter2) <- replace(tax_table(physeqFilter2), is.na(tax_table(physeqFilter2)), " Unidentified")

#Note, wrong identification of Phylum, Class, Order, Family, Genus --> Replace these with correct identification or with nothing but keep ASVs!
WrongID=c("OM190","Candidatus_Campbellbacteria")
for (pattern in WrongID){tax_table(physeqFilter2)[, colnames(tax_table(physeqFilter2))] <- 
  gsub(tax_table(physeqFilter2)[, colnames(tax_table(physeqFilter2))], pattern = pattern, replacement = " Unidentified")}
WrongIDclass=c("Gracilibacteria","Parcubacteria")
for (pattern in WrongIDclass){tax_table(physeqFilter2)[, "Class"] <- 
  gsub(tax_table(physeqFilter2)[, "Class"], pattern = pattern, replacement =  " Unidentified")}
WrongIDorder=c("Gammaproteobacteria_Incertae_Sedis","Candidatus_Kaiserbacteria","Gracilibacteria")
for (pattern in WrongIDorder){tax_table(physeqFilter2)[, "Order"] <- 
  gsub(tax_table(physeqFilter2)[, "Order"], pattern = pattern, replacement =  " Unidentified")}
WrongIDfamily=c("Rhizobiales_Incertae_Sedis","CHAB-XI-27","Candidatus_Kaiserbacteria","Candidatus_Hepatincola",
                "Caenarcaniphilales","Babeliales","Bradymonadales","Gracilibacteria","Saccharimonadales","Kapabacteriales")
for (pattern in WrongIDfamily){tax_table(physeqFilter2)[, "Family"] <- 
  gsub(tax_table(physeqFilter2)[, "Family"], pattern = pattern, replacement =  " Unidentified")}
tax_table(physeqFilter2)[, "Family"] <- gsub(tax_table(physeqFilter2)[, "Family"], pattern = "Oligoflexales", replacement = "Oligoflexaceae")
WrongIDgenus=c("CHAB-XI-27","env.OPS_17","T34","LiUU-11-161","UBA12409",
               "Candidatus_Kaiserbacteria","Fimbriimonadaceae","Caenarcaniphilales","Bradymonadales","Babeliales","Gracilibacteria")
for (pattern in WrongIDgenus){tax_table(physeqFilter2)[, "Genus"] <- 
  gsub(tax_table(physeqFilter2)[, "Genus"], pattern = pattern, replacement =  " Unidentified")}

tax_table(physeqFilter2)["7f1bbb09801e6a094bbb6c84089da68d", "Phylum"] <- "Candidatus_Kaiserbacteria"
ASVphylumList1=c("4e65fbe7d947b881fded12ad64b0c635","28f94f62f88243cf4d36cc83cae8faba", "1fc70f2151946d9a45f322a8460b714f")
for (ASVID in ASVphylumList1){tax_table(physeqFilter2)[ASVID, "Phylum"] <- "Gracilibacteria"}
ASVphylumList2=c("43dfe1403dac33f9f4593da1c049c99a","ad2ca69176d72f3c2260e69716ebd56f","b03c2c6c84a848ac5353fdaaf7d0ddb3",
                 "6bddd141313fcc6078021b6ae597dc75","9a51bdd7e2ba389627e8dfcce957c19e")
for (ASVID in ASVphylumList2){tax_table(physeqFilter2)[ASVID, "Phylum"] <- "Candidatus_Campbellbacteria"}
tax_table(physeqFilter2)["b24bfa72082faec74ac8e287999dd182", "Class"] <- "Acidobacteriia"
tax_table(physeqFilter2)["3d38eef3781d48683af3630ab65a59e4", "Order"] <- "Burkholderiales"
tax_table(physeqFilter2)["b24bfa72082faec74ac8e287999dd182", "Order"] <- "Bryobacterales"
tax_table(physeqFilter2)["c5334bf4b6b64c392ece945829263344", "Order"] <- "Nostocales"
tax_table(physeqFilter2)["06a99d6c9c578768436e60a0049bdbd6", "Order"] <- "Vellionellales"
tax_table(physeqFilter2)["b2e29c33fbbc42cd4cec27d6bc761809", "Order"] <- " Unidentified"
ASVorderList=c("51ee97e240502848b7600fd8591d8e39","a16df1617fab6017e525a43ecda75a8f","b24c3a2da3860563348c3fc81cc6923b")
for (ASVID in ASVorderList){tax_table(physeqFilter2)[ASVID, "Order"] <- "Legionellales"}
tax_table(physeqFilter2)["c5334bf4b6b64c392ece945829263344", "Family"] <- "Calotrichaceae"
tax_table(physeqFilter2)["b24bfa72082faec74ac8e287999dd182", "Family"] <- "Bryobacteraceae"
ASVfamilyList=c("51ee97e240502848b7600fd8591d8e39","a16df1617fab6017e525a43ecda75a8f","b24c3a2da3860563348c3fc81cc6923b")
for (ASVID in ASVfamilyList){tax_table(physeqFilter2)[ASVID, "Family"] <- "Coxiellaceae"}
tax_table(physeqFilter2)["c5334bf4b6b64c392ece945829263344", "Genus"] <- "Calothrix"
ASVgenusList=c("c8d089ea661e8ab776cba132a6b87e00","c9f6c592e5502fcedc04b66177037a1a","bba7f32e41c94f194d68db95f8c43678",
               "9e811b9aa382441089f80386b979b00a","cf21388f722c158f4cfaa678bbe4e8b0","b8d49972a6fd2521bc75e97e158850b3","d295e67d72c220013425157091a56b42")
for (ASVID in ASVgenusList){tax_table(physeqFilter2)[ASVID, "Genus"] <- " Unidentified"}

#Remove singletons (ASVs with only one read) as they are more likely sequencing errors (49)
physeqFiltered<-prune_taxa(taxa_sums(physeqFilter2)>1,physeqFilter2)
physeqFiltered
#Export taxa_table
physeqFilteredDF=as.data.frame(tax_table(physeqFiltered))
# write.table(physeqFilteredDF, "physeqFiltered.txt", sep="\t")


####### Rarefying samples ####### 

#Save rarefaction curve in svg file (Appendix B - Figure B.1)
svg(filename = "FigureB1_AppendixB.svg", width = 15, height = 8)
#Plot rarefaction curve (species discovery curve) in vegan to see if plateau is reached for each sample
rarecurve(t(otu_table(physeqFiltered)), step=50, cex=0.5)
#Note, abline was added after cut-off threshold was determined below.
abline(v=13797, col="blue", lty=2)
#Stop saving code in svg file
dev.off()

#Determine what the minimum coverage is 
min(sample_sums(physeqFiltered))
#Indicates which samples would be removed with certain threshold coverage
#Set threshold as high as possible without losing samples (choose threshold below minimum)
names(sample_sums(physeqFiltered))[which(sample_sums(physeqFiltered)<13797)]
length(sample_sums(physeqFiltered)[which(sample_sums(physeqFiltered)<13797)])

#Rarefy to a custom minimum library size, set random number seed for reproducibility
#if random number seed is not set, all people would get different results (e.g. permutation), default rngseed is 711
physeqFilteredRF<-rarefy_even_depth(physeqFiltered,13797,rngseed = 711)
#203 OTUs were removed because they are no longer present in any sample after random subsampling

#Export taxa_table and otu_table
physeqFilteredRFdf=as.data.frame(tax_table(physeqFilteredRF))
# write.table(physeqFilteredRFdf, "physeqFilteredRF.txt", sep="\t")
ASVfilteredRF=otu_table(physeqFilteredRF)
# write.table(ASVfilteredRF, "ASVfilteredRF.txt", sep="\t")


####### Abundance bar plot - Genus #######

#Merge samples per sampletype-treatment combination 
#(to plot relative abundances per sampletype-treatment instead of per sample)
physeq_merged = merge_samples(physeqFilteredRF, "sampletypetreatment", fun = sum)

#Merges ASVs that have the same taxonomy at the Genus taxanomic rank (keeping the NA values)
physeq_mergedUnique <- tax_glom(physeq_merged, "Genus", NArm=FALSE)

#Absolute count --> relative abundance (ra)
physeq_ra <- transform_sample_counts(physeq_mergedUnique, function(x){round(x/sum(x),digits=4)})

#Shorten the name of "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" to fit in the legend
tax_table(physeq_ra)[, "Genus"] <- gsub(tax_table(physeq_ra)[, "Genus"], 
                                        pattern = "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", 
                                        replacement =  "(Allo-Neo-Para)-Rhizobium")

#Merge unidentified genera for the abundance plot
physeq_ra_unidentified <- subset_taxa(physeq_ra, Genus %in% c(" Unidentified"))
physeq_ra_unidentified_df <- as.data.frame(tax_table(physeq_ra_unidentified)[,5])
physeq_ra_unidentified_list <- row.names(physeq_ra_unidentified_df)
for (ASVID in physeq_ra_unidentified_list){tax_table(physeq_ra)[ASVID,] <- " Unidentified"}
physeq_ra_merged1 <- tax_glom(physeq_ra, "Genus", NArm=FALSE)

#Indicate and merge rare identified taxa: relative abundance < 0.1% (0.0010)
physeq_ra_rare <- prune_taxa(taxa_sums(physeq_ra_merged1)<0.0010,physeq_ra_merged1)
physeq_ra_rare_df <- as.data.frame(tax_table(physeq_ra_rare)[,5])
physeq_ra_rare_list <- row.names(physeq_ra_rare_df)
for (ASVID in physeq_ra_rare_list){tax_table(physeq_ra_merged1)[ASVID,] <- " Rare genera"}
physeq_ra_merged2 <- tax_glom(physeq_ra_merged1, "Genus", NArm=FALSE)

# Table showing all the different genera
table(tax_table(physeq_ra_merged2)[,5])

#Export the stacked bar plot below to an svg file
svg(filename = "Figure3_Manuscript.svg", width = 16, height = 16)
#Make a stacked bar plot using the date from pgut_ra_pruned_rel whereby the colors are different genera
plot_bar(physeq_ra_merged2, fill = "Genus") +
  #Add labels with the genus within the bar plot per color and make sure the text fits the height of the stack.
  geom_fit_text(aes(label = Genus), position = position_stack(vjust = 0.5), size = 15, show.legend = FALSE) +
  #Fit the legend into 6 columns
  guides(fill=guide_legend(ncol=6)) +
  #Set y-label to 'Relative abundance'
  ylab("Relative abundance") +
  #Change the x-axis ticks labels 
  scale_x_discrete(labels= c("gut-control"="Control\nGut","gut-pesticide"="Chlorpyrifos-selected\nGut",
                             "whole-control"="Control\nWhole body","whole-pesticide"="Chlorpyrifos-selected\nWhole body")) +
  #Transform the proportions to percentages on the y-axis and 
  #remove padding around the data that would ensure that 
  #the data were placed some distance away from the axes (no white space between bar and x-axis)
  scale_y_continuous(labels = scales::percent, expand=c(0,0.0), limits=c(0,1)) + 
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set position the legend on the bottom 
    legend.position = "bottom", 
    #Remove the legend title
    legend.title=element_blank(),
    #Set font size of the legend
    legend.text = element_text(size = 15),
    #Set font size of x-axis tick labels
    axis.text.x = element_text(size = 15),
    #Set font size of y-axis tick labels
    axis.text.y = element_text(size = 15),
    #Remove the x-axis label
    axis.title.x = element_blank(),
    #Set font size and margins of y-axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt")),
    #Remove the top axis
    axis.line.x.top = element_blank(),
    #Remove the bottom axis
    axis.line.x.bottom = element_blank(),
    #Remove the right axis
    axis.line.y.right = element_blank(),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  #Manual annotation of text within the stacked bar that is not shown by geom_fit_text 
  #as the height of the stack is not high enough
  annotate(geom="text", x=1, y=0.0708, label="Flavobacterium", color="black", size=2) +
  annotate(geom="text", x=2, y=0.140, label="Flavobacterium", color="black", size=3) 
#Stop saving code in svg file
dev.off()


####### BETA diversity - NMDS ORDINATION - Bray Curtis ####### 

#When we fit NMDS models to microbiome data, we're making no assumptions about how the data correspond to sample groups of interest. 
#Rather, we ordinate samples in n-dimensional space, corresponding to the number of axes we choose to ordinate to, and then overlay our sample metadata on top to look for patterns.
#PERMANOVA using Vegan package: PERMANOVA is a randomisation procedure that will test for the effects of predictors of interest in driving differences in beta diversity/community structure.

#Run line below to run Beta diversity analysis with all ASV
physeq = physeqFilteredRF

#Convert Phyloseq so that Vegan can use it
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
#Convert Sample Data to data frame     
sample_vegan<-as(sample_data(physeq),"data.frame")
#Convert OTU table to abundance matrix
abund_vegan<-vegan_otu(physeq)
#Set seed to get the same outcomst of the permutations
set.seed(10000)
#PERMANOVA using adonis
adonisBray<-adonis(abund_vegan ~ treatment * sampletype, data=sample_vegan, permutations=10000, method="bray")
adonisBray
#Homogeneity assumption
dist<-vegdist(abund_vegan, method="bray")
physeq.disper <- betadisper(dist, sample_vegan$treatment)
permutest(physeq.disper, pairwise = TRUE) 
#P = 0.167 --> ok; assumption for adonis is met for treatment (homogeneous dispersion)
physeq.disper <- betadisper(dist, sample_vegan$sampletype)
permutest(physeq.disper, pairwise = TRUE) 
#P = 0.001 --> not ok; differences in composition within sampletypes (heterogeneous dispersion)

#Anderson MJ, Walsh DCI. PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: What null hypothesis are you testing? 
#Ecological monographs [Internet] 2013; 83: 557. Available from: http://doi.org/10.1890/12-2010.1
#PERMANOVA (which is basically adonis()) was found to be largely unaffected by heterogeneity in Anderson & Walsh's simulations but only for balanced designs.
#For unbalanced designs PERMANOVA was too liberal if the smaller group had greater dispersion, and too conservative if the larger group had greater dispersion.

#Get NDMS scores based on k that determines the stress level (Stress < 0.2 is ok)
nmdsBray2 <- ordinate(physeq, method="NMDS",k=2, distance="bray") 
# Stress = 0.2062988 --> not ok 
goodness(nmdsBray2); stressplot(nmdsBray2) 
nmdsBray3 <- ordinate(physeq, method="NMDS",k=3, distance="bray") 
# Stress = 0.1315289 --> ok 
goodness(nmdsBray3); stressplot(nmdsBray3) 
#Extract scores manually and make your own plot
nmds_scores_test<-data.frame(scores(nmdsBray3))
#Add sample ID
nmds_scores_test$Sample<-rownames(nmds_scores_test)
#Strip out the sample data 
sample_vegan <- as(sample_data(physeq),"data.frame") 
sample_vegan$Sample<-rownames(sample_vegan)
#Add in metadata with treatments
nmds_scores_treatments<-left_join(nmds_scores_test,sample_vegan,"Sample")
#Calculate the average and the 95% confidence interval of NMDS1/NMDS2/NMDS3 and add values in a data frame
#do.call(data.frame, ...) adds the values in a data frame with each calculation as a single column (instead of as attribute) 
#cbind within aggregate lets you run the calculations on multiple variables
#~sampletypetreatment within aggregate is the grouping variable
#by using function(x) in the FUN argument of aggregate you can run multiple calculations (e.g. Avg and Conf)
nmdsBray3 <- do.call(data.frame, aggregate(cbind(NMDS1,NMDS2,NMDS3) ~ sampletypetreatment, 
                          data = nmds_scores_treatments, 
                          FUN= function(x) c(Avg=round(mean(x), digits=2), 
                                             Conf=round(1.96*(sd(x)/sqrt(length(x))), digits=2))))

#Make an NMDS plot of NMDS1 and NMDS2 - used in manuscript
Plot_nmdsBray3 <- ggplot(nmdsBray3, aes(x = NMDS1.Avg, y = NMDS2.Avg)) +
  #Set the horizontal error bars based on the average and the 95% confidence interval
  geom_errorbarh(aes(xmin = (NMDS1.Avg-NMDS1.Conf), xmax = (NMDS1.Avg+NMDS1.Conf)), height = 0) +
  #Set the vertical error bars based on the average and the 95% confidence interval
  geom_errorbar(aes(ymin = (NMDS2.Avg-NMDS2.Conf), ymax = (NMDS2.Avg+NMDS2.Conf)), width = 0) +
  #Set the shape of the points manually
  geom_point(aes(fill = sampletypetreatment), shape = c(22, 21, 22, 21), size = 5) +
  #Set the same scale on the x-axis for both Unifrac plots
  scale_x_continuous(limits=c(-0.9,0.9)) + 
  #Set the fill of the points manually
  scale_fill_manual(
    values = c("White", "White", "Black", "Black"),
    labels = c("Control Gut", "Control Whole body", 
               "Chlorpyrifos-selected Gut", "Chlorpyrifos-selected Whole body")) +
  #Fit the legend into 2 columns, indicate the legend shape and fill
  guides(fill = guide_legend(override.aes = list(shape = c(22, 21, 22, 21)), ncol = 2,
         shape = guide_legend(override.aes = list(fill = c("White", "White", "Black", "Black"))))) + 
  #Set x-label to 'NMDS1'
  xlab(expression("NMDS1")) +
  #Set y-label to 'NMDS2'
  ylab(expression("NMDS2")) +
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set position the legend on the top 
    legend.position = "top", 
    #Remove the legend title
    legend.title = element_blank(),
    #Set font size and margins of the legend
    legend.text = element_text(size = 15, margin = margin(2, 2, 10, 2,"pt")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(30, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt")),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_nmdsBray3

#Make an NMDS plot of NMDS1 and NMDS3 - NOT used in manuscript
Plot_nmdsBray3b <- ggplot(nmdsBray3, aes(x = NMDS1.Avg, y = NMDS3.Avg)) +
  #Set the horizontal error bars based on the average and the 95% confidence interval
  geom_errorbarh(aes(xmin = (NMDS1.Avg-NMDS1.Conf), xmax = (NMDS1.Avg+NMDS1.Conf)), height = 0) +
  #Set the vertical error bars based on the average and the 95% confidence interval
  geom_errorbar(aes(ymin = (NMDS3.Avg-NMDS3.Conf), ymax = (NMDS3.Avg+NMDS3.Conf)), width = 0) +
  #Set the shape of the points manually
  geom_point(aes(fill = sampletypetreatment), shape = c(22, 21, 22, 21), size = 5) +
  #Set the same scale on the x-axis for both Unifrac plots
  scale_x_continuous(limits=c(-0.9,0.9)) + 
  #Set the fill of the points manually
  scale_fill_manual(
    values = c("White", "White", "Black", "Black"),
    labels = c("Control Gut", "Control Whole body", 
               "Chlorpyrifos-selected Gut", "Chlorpyrifos-selected Whole body")) +
  #Fit the legend into 2 columns, indicate the legend shape and fill
  guides(fill = guide_legend(override.aes = list(shape = c(22, 21, 22, 21)), ncol = 2,
                             shape = guide_legend(override.aes = list(fill = c("White", "White", "Black", "Black"))))) + 
  #Set x-label to 'NMDS1'
  xlab(expression("NMDS1")) +
  #Set y-label to 'NMDS3'
  ylab(expression("NMDS3")) +
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set position the legend on the top 
    legend.position = "top", 
    #Remove the legend title
    legend.title = element_blank(),
    #Set font size and margins of the legend
    legend.text = element_text(size = 15, margin = margin(2, 2, 10, 2,"pt")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(30, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt")),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_nmdsBray3b

#Make an NMDS plot of NMDS2 and NMDS3 - NOT used in manuscript
Plot_nmdsBray3c <- ggplot(nmdsBray3, aes(x = NMDS2.Avg, y = NMDS3.Avg)) +
  #Set the horizontal error bars based on the average and the 95% confidence interval
  geom_errorbarh(aes(xmin = (NMDS2.Avg-NMDS2.Conf), xmax = (NMDS2.Avg+NMDS2.Conf)), height = 0) +
  #Set the vertical error bars based on the average and the 95% confidence interval
  geom_errorbar(aes(ymin = (NMDS3.Avg-NMDS3.Conf), ymax = (NMDS3.Avg+NMDS3.Conf)), width = 0) +
  #Set the shape of the points manually
  geom_point(aes(fill = sampletypetreatment), shape = c(22, 21, 22, 21), size = 5) +
  #Set the same scale on the x-axis for both Unifrac plots
  scale_x_continuous(limits=c(-0.9,0.9)) + 
  #Set the fill of the points manually
  scale_fill_manual(
    values = c("White", "White", "Black", "Black"),
    labels = c("Control Gut", "Control Whole body", 
               "Chlorpyrifos-selected Gut", "Chlorpyrifos-selected Whole body")) +
  #Fit the legend into 2 columns, indicate the legend shape and fill
  guides(fill = guide_legend(override.aes = list(shape = c(22, 21, 22, 21)), ncol = 2,
                             shape = guide_legend(override.aes = list(fill = c("White", "White", "Black", "Black"))))) + 
  #Set x-label to 'NMDS2'
  xlab(expression("NMDS2")) +
  #Set y-label to 'NMDS3'
  ylab(expression("NMDS3")) +
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set position the legend on the top 
    legend.position = "top", 
    #Remove the legend title
    legend.title = element_blank(),
    #Set font size and margins of the legend
    legend.text = element_text(size = 15, margin = margin(2, 2, 10, 2,"pt")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(30, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt")),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_nmdsBray3c


####### BETA diversity - NMDS ORDINATION - Weighted Unifrac ####### 
#https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html

physeq = physeqFilteredRF
metadf <- data.frame(sample_data(physeq))
unifrac.dist <- UniFrac(physeq, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)
#Set seed to get the same outcome of the permutations
set.seed(10000)
#PERMANOVA using adonis
adonisUnifrac <- adonis(unifrac.dist ~ treatment * sampletype, data = metadf, permutations=10000)
adonisUnifrac
#Homogeneity assumption
physeq.disper <- betadisper(unifrac.dist, metadf$treatment)
permutest(physeq.disper, pairwise = TRUE) 
#P = 0.28 --> ok; assumption for adonis is met for treatment (homogeneous dispersion)
physeq.disper <- betadisper(unifrac.dist, metadf$sampletype)
permutest(physeq.disper, pairwise = TRUE) 
#P = 0.001 --> not ok; differences in composition within sampletypes (heterogeneous dispersion)

#Anderson MJ, Walsh DCI. PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: What null hypothesis are you testing? 
#Ecological monographs [Internet] 2013; 83: 557. Available from: http://doi.org/10.1890/12-2010.1
#PERMANOVA (which is basically adonis()) was found to be largely unaffected by heterogeneity in Anderson & Walsh's simulations but only for balanced designs.
#For unbalanced designs PERMANOVA was too liberal if the smaller group had greater dispersion, and too conservative if the larger group had greater dispersion.

##Get NDMS scores (Stress < 0.20)
nmdsUnifrac <- ordinate(physeq, method="NMDS", distance="unifrac", weighted=TRUE) # Stress = 0.05070283 --> ok
goodness(nmdsUnifrac); stressplot(nmdsUnifrac) 
#Extract scores manually and make your own plot
nmdsUni_scores_test<-data.frame(scores(nmdsUnifrac))
#Add sample ID
nmdsUni_scores_test$Sample<-rownames(nmdsUni_scores_test)
#Strip out the sample data 
sample_vegan <- as(sample_data(physeq),"data.frame") 
sample_vegan$Sample<-rownames(sample_vegan)
#Add in metadata
nmdsUni_scores_treatments<-left_join(nmdsUni_scores_test,sample_vegan,"Sample")
#Calculate the average and the 95% confidence interval of NMDS1/NMDS2/NMDS3 and add values in a data frame
#do.call(data.frame, ...) adds the values in a data frame with each calculation as a single column (instead of as attribute)
#cbind within aggregate lets you run the calculations on multiple variables
#~sampletypetreatment within aggregate is the grouping variable
#by using function(x) in the FUN argument of aggregate you can run multiple calculations (e.g. Avg and Conf)
nmdsUni <- do.call(data.frame, aggregate(cbind(NMDS1,NMDS2) ~ sampletypetreatment, 
                      data = nmdsUni_scores_treatments, 
                      FUN= function(x) c(Avg=round(mean(x), digits=2), 
                                         Conf=round(1.96*(sd(x)/sqrt(length(x))), digits=2))))

#Make an NMDS plot
Plot_nmdsUnifrac <- ggplot(nmdsUni, aes(x = NMDS1.Avg, y = NMDS2.Avg)) +
  #Set the horizontal error bars based on the average and the 95% confidence interval
  geom_errorbarh(aes(xmin = (NMDS1.Avg-NMDS1.Conf), xmax = (NMDS1.Avg+NMDS1.Conf)), height = 0) +
  #Set the vertical error bars based on the average and the 95% confidence interval
  geom_errorbar(aes(ymin = (NMDS2.Avg-NMDS2.Conf), ymax = (NMDS2.Avg+NMDS2.Conf)), width = 0) +
  #Set the shape of the points manually
  geom_point(aes(fill = sampletypetreatment), shape = c(22, 21, 22, 21), size = 5) +
  #Set the same scale on the x-axis for both Unifrac plots
  scale_x_continuous(limits=c(-0.4,0.4)) +
  #Set the fill of the points manually
  scale_fill_manual(
    values = c("White", "White", "Black", "Black"),
    labels = c("Control Gut", "Control Whole body", "Chlorpyrifos-selected Gut", 
               "Chlorpyrifos-selected Whole body")) +
  #Fit the legend into 2 columns, indicate the legend shape and fill
  guides(fill = guide_legend(override.aes = list(shape = c(22, 21, 22, 21)), ncol = 2,
         shape = guide_legend(override.aes = list(fill = c("White", "White", "Black", "Black"))))) + 
  #Set x-label to 'NMDS1'
  xlab("NMDS1") +
  #Set y-label to 'NMDS2'
  ylab(expression("NMDS2")) +
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set position the legend on the top 
    legend.position = "top", 
    #Remove the legend title
    legend.title = element_blank(),
    #Set font size and margins of the legend
    legend.text = element_text(size = 15, margin = margin(2, 2, 10, 2,"pt")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(30, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt")),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_nmdsUnifrac


####### ALPHA diversity - Shannon index and Faith's phylogenetic index ####### 

#Determine Shannon index 
Shannon<-estimate_richness(physeqFilteredRF,measures=c("Shannon"))
?estimate_richness

#Determine Faith's phylogenetic index
OTUphyseqGut=as.data.frame(physeqFilteredRF@otu_table)
TREEphyseqGut=physeqFilteredRF@phy_tree 
TREEphyseqGut #Rooted tree
FaithPD=pd(t(OTUphyseqGut), TREEphyseqGut,include.root=T)

#Add Richness Onto Our Metadata
Shannon$sampleID<-rownames(Shannon)
FaithPD$sampleID<-rownames(FaithPD)
alphaRF=sample_data(physeqFilteredRF)
alphaRF$sampleID<-rownames(alphaRF)
alphaRF<-left_join(Shannon,alphaRF,"sampleID")
alphaRF<-left_join(FaithPD,alphaRF,"sampleID")
alphaRF

#Statistics Shannon index
set_sum_contrasts() 
lmShannon<-lmer(Shannon ~ treatment * sampletype + (1|treatment:population), data=alphaRF)
Anova(lmShannon, type=3)
ranova(lmShannon)
# r.squaredGLMM(lmShannon)
#Marginal r2 (R2m) is our % variance explained just due to the fixed effects. 
#Conditional r2 (R2c) is variance explained by both fixed and random effects. 
shapiro.test(residuals(lmShannon))
hist(resid(lmShannon)) 
qqnorm(resid(lmShannon))
qqline(resid(lmShannon))
leveneTest(Shannon~ treatment * sampletype, data = alphaRF)
aggregate(Shannon~ treatment * sampletype, data = alphaRF, var) 

SEplotShannon <- summary(emmeans(lmShannon, ~ treatment * sampletype, type = "response"))
#Plot shannon diversity
Plot_Shannon <- ggplot(SEplotShannon, aes(x = treatment, y = emmean, group = sampletype)) +
  #Set white/black instead of red/green as color
  scale_fill_manual(values = c("white", "black")) +
  #Set the error bars based on the average and standard error
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  #Set the shape of the points manually
  geom_point(aes(fill = treatment), shape = c(22, 22, 21, 21), size = 5, 
             position = position_dodge(.5), color="black") +
  #Set x-label to 'Selection treatment'
  xlab("Selection treatment") +
  #Set y-label to 'Shannon index'
  ylab(expression("Shannon index")) +
  #Set x-ticks labels
  scale_x_discrete(labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")) +
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set the legend position
    legend.position = "none",
    #Remove the legend title
    legend.title = element_blank(),
    # (labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(25, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt")),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_Shannon


#Statistics Faith's phylogenetic index
set_sum_contrasts() 
lmFaithsPD<-lmer(PD ~ treatment * sampletype + (1|treatment:population), data=alphaRF)
Anova(lmFaithsPD, type=3)
ranova(lmFaithsPD)
# r.squaredGLMM(lmFaithsPD)
#Marginal r2 (R2m) is our % variance explained just due to the fixed effects. 
#Conditional r2 (R2c) is variance explained by both fixed and random effects. 
shapiro.test(residuals(lmFaithsPD))
hist(resid(lmFaithsPD)) 
qqnorm(resid(lmFaithsPD))
qqline(resid(lmFaithsPD))
leveneTest(PD~ treatment * sampletype, data = alphaRF)
aggregate(PD~ treatment * sampletype, data = alphaRF, var) 

SEplotPD <- summary(emmeans(lmFaithsPD, ~ treatment * sampletype, type = "response"))
Plot_PD <- ggplot(SEplotPD, aes(x = treatment, y = emmean, group = sampletype, col=sampletype, shape=sampletype)) +
  #Set white/black instead of red/green as color
  scale_fill_manual(values = c("white", "black")) +
  #Set the error bars based on the average and standard error
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0, 
                position = position_dodge(.5), color="black") +
  #Set the shape of the points manually
  geom_point(aes(fill = treatment), shape = c(22, 22, 21, 21), size = 5, 
             position = position_dodge(.5), color="black") +
  #Set x-label to 'Selection treatment'
  xlab("Selection treatment") +
  #Set y-label to 'Faith's phylogenetic index'
  ylab(expression("Faith's phylogenetic index")) +
  #Set x-ticks labels
  scale_x_discrete(labels = c("control"="Control", "pesticide"="Chlorpyrifos-selected")) +
  #Set the classic theme (no gray panel background and no major/minor grids)
  theme_classic() +
  #Make further adjustments to the lay-out
  theme(
    #Set the legend position
    legend.position = "none",
    #Remove the legend title
    legend.title=element_blank(),
    #Set font size, color and margins of the x-axis tick labels
    axis.text.x = element_text(size = 15, margin = margin(10, 2, 2, 2,"pt"), colour = "black"),
    #Set font size, color and margins of the y-axis tick labels
    axis.text.y = element_text(size = 15, margin = margin(2, 10, 2, 2,"pt"), colour = "black"),
    #Set font size and margins of the x axis label
    axis.title.x = element_text(size = 20, margin = margin(25, 2, 2, 2, "pt")),
    #Set font size and margins of the y axis label
    axis.title.y = element_text(size = 20, margin = margin(2, 30, 2, 2, "pt")),
    #Set plot margins
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
Plot_PD


#Combine the four ggplots in one figure
#http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/

svg(filename = "Figure2_Manuscript.svg", width = 16, height = 21)
Figure2 <- ggarrange(Plot_Shannon, Plot_PD, Plot_nmdsBray3, Plot_nmdsUnifrac, 
                     legend = "none", ncol = 2, nrow = 3, labels = c("A","B","C","D"))
Figure2
dev.off()
