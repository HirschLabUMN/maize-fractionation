#!/usr/bin/Rscript
library(DBI)
library(RSQLite)
library(tidyverse)
library(magrittr)
##ABOUT: This script is for testing different exploratory plots

# perl make--gene-matrix.pl -d /home/hirschc1/shared/projects/fractionation/cache/rnaseq/ -o out.test
# lete from normCounts;
# .separator "\t"
# .import /home/hirschc1/shared/projects/fractionation/munge/out.test normCounts

setwd('/home/hirschc1/shared/projects/fractionation/cache/rnaseq')
con = dbConnect(SQLite(), dbname="/home/hirschc1/shared/projects/fractionation/data/sqlite/all_info.db")

#########################################################################
# Before we can do anything we need to get metadata and expression info #
#########################################################################

# Get metadata on assignments
assignMeta <- dbReadTable(con, "assignments")
gff <- dbReadTable(con, "gff")
countData <- dbReadTable(con, "normCounts",header=TRUE)
countData <- countData[2:nrow(countData),]
fusions <- read.delim("/home/hirschc1/shared/projects/fractionation/cache/coge/fusions.txt",header=F)


###############################################################################
#                             ON/OFF Analysis                                 #
###############################################################################


retDupMeta <- assignMeta %>%
	mutate( syntelog = case_when(
	  grepl("Zm00001d", cognate) == TRUE ~ paste(cognate,gene,sep="_"),
	  grepl("Zm00008a", cognate) == TRUE ~ paste(gene,cognate,sep="_"))) %>%
	filter(grepl("Zm",syntelog),grepl("Zm",gene)) %>%
  select(gene,syntelog,duplicate,cognate,reciprocal,subgenome,genotype)

# Ignore fusion lines
retDupMeta <- anti_join(retDupMeta,fusions,by=c("gene"="V1"))

countMeta <- inner_join(countData,retDupMeta, by=c("gene"="syntelog"))

countMeta %<>%
	rename(id = gene, gene = gene.y) %>%
	filter(
    	grepl("Zm", duplicate),
   		grepl("Zm", cognate),
   	 	grepl("Zm", reciprocal)
  	) %>%
	select(-duplicate,-cognate,-reciprocal)


# Get the ids of syntelogs in which one of the homologs is not chormosomal so I can avoid these
chrORscaffolds <- inner_join(countMeta,gff, by=c("gene"="gene")) %>%
	filter( chr != 1 & chr != 2 & chr != 3 & chr != 4 & chr != 5 & chr != 6 & chr != 7 & chr != 8 & chr != 9 & chr != 10 ) %>%
	select(id)

countMeta <- anti_join(countMeta, chrORscaffolds, by=c("id"="id"))

# Separate the counts for time being
bcounts <- countMeta %>%
	filter(genotype == "B73") %>%
	select(1:5,7,14:16) %>%
	rename(bd = b73Bd, cp = b73Cp, gk = b73Gk, rt = b73Rt, st = b73St) %>%
	rowwise() %>%
	mutate(count = sum( bd> 3, cp> 3, gk> 3, rt> 3, st> 3 ))

pcounts <- countMeta %>%
	filter(genotype == "PH207") %>%
	select(1,8:11,13:16) %>%
	rename(bd = ph207Bd, cp = ph207Cp, gk = ph207Gk, rt = ph207Rt, st = ph207St) %>%
	rowwise() %>%
	mutate(count = sum( bd> 3, cp> 3, gk> 3, rt> 3, st> 3 ))

# Now combine together again
bpcounts <- rbind(bcounts,pcounts) %>% as.data.frame()
# Put things in proper class
for(i in 2:6) {
		bpcounts[,i] <- as.numeric(bpcounts[,i])
}
for(i in 8:10) {
		bpcounts[,i] <- as.factor(bpcounts[,i])
}

# How many maize1 syntelogs and maize2 syntelog pairs have expression from B73 or PH207 ??
bpcounts %>%
	filter(count != 0) %>%  # here count is the number of tissues >3 normalized counts
	group_by(id) %>%
	filter(row_number()==1) %>%
	group_by(subgenome) %>%
	tally()

# 	subgenome     n
# 1    maize1  4580 / 4834 # comment out the filter line to get denom
# 2    maize2  4491 / 4807
# 					=  9071 / 9641

# Get the number of tissues that B73 and PH207 genes are expressed in (above 3 counts).
bpcounts_by_tissue <- bpcounts %>%
	select(-gene, -bd, -cp, -gk, -rt, -st, -subgenome) %>%
	gather(count, value, -id, -genotype ) %>%
	spread(genotype,value) %>%
	mutate(B73 = factor(B73), PH207=factor(PH207)) %>%
	rename(B73count = B73, PH207count = PH207) %>%
	select(-count)

# Now let's go back and treat each syntelog and tissue combination as separate datapoints.
bpcounts_long <- bpcounts %>%
	select(-gene,-count) %>%
	gather(tissue, value, -id, -genotype, -subgenome) %>%
	spread(genotype,value) %>%
	inner_join(bpcounts_by_tissue, by=c("id"="id"))

ELSE <- TRUE  # this must be set for dplyr if else statement to run.
CUTOFF <- 8.485281  # this value equates to a log2FC > 1.5 compared to 3.

summarized_bpcounts <- bpcounts_long %>%
	filter( (B73 > CUTOFF & PH207 <= 3) | (PH207 > CUTOFF & B73 <= 3) ) %>%
	mutate(B73expr = case_when(
			B73 > PH207 ~ 1,
			ELSE ~ 0),
		   PH207expr = case_when(
			PH207 > B73 ~ 1,
			ELSE ~ 0)) %>%
	select(-B73, -PH207)

final_bpcounts <- summarized_bpcounts %>%
	group_by(id) %>%
	mutate(B73on = sum(B73expr), PH207on = sum(PH207expr)) %>%
	select(-B73expr, -PH207expr, -tissue) %>%
	filter(row_number()==1) %>%
	as.data.frame()

# OK, now let's get some summary numbers
final_bpcounts <- final_bpcounts %<>%
	mutate(direction= case_when(
		B73on == 0 ~ "PH207only",
		PH207on == 0 ~ "B73only",
		ELSE ~ "variable"))

final_bpcounts %>%
	group_by(direction) %>%
	tally()
#   direction     n
# 1   B73only  604
# 2 PH207only  760
# 3  variable  97

# Everywhere where expressed, expressed log2FC > 1.5 than the opposite genotype
final_bpcounts %>%
	filter(direction=="PH207only" & B73count==0) %>%
	group_by(subgenome) %>%
	tally()
#   subgenome     n
# 1    maize1   93
# 2    maize2   91
#      total   184

final_bpcounts %>%
	filter(direction=="B73only" & PH207count==0) %>%
	group_by(subgenome) %>%
	tally()
#   subgenome     n
# 1    maize1    71
# 2    maize2    88
#      total    159

# Total maize1: 93 + 91 = 184
# Total maize2: 71 + 88 = 159
# Total ON/OFF:         = 343
# Proportion ON/OFF: 343 / 9071 = 0.038 or ~3.8%

# Now, we want to see how many of the ON/OFF genes are false positives
# using the maize expression atlas to check
atlas <- read.delim(file="/home/maize/shared/resources/rna-seq/Final_Gene_Atlas_AGPv2_FPKM_matrix_Gene_Nov2014_v4converted.txt",header=F)

atlas %<>%
	select(1,6:84) %>%
	gather(Tissue,FPKM,-V1) %>%
	group_by(V1) %>%
	summarise(mean=mean(FPKM),aboveOne=sum(FPKM > 1))

# How many B73 genes are expressed in these tissues ?
B73_off_genes <- final_bpcounts %>%
	filter(direction=="PH207only" & B73count==0) %>%
	separate(id, c("b73","ph207"), "_", remove=FALSE) %>%
	inner_join(atlas, by=c("b73"="V1")) %>%
	arrange(as.numeric(mean))

# 158 of 184 gene models converted
# 45 have average FPKM across tissues less than 1
# 15 are never expressed greater than FPKM of 1 in any of the tissues

# The following is for general syntenic vs. non-syntenic analyses and was not included in final manuscript:

# Get RPKM values
bRPKM <- dbReadTable(con, "b73RPKM")
pRPKM <- dbReadTable(con, "ph207RPKM")
# Convert character columns to numeric
for(i in 2:7) {
		bRPKM[,i] <- as.numeric(as.character(bRPKM[,i]))
}
for(i in 2:7) {
		pRPKM[,i] <- as.numeric(as.character(pRPKM[,i]))
}
# bind
bpRPKM <- rbind(bRPKM,pRPKM)

# Omit seedling tissue, b/c it seems to globally elevated in B73-v-PH207
bpRPKM %<>% select(-sd)

#########################################################################
# Now we need to get the number of tissues that a genes is expressed in #
#########################################################################

# This is the tissue count if we only consider RPKM as 0 as non-expressed
tissuesExpr <- bpRPKM %>%
  gather("Tissue", "RPKM",2:6) %>%
  filter( RPKM > 0 ) %>%
  group_by(gene) %>%
  count() %>%
  rename(count=n) %>%
  as.data.frame()
# Now we need to get those genes not expressed in any tissue
notissuesExpr <- anti_join(bpRPKM,tissuesExpr,by=c("gene"="gene")) %>%
  select(gene) %>%
  mutate(count=0) %>%
  as.data.frame()
# Combine genes expressed in 1 or more tissues with those not expressed
tissueCount <- rbind(tissuesExpr,notissuesExpr) %>% as.data.frame()

####################################################################################
# Now we need to catergize these tissue counts into: maize1, maize2, non-syntenic  #
####################################################################################

bpRPKMwMeta <- inner_join(bpRPKM,assignMeta, by=c("gene"="gene")) %>%
	select(1:6,subgenome,genotype) %>%
	group_by(gene) %>%
	filter(row_number() == 1)

# First, let's look at syntenic genes
syntenicDat <- inner_join(bpRPKMwMeta,tissueCount,by=c("gene"="gene")) %>%
  select(gene,subgenome,genotype,count) %>%
	as.data.frame()

# Same thing for non-syntenic genes
nonsyntenicDat <- anti_join(tissueCount,bpRPKMwMeta,by=c("gene"="gene"))

# Ok, now let's do some sanity checking
syntenicDat %>%
	group_by(subgenome,genotype) %>%
	tally()
# subgenome genotype     n
# 1    maize1      B73 14589
# 2    maize1    PH207 13937
# 3    maize2      B73  9219
# 4    maize2    PH207  8885

# Yes, these numbers match those in the master file exactly when awk for genes and use uniq

# We have to do a little extra to get the same facet information we have for the syntenic genes
nonsyntenicDat <- nonsyntenicDat %>%
  mutate( genotype = case_when(
	  grepl("Zm00001d", gene) == TRUE ~ "B73",
	  grepl("Zm00008a", gene) == TRUE ~ "PH207"),
	  subgenome="non-syntenic") %>%
  select(gene,subgenome,genotype,count)

# Combine syntenic and non-syntenic back together
allDat <- rbind(syntenicDat,nonsyntenicDat)
#   gene            subgenome   genotype count
# 1 Zm00001d000001    maize1      B73     1

# Now let's add some expression values to this list
allDatwExpr <- inner_join(allDat,bpRPKM,by=c("gene"="gene"))

# Filtering out lines with no expression and plotting global expression distribution
allDatwExpr %>%
  rowwise() %>%
  filter( (bd > 0 | cp > 0 | gk > 0 | rt > 0 | st > 0) ) %>%
  mutate(meanExpr=mean(c(bd,cp,gk,rt,st)),log2RPKM=log2(meanExpr)) %>%
  select(-bd,-cp,-gk,-rt,-st) %>%
  ggplot(aes(log2RPKM,colour=subgenome)) +
    geom_density(aes(linetype=subgenome), size=1.15) +
    facet_grid(genotype~.) +
    theme_classic() +
    scale_colour_manual(values=c("black", "#969696", "#bdbdbd")) +
    scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
    lims(x=c(-10,15)) +
    ggsave(file="/panfs/roc/groups/14/hirschc1/shared/projects/fractionation/graphs/global-expression-facet-notissue.pdf", device="pdf")

######################################################################################
######################################################################################

#################################################################
# Ok, now let's focus on getting some basic expression metrics  #
#################################################################

# Let's start with the proportion of syntenic and nonsyntenic transcription

# Now I'm going to call anything below RPKM value of 1 to be non-expresed.
#allDatwExpr <- inner_join(allDat,bpRPKM,by=c("gene"="Gene")) %>%
#   mutate_each(funs(replace(., . < 1, 0.0)), -gene, -subgenome, -genotype, -count)

non_trans <- allDatwExpr %>%
  filter(count==0) %>%
  group_by(subgenome,genotype) %>%
  tally() %>%
  rename(Non_transcribed=n) %>%
  as.data.frame()

# Do the same for the number of genes transcribed
trans <- allDatwExpr %>%
  filter(count > 0) %>%
  group_by(subgenome,genotype) %>%
  tally() %>%
  rename(Transcribed=n) %>%
  as.data.frame()

# Combine these back together
dat<-cbind(non_trans,trans$Transcribed)
colnames(dat) <- c("Subgenome", "Genotype","Non_Transcribed", "Transcribed")

##################################################
#  Bar-chart of transcribed vs. non-transcribed  #
##################################################

dat %>% gather("Transcription_Status","Count",3:4) %>%
  ggplot(aes(Subgenome,Count,fill=Transcription_Status)) +
    geom_bar(stat="identity") +
    facet_grid(Genotype~.) +
    theme_classic() +
    scale_fill_manual(values=c("darkgrey", "lightgrey")) +
		scale_y_continuous(expand=c(0,0)) +
    coord_flip() +
    ggsave(file="transcription-status.pdf", device="pdf", width = 7.25, height = 3.5, units="in")

# For getting the actual proportion...
dat %>% mutate(Proportion_Transcribed=Transcribed/(Non_Transcribed+Transcribed))

# Subgenome Genotype Non_Transcribed Transcribed Proportion_Transcribed
# 1       maize1      B73             487       14102              0.9666187
# 2       maize1    PH207             487       13450              0.9650570
# 3       maize2      B73             344        8875              0.9626858
# 4       maize2    PH207             351        8534              0.9604952
# 5 non-syntenic      B73            7199        8317              0.5360273
# 6 non-syntenic    PH207            8088        9647              0.5439526

######################################################################################
######################################################################################
######################################################################################

#################################################################
#           Now, let's get global expression trends             #
#################################################################

# Filter the expression data to only consider transcribed genes
# and get the log expression and coefficient of variation
avrExpr <- allDatwExpr %>%
  filter( (bd > 0 | cp > 0 | gk > 0 | rt > 0 | st > 0) ) %>%
  rowwise() %>%
  mutate(avrExpr= mean(c(bd,cp,gk,rt,st)),
         std=sd(c(bd,cp,gk,rt,st)),
         cv=(std/avrExpr),
         logExpr=(log2(avrExpr)))

# Make a box plot of expression values ( although this isn't very informative )
avrExpr %>%
  ggplot(aes(subgenome,logExpr,fill=subgenome)) +
  geom_boxplot() +
  facet_grid(genotype~.) +
  theme_classic() +
  scale_fill_manual(values=c("white","lightgrey","darkgrey")) +
  ggsave(file="global-expression-boxplot.pdf", device="pdf")

# Boxplot isn't informative but what are the actual numbers?
avrExpr %>%
  as.data.frame() %>%
  filter(logExpr>0) %>%
  group_by(genotype,subgenome) %>%
  summarise(meanCV=mean(cv),meanExpr=mean(logExpr),meanSd=mean(std))

# 	genotype    subgenome    meanCV meanExpr   meanSd
# 1      B73       maize1 0.8806316 3.943440 46.60283
# 2      B73       maize2 0.8369492 3.923629 38.58163
# 3      B73 non-syntenic 0.9283515 3.353174 38.29810
# 4    PH207       maize1 0.8066392 4.057550 44.84540
# 5    PH207       maize2 0.7607184 4.031015 38.79319
# 6    PH207 non-syntenic 0.8333154 3.200119 27.47435

# This is to get the proportion for the graph below
exprCounts <- avrExpr %>%
  select(-bd,-cp,-gk,-rt,-st) %>%
  ungroup() %>%
  filter(count != 0 ) %>%
  group_by(genotype,subgenome)

exprCounts %>%
  tally()

# 	genotype    subgenome     n
# 1      B73       maize1 14102
# 2      B73       maize2  8875
# 3      B73 non-syntenic  8317
# 4    PH207       maize1 13450
# 5    PH207       maize2  8534
# 6    PH207 non-syntenic  9647


avrExpr %>%
  group_by(genotype,subgenome,count) %>%
  tally() %>%
  rowwise() %>%
  mutate(denom=case_when(
    (genotype == "B73" & subgenome == "maize1" ) ~ 14102,
    (genotype == "B73" & subgenome == "maize2" ) ~ 8875,
    (genotype == "B73" & subgenome == "non-syntenic" ) ~ 8317,
    (genotype == "PH207" & subgenome == "maize1" ) ~ 13450,
    (genotype == "PH207" & subgenome == "maize2" ) ~ 8534,
    (genotype == "PH207" & subgenome == "non-syntenic" ) ~ 9647),
    prop=n/denom) %>%
  filter(count != 0 ) %>%
  ggplot(aes(subgenome,prop,fill=as.factor(count))) +
    geom_bar(stat="identity",position="stack") +
    facet_grid(~genotype) +
    theme_classic() +
    scale_y_continuous(expand=c(0,0)) +
    coord_flip() +
    scale_fill_manual(values=c("#f7f7f7", "#d9d9d9", "#bdbdbd", "#969696", "#636363", "#252525")) +
    ggsave(file="/panfs/roc/groups/14/hirschc1/shared/projects/fractionation/graphs/tissue-wise-expression_barchart.pdf", device="pdf", width = 7.25, height = 3.5, units="in")

# So, if we consider expression between m1, m2, and non-syntenic genes on the
# the levels of the same # of tissues with expression than m1 is expressed higher
# than m2 and non-syntenic genes. Although the coefficient of variation also follows this
# trend.


######################################################################################
######################################################################################
######################################################################################


# That analysis was for all genes regardless of fractionation class, but now
# let's focus on comparing homologous genes between the two genotypes.
# For now we are only considering cases where *both* maize1 and maize2 sets are retained.

# Get the retained duplicate metadata
retDupMeta <- assignMeta %>%
	mutate( syntelog = case_when(
	  grepl("Zm00001d", cognate) == TRUE ~ paste(cognate,gene,sep="_"),
	  grepl("Zm00008a", cognate) == TRUE ~ paste(gene,cognate,sep="_"))) %>%
	filter(grepl("Zm",syntelog),grepl("Zm",gene)) %>%
  select(gene,syntelog,duplicate,cognate,reciprocal)

# I only want to consider lines in the master file that weren't involved in a duplicate fusion so I need to filter these out.
# This requires some ad hoc filtering outside of R:
# grep "fusion" b73-ph207-synteny-bt2corrected.txt | cut -f 2-5 | sed 's/\t/\n/g' | awk '$1 ~/Zm/' | sort | uniq > fusion-genes.txt
uniqDup <- anti_join(retDupMeta,fusions,by=c("gene"="V1"))
dat <- inner_join(allDatwExpr,uniqDup,by=c("gene"="gene"))

retDupDat <- dat %>%
  filter(
    grepl("Zm", duplicate),
    grepl("Zm", cognate),
    grepl("Zm", reciprocal)
  ) %>%
  select(gene,cognate,syntelog,subgenome)

# tack on the B73 and PH207 expression values
retDupDatwExpr <- inner_join(retDupDat,bRPKM,by=c("gene"="gene")) %>%
	select(-sd) %>%
  rename(b_bd=bd,b_cp=cp,b_gk=gk,b_rt=rt,b_sd=sd,b_st=st) %>%
	inner_join(pRPKM,by=c("cognate"="gene")) %>%
	rename(p_bd=bd,p_cp=cp,p_gk=gk,p_rt=rt,p_sd=sd,p_st=st) %>%
	select(-gene,-cognate,-b_sd,-p_sd)

# Now, we can summarise expression for both B73 and PH207
summarizedDat <- retDupDatwExpr %>%
  rowwise() %>%
  mutate(avrExprB= mean(b_bd,b_cp,b_gk,b_rt,b_st),
         avrExprP= mean(p_bd,p_cp,p_gk,p_rt,p_st),
         logExprB= log2(avrExprB+0.25),
         logExprP= log2(avrExprP+0.25),
         sd_B= sd(c(b_bd,b_cp,b_gk,b_rt,b_st)),
         sd_P= sd(c(p_bd,p_cp,p_gk,p_rt,p_st)),
         cv_B= sd_B/avrExprB,
         cv_P= sd_P/avrExprP)  %>%
  ungroup()

summarizedDat %>%
  filter(avrExprB > 0) %>%
  group_by(subgenome) %>%
  summarise(meancv=mean(cv_B),meanExpr=mean(logExprB))

#   subgenome   meancv meanExpr
# 1    maize1 5.120491 3.391747
# 2    maize2 5.121315 3.120274

summarizedDat %>%
  filter(avrExprP > 0) %>%
  group_by(subgenome) %>%
  summarise(meancv=mean(cv_P),meanExpr=mean(logExprP))

#   subgenome   meancv meanExpr
# 1    maize1 5.004635 3.535159
# 2    maize2 4.669658 3.289134

# When B73 is expressed higher than PH207, what is the average ratio of expression?
summarizedDat %>%
	filter(logExprB > logExprP) %>%
	select(logExprB,logExprP) %>%
	mutate(ratioBP = logExprB/logExprP) %>%
	summarise(avrRatio = mean(ratioBP))

#   avrRatio
# 1 2.334953

# When PH207 is expressed higher than B73, what is the average ratio of expression?
summarizedDat %>%
	filter(logExprP > logExprB) %>%
	select(logExprB,logExprP) %>%
	mutate(ratioBP = logExprP/logExprB) %>%
	summarise(avrRatio = mean(ratioBP))

#   avrRatio
# 1 1.026389
<<<<<<< HEAD
=======

##################################################################################
#                             Homeolog Comparison                                #
##################################################################################

retDupMeta <- assignMeta %>%
	mutate( syntelog = case_when(
	  grepl("Zm00001d", duplicate) == TRUE ~ ifelse(subgenome=="maize1",paste(duplicate,gene,sep="_"),paste(gene,duplicate,sep="_")),
	  grepl("Zm00008a", duplicate) == TRUE ~ ifelse(subgenome=="maize1",paste(duplicate,gene,sep="_"),paste(gene,duplicate,sep="_"))
	)) %>%
	filter(grepl("Zm",syntelog),grepl("Zm",gene)) %>%
  select(gene,syntelog,duplicate,cognate,reciprocal,subgenome)

# I only want to consir lines in the master file that weren't involved in a duplicate fusion so I need to filter these out. This requires some ad hoc filtering outsi of R
#grep "fusion" b73-ph207-synteny-bt2corrected.txt | cut -f 2-5 | sed 's/\t/\n/g' | awk '$1 ~/Zm/' | sort | uniq > fusion-genes.txt
fusions <- read.delim("/home/hirschc1/shared/projects/fractionation/cache/coge/fusions.txt",header=F)
uniqDup <- anti_join(retDupMeta,fusions,by=c("gene"="V1"))
dat <- inner_join(allDatwExpr,uniqDup,by=c("gene"="gene"))

retDupDat <- dat %>%
  select(-gene,-duplicate,-cognate,-reciprocal,-count,-subgenome.y) %>%
	rename(subgenome=subgenome.x) %>%
  gather(tissue, value, -syntelog, -genotype, -subgenome) %>%
  filter(value != "NA") %>%
  spread(subgenome,value) %>%
  group_by(syntelog) %>%
  tbl_df()

hlogDat <- retDupDat %>%
	mutate(m1Expr= log2(maize1+0.25),
				 m2Expr= log2(maize2+0.25)) %>%
	filter(m1Expr > -1.23 | m2Expr > -1.23) %>%
	mutate(ratio= abs(m1Expr/m2Expr),
				 status= case_when(
					 ratio>= 1.5 & m1Expr > m2Expr ~ "m1dom",
					 ratio>= 1.5 & m2Expr > m1Expr ~ "m2dom",
					 ratio<1.5 & m1Expr > m2Expr ~ "m1fav",
					 ratio<1.5 & m2Expr > m1Expr ~ "m2fav",
					 m1Expr==m2Expr ~ "equal")
	) %>%
	select(-maize1, -maize2) %>%
	group_by(status,genotype)

hlogDat %>% group_by(syntelog) %>% filter(row_number()==1) %>% group_by(genotype) %>% tally()
#   genotype     n
# 1      B73  5112
# 2    PH207  4774

hlogDat$status <- factor(hlogDat$status, levels=c("m1fav","m2fav","m1dom","m2dom","equal"))
ggplot(hlogDat,aes(tissue,fill=status)) +
  geom_bar(position="dodge", stat="count") +
  facet_grid(genotype~.) +
  theme_classic() +
  scale_fill_manual(values=c("#ca0020","#f4a582","#0571b0","#92c5","darkgrey")) +
  ggsave(file="retained-duplicate-expression-plot.pdf", vice="pdf")

	retDupDeDat$status <- factor(retDupDeDat$status, levels=c("B73dom","PH207dom","B73fav","PH207fav","Equal"))
	ggplot(retDupDeDat,aes(tissue,fill=status)) +
	  geom_bar(position="dodge", stat="count") +
	  facet_grid(subgenome~.) +
	  theme_classic() +
	  scale_fill_manual(values=c("#ca0020","#f4a582","#0571b0","#92c5","darkgrey")) +
	  ggsave(file="retained-duplicate-expression-plot2.pdf", vice="pdf")

hlogDat %>%
	filter(genotype == "B73" ) %>%
	group_by(status) %>%
	summarise (n = n()) %>%
	mutate(freq = n / sum(n))
# 	status     n       freq
# 1  m1dom  6003 0.21208267
# 2  m1fav 10154 0.35873521
# 3  m2dom  1292 0.04564565
# 4  m2fav 10856 0.38353648


hlogDat %>%
	filter(genotype == "PH207" ) %>%
	group_by(status) %>%
	summarise (n = n()) %>%
	mutate(freq = n / sum(n))
# 	status     n       freq
# 1  m1dom  5638 0.21032605
# 2  m1fav  9765 0.36428412
# 3  m2dom  1100 0.04103559
# 4  m2fav 10303 0.38435425


	mutate(freq=n/sum(n))


hlogDat %>% tally()
# status genotype     n
# 1  m1dom      B73  6003
# 2  m1dom    PH207  5638
# 3  m1fav      B73 10154
# 4  m1fav    PH207  9765
# 5  m2dom      B73  1292
# 6  m2dom    PH207  1100
# 7  m2fav      B73 10856
# 8  m2fav    PH207 10303
# Total == 55,111



	bd=(m1bd/m2bd)-1,
    cp=(m1cp/m2cp)-1,
    gk=(m1gk/m2gk)-1,
    rt=(m1rt/m2rt)-1,
    sd=(m1sd/m2sd)-1,
    st=(m1st/m2st)-1) %>%
  select(id,bd,cp,gk,rt,sd,st) %>%
  gather(Tissue, Value, 2:7) %>%
  mutate(status=case_when(
    (Value >= 1.5) ~ "m1dom",
    (Value <= -0.6) ~ "m2dom",
    (Value < 0 & Value > -1.5) ~ "m2fav",
    (Value > 0 & Value < 1.5) ~ "m1fav",
    (Value == 0) ~ "equal")) %>%)


categories <- c("a","a","b","b","c","c","c","c","d","d","d","d","e","e","e","e")
values <- c(5950,5842,8823,3560,17,0,0,2,159,107,96,57,104,31,46,18)
df <- as.data.frame(cbind(categories,values))

ggplot(df,aes(categories,values,fill=categories)) +
	geom_bar(position="dodge", stat="intity") +
	theme_classic() +
	scale_fill_manual(values=c("#ca0020","#f4a582","#0571b0","#92c5","darkgrey")) +
	ggsave(file="relative-frequency.pdf", vice="pdf")


###############################################################################
#                             ON/OFF Analysis                                 #
###############################################################################


retDupMeta <- assignMeta %>%
	mutate( syntelog = case_when(
	  grepl("Zm00001d", cognate) == TRUE ~ paste(cognate,gene,sep="_"),
	  grepl("Zm00008a", cognate) == TRUE ~ paste(gene,cognate,sep="_"))) %>%
	filter(grepl("Zm",syntelog),grepl("Zm",gene)) %>%
  select(gene,syntelog,duplicate,cognate,reciprocal,subgenome,genotype)

# Ignore fusion lines  
retDupMeta <- anti_join(retDupMeta,fusions,by=c("gene"="V1"))
retDupMeta <- anti_join(retDupMeta,fusions,by=c("cognate"="V1"))

countMeta <- inner_join(countData,retDupMeta, by=c("gene"="syntelog"))

countMeta %<>% 
	rename(id = gene, gene = gene.y) %>% 
	filter(
    	grepl("Zm", duplicate),
   		grepl("Zm", cognate),
   	 	grepl("Zm", reciprocal)
  	) %>% 
	select(-duplicate,-cognate,-reciprocal)


# Get the ids of syntelogs in which one of the homologs is not chormosomal so I can avoid these  
chrORscaffolds <- inner_join(countMeta,gff, by=c("gene"="gene")) %>%
	filter( chr != 1 & chr != 2 & chr != 3 & chr != 4 & chr != 5 & chr != 6 & chr != 7 & chr != 8 & chr != 9 & chr != 10 ) %>%
	select(id)

countMeta <- anti_join(countMeta, chrORscaffolds, by="id")
for(i in 2:13) {
		countMeta[,i] <- as.numeric(countMeta[,i])
}

# Separate the counts for time being
bcounts <- countMeta %>%
	filter(genotype == "B73") %>%
	select(1:7,14:16) %>%
	rename(bd = b73Bd, cp = b73Cp, gk = b73Gk, rt = b73Rt, sd = b73Sd, st = b73St) %>%
	rowwise() %>% 
	mutate(count = sum( bd > 3, cp > 3, gk > 3, rt > 3, sd > 3, st > 3))
pcounts <- countMeta %>%
	filter(genotype == "PH207") %>%
	select(1,8:16) %>%
	rename(bd = ph207Bd, cp = ph207Cp, gk = ph207Gk, rt = ph207Rt, sd = ph207Sd, st = ph207St) %>%
	rowwise() %>% 
	mutate(count = sum( bd > 3, cp > 3, gk > 3, rt > 3, sd > 3, st > 3))

# Now combine together again
bpcounts <- rbind(bcounts,pcounts) %>% as.data.frame()
# Put things in proper class  
for(i in 2:7) {
		bpcounts[,i] <- as.numeric(bpcounts[,i])
}
for(i in 9:11) {
		bpcounts[,i] <- as.factor(bpcounts[,i])
}

# Get the # of tissues that B73 and PH207 genes are expressed in (above 3 counts).
bpcounts_by_tissue <- bpcounts %>%
	select(-gene, -bd, -cp, -gk, -rt, -sd, -st, -subgenome) %>%
	gather(count, value, -id, -genotype ) %>%
	spread(genotype,value) %>%
	rename(B73count = B73, PH207count = PH207) %>%
	select(-count)

# Now let's go back and treat each syntelog and tissue combination as separate datapoints. 
bpcounts_long <- bpcounts %>%
	select(-gene,-count) %>%
	gather(tissue, value, -id, -genotype, -subgenome) %>%
	spread(genotype,value) %>%
	inner_join(bpcounts_by_tissue, by="id") 
#
ELSE <- TRUE
CUTOFF <- 8.485281

summarized_bpcounts <- bpcounts_long %>% 
	filter( (B73 > CUTOFF & PH207 <= 3) | (PH207 > CUTOFF & B73 <= 3) ) %>%
	mutate(B73expr = case_when(
			B73 > PH207 ~ 1, 
			ELSE ~ 0),
		   PH207expr = case_when(
			PH207 > B73 ~ 1,
			ELSE ~ 0)) %>% 
	select(-B73, -PH207) 

# B73/PH207 on - # of times B73 gene is ON and PH207 gene is OFF across tissues
final_bpcounts <- summarized_bpcounts %>%
	group_by(id,subgenome,B73count,PH207count) %>% 
	summarise(B73on = sum(B73expr), PH207on = sum(PH207expr)) %>% 
	as.data.frame()

# OK, now let's get some summary numbers 
final_bpcounts <- final_bpcounts %<>%
	mutate(direction= case_when(
		B73on == 0 ~ "PH207only",
		PH207on == 0 ~ "B73only", 
		ELSE ~ "variable")) %>%
	filter(B73count==0 | PH207count==0)

final_bpcounts %>% 
	group_by(direction) %>% tally()
#   direction     n
# 1   B73only  186
# 2 PH207only  100

head(bpcounts_long) # might be able to do something from here

# For determining the influence of seedling tissue:  
bpcounts_long %>% 
	filter( (B73 > CUTOFF & PH207 <= 3) | (PH207 > CUTOFF & B73 <= 3) )  %>% 
	semi_join( final_bpcounts, by=c("id"="id") ) %>%
	filter(B73count==1 | PH207count==1) %>%
	mutate(category= case_when(
		B73 > PH207 ~ "B73ON", 
		ELSE ~ "PH207ON")) %>% 
	group_by(tissue,category) %>% tally()

#  1     bd    B73ON     2
#  2     cp    B73ON     2
#  3     cp  PH207ON     1
#  4     gk    B73ON    11
#  5     gk  PH207ON     9
#  6     rt    B73ON     2
#  7     sd    B73ON    33
#  8     sd  PH207ON     1
#  9     st    B73ON     2
# 10     st  PH207ON     2

# Summary 
final_bpcounts %>%
	filter(direction=="PH207only") %>%
	group_by(PH207count) %>%
	tally()
#   PH207count     n
# 1          1    13
# 2          2    14
# 3          3    10
# 4          4     9
# 5          5    10
# 6          6    44

	
final_bpcounts %>%
	filter(direction=="B73only") %>%
	group_by(B73count) %>%
	tally()
# # A tibble: 6 x 2
#   B73count     n
# 1        1    52
# 2        2    53
# 3        3    23
# 4        4    19
# 5        5    10
# 6        6    29

# Output for Promoter analysis
onoff_genes <- final_bpcounts %>%
	separate(id, c("b73","ph207"), "_", remove=FALSE) %>% 
	select(subgenome,id,b73,ph207)
# Write to outfile to analyze promoter seqs
onoff_genes %>%
	write.table(file="/home/hirschc1/shared/projects/fractionation/cache/tss/on-off-genelist.txt", sep="\t",row.names=FALSE,quote=FALSE)

# We also want a control set of genes to compare these promoter similarities to. 
# a good control should be expressed in at lest one tissue above the cutoff

# Now get all other maize1 and maize2 genes 
non_onoff_genes <- anti_join(countMeta,onoff_genes, by=c("id"="id")) %>% 
						separate(id, c("b73","ph207"), "_", remove=FALSE) 

# REMEMBER TO FILTER OUT DUPLICATES
non_onoff_genes %>% 
	filter( (b73Bd > CUTOFF | b73Cp > CUTOFF | b73Gk > CUTOFF | b73Rt > CUTOFF | b73Sd > CUTOFF | b73St > CUTOFF) &  
		    (ph207Bd > CUTOFF | ph207Cp > CUTOFF | ph207Gk > CUTOFF | ph207Rt > CUTOFF | ph207Sd > CUTOFF | ph207St > CUTOFF) 
		   ) %>%
	select(subgenome,id,b73,ph207) %>%
	write.table(file="/home/hirschc1/shared/projects/fractionation/cache/tss/control-genelist.txt",sep="\t", row.names=FALSE,quote=FALSE)

