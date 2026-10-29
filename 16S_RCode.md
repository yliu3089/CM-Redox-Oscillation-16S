install.packages("devtools")

library(devtools)



devtools::install\_github("benjjneb/dada2")

library("dada2")

library(usethis)

remove.packages("rlang")

install.packages("rlang")

library(rlang)



remove.packages("vctrs")

install.packages("vctrs")

install.packages("remotes")

remotes::install\_github("adrientaudiere/MiscMetabar")



if (!require("BiocManager", quietly = TRUE))

&#x20; install.packages("BiocManager")



BiocManager::install("phyloseq")



install.packages("remotes")

remotes::install\_github("adrientaudiere/MiscMetabar")

library("MiscMetabar")



install.packages("circlize")

library(circlize)

\# Source to install these are Biocondcutor/CRAN/GitHub

library(BiocManager)

BiocManager::install("microbiome")

library(microbiome) # data analysis and visualisation # Bioconductor

library(phyloseq) # also the basis of data object. Data analysis and visualisation # Bioconductor

library(RColorBrewer) # nice color options # CRAN

library(ggpubr) # publication quality figures, based on ggplot2 # CRAN

library(DT) # interactive tables in html and markdown # CRAN

library(data.table) # alternative to data.frame # CRAN

library(dplyr) # d # CRAN



install.packages("devtools")

devtools::install\_github("microsud/microbiomeutilities")



library(microbiomeutilities) # GitHub

library(tidyr) # CRAN

library(ggplot2)



library(tidyverse)

library(readxl)

library(ape)



ps <- read\_phyloseq(otu.file = "Fill with otu file", 

&#x20;                   taxonomy.file = "Fill with taxonomy file", 

&#x20;                   metadata.file = "Redox Oscillation\_metadata.csv", 

&#x20;                  ) 



ps@tax\_table<-ps@tax\_table\[,-1]

print(ps)



any(taxa\_sums(ps)==0)

ps1 <- prune\_taxa(taxa\_sums(ps) !=0, ps)

print(ps1)

summarize\_phyloseq(ps)

summarize\_phyloseq(ps1)

DT::datatable(tax\_table(ps1))



taxa\_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))



tax\_table(ps1)\[,colnames(tax\_table(ps1))] <- gsub(tax\_table(ps1)\[,colnames(tax\_table(ps1))],pattern="\[a-z]\_\_",replacement="")



tax\_table(ps1)\[,colnames(tax\_table(ps1))] <- gsub(tax\_table(ps1)\[,colnames(tax\_table(ps1))],pattern="=\*",replacement="")



\# remove astriks

tax\_table(ps1)\[,colnames(tax\_table(ps1))] <- gsub(tax\_table(ps1)\[,colnames(tax\_table(ps1))],pattern="\[\*]",replacement="")



\# remove \~

tax\_table(ps1)\[,colnames(tax\_table(ps1))] <- gsub(tax\_table(ps1)\[,colnames(tax\_table(ps1))],pattern="\[\~]",replacement="")



ps1.a <- ps1



tax\_table(ps1.a)\[tax\_table(ps1.a) == "k\_\_<empty>"] <- "k\_\_"

tax\_table(ps1.a)\[tax\_table(ps1.a) == "p\_\_<empty>"] <- "p\_\_"

tax\_table(ps1.a)\[tax\_table(ps1.a) == "c\_\_<empty>"] <- "c\_\_"

tax\_table(ps1.a)\[tax\_table(ps1.a) == "o\_\_<empty>"] <- "o\_\_"

tax\_table(ps1.a)\[tax\_table(ps1.a) == "f\_\_<empty>"] <- "f\_\_"

tax\_table(ps1.a)\[tax\_table(ps1.a) == "g\_\_<empty>"] <- "g\_\_"



tax\_table(ps1.a)\[tax\_table(ps1.a) == "o\_\_Unknown\_Order"] <- "o\_\_"

tax\_table(ps1.a)\[tax\_table(ps1.a) == "c\_\_Unknown\_Class"] <- "c\_\_"

tax\_table(ps1.a)\[tax\_table(ps1.a) == "f\_\_Unknown\_Family"] <- "f\_\_"



\# some more ambiguities 

tax\_table(ps1.a)\[tax\_table(ps1.a) == "o\_\_Unknown\_Order"] <- "o\_\_"

tax\_table(ps1.a)\[tax\_table(ps1.a) == "c\_\_Unknown\_Class"] <- "c\_\_"

tax\_table(ps1.a)\[tax\_table(ps1.a) == "f\_\_Unknown\_Family"] <- "f\_\_"



ps1 <- ps1.a

unique(tax\_table(ps1.a)\[,2])



\# Remove ASVs not classified to any Kingdom/Domain

ps1.as <- subset\_taxa(ps1, Kingdom !="NA")

\# Remove ZP30r2 sample as it has only one read

ps1.as <- prune\_samples(sample\_sums(ps1.as) >50, ps1.as)



\# Convert counts to compositional

ps1.rel <- microbiome::transform(ps1, "compositional")

ps1.rel@tax\_table<-ps1.rel@tax\_table\[,-1]

ps1.phy.Os <- aggregate\_top\_taxa2(ps1.rel, top=10 ,level="Genus")







\###Because the 16S rRNA amplicon sequences in this study were sequenced in different batches, we ran the above code separately for each batch

\###the resulting datasets were then merged into a single phyloseq object (Red\_Oscillation), which includes 16S rRNA amplicon data from each batch, each interval

Red\_Oscillation <- merge\_phyloseq(ps1.phy.Os, ps1.phy.oxic5)

Red\_Oscillation <- aggregate\_top\_taxa2(Red\_Oscillation, top=10 ,level="Genus")

otu\_table(Red\_Oscillation)<-otu\_table(otu\_table\_un,taxa\_are\_rows = TRUE)



\##Create a named vector with the replacements

replacements <- c("Anoxic\_1" = "Anoxic\\nInterval 1",  

&#x20;                 "Oxic\_1" = "Oxic\\nInterval 1",  

&#x20;                 "Anoxic\_2" = "Anoxic\\nInterval 2",

&#x20;                 "Oxic\_2" = "Oxic\\nInterval 2",

&#x20;                 "Oxic\_5\_1" = "Oxic\\nInterval 5-1",

&#x20;                 "Oxic\_5\_2" = "Oxic\\nInterval 5-2",

&#x20;                 "Anoxic\_Control" = "Anoxic\\nControl",

&#x20;                 "Oxic\_Control" = "Oxic\\nControl")



\# Extract current sample names

current\_sample\_names <- sample\_names(Red\_Oscillation)



\# Replace sample names using the named vector

new\_sample\_names <- replacements\[current\_sample\_names]



\# Assign the new sample names back to the phyloseq object

sample\_names(Red\_Oscillation) <- new\_sample\_names



\# Verify the changes

sample\_names(Red\_Oscillation)



\### Or at % scale

Red\_Os\_Plot <- plot\_composition(Red\_Oscillation) + 

&#x20; geom\_bar(stat = "identity") +

&#x20; xlab("Sample")+

&#x20; ylab("Relative abundance \\nat the genus level (%)")+

&#x20; theme(axis.text.x = element\_text(angle = 0))+

&#x20; scale\_fill\_manual("Genus",values=c("Anoxic\_1","Oxic\_1","Anoxic\_2",

&#x20;                                    "Oxic\_2","Oxic\_5\_1","Oxic\_5\_2",

&#x20;                                    "Anoxic\_Control","Oxic\_Control"))+

&#x20; theme(panel.grid.major = element\_blank(), 

&#x20;       panel.grid.minor = element\_blank(),

&#x20;       panel.background = element\_blank(),

&#x20;       legend.text=element\_text(size=2),

&#x20;       axis.title=element\_text(size=26),

&#x20;       axis.text.y = element\_text(size = 30),

&#x20;       axis.text.x = element\_text(size = 18,hjust=0.5),

&#x20;       axis.title.x = element\_blank(),

&#x20;       legend.title= element\_text(size=15))+

&#x20; theme(legend.text=element\_text(size=15)

&#x20; )+

&#x20; scale\_x\_discrete(limits = c("Anoxic\\nInterval 1","Oxic\\nInterval 1","Anoxic\\nInterval 2",

&#x20;                             "Oxic\\nInterval 2","Oxic\\nInterval 5-1","Oxic\\nInterval 5-2",

&#x20;                             "Anoxic\\nControl","Oxic\\nControl"))+

&#x20; scale\_y\_continuous(expand = c(0, 0), limits = c(0, 1.02),

&#x20;                    breaks = scales::pretty\_breaks(n = 4))+

&#x20; guides(fill = guide\_legend(ncol = 1))+

&#x20; scale\_fill\_manual("Genus",values=c("red2","blue3",

&#x20;                                    "khaki1","coral4","palegreen1",

&#x20;                                    "green1","gold",'hotpink1',"lightblue1",

&#x20;                                    "grey3"))





