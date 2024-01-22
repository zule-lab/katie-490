###Kathleen Philp BIOL 490 

##Statistics Analysis using GHLT (region not included as interaction)

Philp490Trees <- read.csv("input/MasterDataSheetNDGTrees_DBH5_28-1-2021.csv")
head(Philp490Trees)
dim(Philp490Trees)
str(Philp490Trees)
View(Philp490Trees)

#### Comparing native tree ABUNDANCES across land-use types

library(plyr)

###SLL

#Number of trees native to SLL per land-use type
SLLData.Nat. <- subset(Philp490Trees, Native_SLL == "Y", select = c(Site_Type, Site_Code, Native_SLL))
View(SLLData.Nat.)
SLL.Treespersite <- count(SLLData.Nat., c("Site_Type", "Site_Code"))
View(SLL.Treespersite)
class(SLL.Treespersite)

#Adding column of total trees per site
TOT.treespersite <- count(Philp490Trees, c("Site_Code"))
View(TOT.treespersite)
SLL.Treespersite2 <- merge(SLL.Treespersite, TOT.treespersite, by = "Site_Code")
View(SLL.Treespersite2)
names(SLL.Treespersite2)[names(SLL.Treespersite2)=="freq.x"] <- "Abundance_Native"
names(SLL.Treespersite2)[names(SLL.Treespersite2)=="freq.y"] <- "Abundance_Total"
View(SLL.Treespersite2)


#Adding column of the proportion of native trees per site
SLL.Treespersite2[,"Prop_Native"] <- SLL.Treespersite2$Abundance_Native/SLL.Treespersite2$Abundance_Total

#Adding a column for region
SLL.Treespersite2$Region <- "SLL"
View(SLL.Treespersite2)

###ETF

#Number of trees native to ETF per land-use type
ETFData.Nat. <- subset(Philp490Trees, Native_ETF == "Y", select = c(Site_Type, Site_Code, Native_ETF))
View(ETFData.Nat.)
ETF.Treespersite <- count(ETFData.Nat., c("Site_Type", "Site_Code"))
View(ETF.Treespersite)

#Adding column of total trees per site
TOT.treespersite <- count(Philp490Trees, c("Site_Code"))
View(TOT.treespersite)
ETF.Treespersite2 <- merge(ETF.Treespersite, TOT.treespersite, by = "Site_Code")
View(ETF.Treespersite2)
names(ETF.Treespersite2)[names(ETF.Treespersite2)=="freq.x"] <- "Abundance_Native"
names(ETF.Treespersite2)[names(ETF.Treespersite2)=="freq.y"] <- "Abundance_Total"
View(ETF.Treespersite2)


#Adding column of the proportion of native trees per site
ETF.Treespersite2[,"Prop_Native"] <- ETF.Treespersite2$Abundance_Native/ETF.Treespersite2$Abundance_Total

#Adding a column for region
ETF.Treespersite2$Region <- "ETF"
View(ETF.Treespersite2)

#Combining SLL & ETF

Data.Trees <- rbind(SLL.Treespersite2, ETF.Treespersite2)

#Making sure variables as.factor
Data.Trees$Site_Type <- as.factor(Data.Trees$Site_Type)
Data.Trees$Region <- as.factor(Data.Trees$Region)

#Running the binomial GLM

glm.Data.Trees <- glm(Prop_Native~Site_Type*Region, binomial, weights = Abundance_Total,
                      data = Data.Trees)
install.packages("car")
library(car)
Anova(glm.Data.Trees)

#Pairwise comparisons
install.packages("emmeans")
library(emmeans)

lsmeans(glm.Data.Trees, pairwise~Site_Type|Region, adjust = "tukey")


###****************Clear environment***************

#####Comparing native SPECIES RICHNESS across the four land-use types

Philp490Trees <- read.csv("input/MasterDataSheetNDGTrees_DBH5_28-1-2021.csv")

#Total number of species across sites
TOT.sppspersite <- count(Philp490Trees, c("Site_Code","Species_Code"))
View(TOT.sppspersite)
library(tidyr)
wide.TOT.sppspersite <- spread(TOT.sppspersite, Species_Code, freq)
View(wide.TOT.sppspersite)
wide.TOT.sppspersite[is.na(wide.TOT.sppspersite)] <- 0
TOTrichness <- ddply(wide.TOT.sppspersite, ~Site_Code, function(x){data.frame((RICHNESS=sum(x[-(1:2)]>0)))})
View(TOTrichness)
names(TOTrichness)[names(TOTrichness)=="X.RICHNESS...sum.x...1.2.....0.."] <- "Species_Total"
View(TOTrichness)


###SLL

#Number of species native to SLL across sites
SLLData.Nat.spp <- subset(Philp490Trees, Native_SLL == "Y", select = c(Site_Type, Site_Code, Species_Code))
View(SLLData.Nat.spp)
SLL.sppspersite <- count(SLLData.Nat.spp, c("Site_Type", "Site_Code","Species_Code"))
View(SLL.sppspersite)

wide.SLL.sppspersite <- spread(SLL.sppspersite, Species_Code, freq)
View(wide.SLL.sppspersite)
wide.SLL.sppspersite[is.na(wide.SLL.sppspersite)] <- 0
SLLrichness <- ddply(wide.SLL.sppspersite, ~Site_Type+Site_Code, function(x){data.frame((RICHNESS=sum(x[-(1:2)]>0)))})
View(SLLrichness)
names(SLLrichness)[names(SLLrichness)=="X.RICHNESS...sum.x...1.2.....0.."] <- "Species_Native"
View(SLLrichness)

#Merging spp per site and the total richness
SLLrichness2 <- merge(SLLrichness, TOTrichness, by = "Site_Code")
View(SLLrichness2)

#Adding column of the proportion of native species per site
SLLrichness2[,"Prop_Native"] <- SLLrichness2$Species_Native/SLLrichness2$Species_Total

#Adding a column for region
SLLrichness2$Region <- "SLL"
View(SLLrichness2)


####ETF

#Number of species native to ETF per land-use type
ETFData.Nat.spp <- subset(Philp490Trees, Native_ETF == "Y", select = c(Site_Type, Site_Code, Species_Code))
View(ETFData.Nat.spp)
ETF.sppspersite <- count(ETFData.Nat.spp, c("Site_Type", "Site_Code","Species_Code"))
View(ETF.sppspersite)

wide.ETF.sppspersite <- spread(ETF.sppspersite, Species_Code, freq)
View(wide.ETF.sppspersite)
wide.ETF.sppspersite[is.na(wide.ETF.sppspersite)] <- 0
ETFrichness <- ddply(wide.ETF.sppspersite, ~Site_Type+Site_Code, function(x){data.frame((RICHNESS=sum(x[-(1:2)]>0)))})
View(ETFrichness)
names(ETFrichness)[names(ETFrichness)=="X.RICHNESS...sum.x...1.2.....0.."] <- "Species_Native"
View(ETFrichness)

#Merging spp per site and the total richness
ETFrichness2 <- merge(ETFrichness, TOTrichness, by = "Site_Code")
View(ETFrichness2)

#Adding column of the proportion of native species per site
ETFrichness2[,"Prop_Native"] <- ETFrichness2$Species_Native/ETFrichness2$Species_Total

#Adding a column for region
ETFrichness2$Region <- "ETF"
View(ETFrichness2)


#Combining SLL & ETF

Data.Species <- rbind(SLLrichness2, ETFrichness2)
View(Data.Species)

#Making sure variables as.factor
Data.Species$Site_Type <- as.factor(Data.Species$Site_Type)
Data.Species$Region <- as.factor(Data.Species$Region)

#Running the binomial GLM

glm.Data.Species <- glm(Prop_Native~Site_Type*Region, binomial, weights = Species_Total,
                      data = Data.Species)
library(car)
Anova(glm.Data.Species)

#Pairwise comparisons
install.packages("emmeans")
library(emmeans)

lsmeans(glm.Data.Species, pairwise~Site_Type|Region, adjust = "tukey")

