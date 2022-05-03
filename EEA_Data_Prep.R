library(plater) 
library(readxl)
library(tidyverse)
library(writexl)
library(stringr)
library(broom)
library(naniar)
library(AICcmodavg)
library(gridExtra)
library(ggpubr)
library(here)

# EEA Data Prep - Table of contents 
## 1 - Use Plate Data Prep script to get data from excel sheets and through plater
## 2 - Cleaning EEA Data (Renaming, etc.)
## 3 - Calculate Activities- see details below in that section
## 4 - Data Analysis


## This script relies upon consistent naming conventions for your plates! You can pick whatever order you'd like, but you should have these attributes: 
### 1) Plate reads from each batch should all be stored in 1 excel file, with each plate on different tabs. Different batches should be in different excel files. 
### 2) The name of each tab should be plot #/name, or Bpl (B plate) 
### 3) Batch number/name should be in the name of the excel file, in a consistent space such that it can be pulled out by index #

## For example, my plates are named 2021_PARCE_EEA_Batch7_EL (Thats YEAR_EXPERIMENT_ASSAY_BATCH_INITIALS). You can use any consistent format! 

### We use pipes a lot in here (%>%). If you're unfamiliar, look it up here! https://r4ds.had.co.nz/pipes.html

# 1 - Plate Data Prep 
## Use ReadInPlateData script to get your raw A and B plates read into R (read them in together!)
## https://github.com/emmalink1/ReadInPlateData
## You should now have a dataset with 3 columns: 
  ## Plate, chr type, with plate number on the end of the batch worksheet name
  ## wells, chr type; 
  ## and plate reads, int type, probably named something weird

## Make sure to save data your read-in data at this point so you can go back to it if need be. 







# 2 - Cleaning EEA Data (Renaming, etc.)

## You should now have a dataset with 3 columns: Plate, chr, with plate number on the end of the batch worksheet name; wells, chr; and plate reads, int, probably named something weird 
## Now, we're going to rename that plate read column, the Plate column, and take the batch # and plate # out of the file extension

EEAData2021 <- dplyr::rename(EEAData2021, PlateRead = ...2)

EEAData2021 <- dplyr::rename(EEAData2021, "FileExtension" = "Plate")

## note, I'm able to take the batch # and plate # out this way because my FileExtensionFormat is consistent; plot # is the last 3 digits and batch the last
## If your consistent format is different, you will need to pull different indexes of the string to extract your plot and batch information

### pulls out string of plot number or "Bpl" for blanks plate
EEAData2021$Plot <- str_sub(EEAData2021$FileExtension, -3, -1)
### pulls out string of batch
EEAData2021$Batch <- str_sub(EEAData2021$FileExtension, 21, 21)

##Double check how many plates you ran in each batch (did they all get in?) and whether you still have any duplicates you meant to throw out 
EEAData2021 %>% count(Batch)
PlotList <- EEAData2021 %>% count(Plot)


## Now, we need to define what each well is based on the well map (provided blank template A and B). Change to your relative path names.
###Read in templates of well plates 
EEATemplateA <- read_plate("/Users/emmalink/Documents/R/PARCE/Analysis 2020/EEA_Analysis /blank_template_a.csv")
EEATemplateA <- EEATemplateA %>% 
  add_column(TYPE = "A")

EEATemplateB <- read_plate("/Users/emmalink/Documents/R/PARCE/Analysis 2020/EEA_Analysis /blank_template_b.csv")
EEATemplateB <- EEATemplateB %>% 
  add_column(TYPE = "B")

## Now we split our data, merge both types with relative well plate templates, and rejoin 
### Add blank column where type of plate will go
EEAData2021 <- EEAData2021 %>% 
  add_column(TYPE = "NA")

### I can tell plate type based on whether the plot field has a plot number (a plate), or "Bpl" (my naming convention, yours may be different). This 4 loop iterates over all items in dataset and assigns a type (a or B) based on the "plot" field
for (i in 1:nrow(EEAData2021)){
  if (EEAData2021[i, "Plot"]  == "Bpl") {
    EEAData2021[i, "TYPE"] = "B"
  } else {
    EEAData2021[i, "TYPE"] = "A"
  }
}
  
EEAData2021$TYPE <- as.factor(EEAData2021$TYPE)

#Re-merge keyed datasets
EEAData2021_KeyedA <- merge(EEAData2021, EEATemplateA, by = c("TYPE", "Wells"))
EEAData2021_KeyedB <- merge(EEAData2021, EEATemplateB, by = c("TYPE", "Wells"))
EEAData2021_Keyed <- rbind(EEAData2021_KeyedA, EEAData2021_KeyedB)


## merge with metadata (we need plate read times and weights)
## my metadata file has depth, weight, time slurry add, time NaOH add, time read, moisture fraction, dry weight equivalent, and notes (EEA_PARCE_metadata2021.csv)
EEA_PARCE_Metadata_2021 <- read.csv("/Users/emmalink/Documents/R/PARCE/Analysis 2020/EEA_Analysis /EEA_PARCE_Metadata_2021.csv")
EEAData2021_WMetadata <- merge(EEAData2021_Keyed, EEA_PARCE_Metadata_2021, by = c("Plot", "Batch"), all.x = TRUE)

##save data at this intermediate point in an intermediate data files folder 
setwd("/Users/emmalink/Documents/R/PARCE/Analysis 2020/EEA_Analysis /EEA_IntermediateDataFiles")
write.csv(EEAData2021_WMetadata, file = "2021DataWMetadata_NoCalculationsDone.csv")


## Now you should have a datasheet with columns Plot, Batch, TYPE, Wells, FileExtension, PlateRead, Template, Depth, Weight, time slurry add, time NaOH add, Time Read, Moisture Fraction, dry weight equivalent, and notes 









# 3 - Calculations (Divide up the dataset at this point)

##Now, ready to calculate activities for each plate :) The workflow is to
## 1) Calculate means of reads for each plate, for each template type, and inspect for unreasonable variability 
## 3) Calculate standard curve for A plate 
## 4) Calculate standard curve for B plate 
## 5) Calculate quench coefficient (slope A Plate STDS / slope B plate STDS)  
## 6) Calculate net fluorescence 
## 7) Calculate emission coefficient 
## 8) Calculate activity and convert to desired units 

### We are generally going to use the workflow of group_by factors and then perform a calculation   

##1) calculate Means for each template type (read: column) on each plate 

plateMeans <- EEAData2021_WMetadata %>%
  group_by(TYPE, Plot, Batch, Template) %>%
  mutate(Template_Mean = mean(PlateRead), n = n(), sd = sd(PlateRead), cv = sd(PlateRead)/mean(PlateRead)*100, .groups = "keep") %>%
  ungroup()

## Consider saving plateMeans as an intermediate data for QAQC reporting

## Optional QAQC evaluation:
##This part should be replaced with a function but I didn't get to it yet. Inspect averages for bad CV's (>10) and consider removing data points to improve CV. Note:I generally DONT do this level of data cleanup until after I've done initial data analaysis, if it looks like there are significant trends that a cleaner dataset might show better. If there's a really bad plate I'll rerun it
badAverages <- plateMeans %>%
  filter(Template %in% c("quench_MUB0.16", "quench_MUB1.25", "quench_MUB0.625", "quench_MUB2.5"))
##See below on string lists to look into cv of Buffer, Hombl, etc. 
badAveragesIndex <- which(badAverages$cv > 10)

badAveragesNames <- badAverages[badAveragesIndex, 2]
BadMUB <- unique(badAveragesNames)
view(unique)

##also run above with : 
## filter(Template %in% "Buf")
## filter(Template %in% "Hombl")
## filter(Template %in% c("Cello_sub_blank", "P_sub_blank", "NAG_sub_blank", "BG_sub_blank"))
## filter(Template %in% c("quench_MUB0.16", "quench_MUB1.25", "quench_MUB0.625", "quench_MUB2.5"))

###Note that this type data is really messy. These are the issues that I saw with my 2021 dataset: 
##buffer has outliers in : 105, 107, 114, 116, 201, 202, 203, 208, 209, 210, 213, 214, 215, 216, 302, 305, 306, 308, 313, 314, 401, 406, 408, 410, 411, 413, 415, 416, 503, 504, 505, 507, 508, 515
##Nearly every single plate read has bad cv for quench_MUB (74 of them :/)
##Almost half have high substrate blanks
##Also almost half have high hombl


## 2) Calculate Standard Curve for A plates 
### Note that we use the values in the presence of the homogenate (uncorrected) for this standard curve 

##Here, add quench concentration as a variable, and pull out the quench data into another dataset
APlateSTD <- plateMeans %>%
  filter(Template %in% c("quench_MUB0.16", "quench_MUB0.625", "quench_MUB1.25", "quench_MUB2.5"))

APlateSTD$Concentration <- str_sub(APlateSTD$Template, -3, -1) 
APlateSTD$Concentration <- str_replace(APlateSTD$Concentration, "625", ".625")
APlateSTD$Concentration <- str_replace(APlateSTD$Concentration, "\\.25", "1.25")

APlateSTD$Concentration <- as.numeric(APlateSTD$Concentration)

##Inspect the standard curve visually just to see distrubution across A plates 
ASTDSlopeScatter <- APlateSTD %>% 
  ggplot(aes(x=Concentration, y = Template_Mean)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(color = Plot))+
  theme(legend.position = 'none')

ASTDSlopeScatter

##So I've found that this distribution generally looks terrible...It's good to check if there are any patterns across your batches and whatnot and think about what that might mean, but I think this is just kinda what it is.
##The one plate where all the values are 0 is 316, already knew that one needs to be thrown out

##Use group_by and run model in mutate call
##Here, we individually pull out Intercept, Slope, and R^2 by index in the tidy and glance objects. could go back for anything else  

APlateSTD <- APlateSTD %>%
  group_by(Plot, Batch) %>%
  mutate(AIntercept = tidy(lm(Template_Mean ~ Concentration))[1,2], ASlope = tidy(lm(Template_Mean ~ Concentration))[2,2], AR2 <- glance(lm(Template_Mean ~ Concentration))[1,1])

##For the A Plates we have no hard and fast R^2 cutoff but inspect and consider pulling outliers
which(APlateSTD$r.squared < 0.985)


## 3) Calculate standard curve for B plates 
### Note that we use the values in the presence of the buffer (uncorrected) for this standard curve 

##take out B plate with standards, get standard concentrations out and change to a number 
BPlateSTD <- plateMeans %>%
  filter(Template %in% c("MUB0.625", "MUB0.16", "MUB1.25", "MUB2.5"))

BPlateSTD$Concentration <- str_sub(BPlateSTD$Template, -4, -1)
BPlateSTD$Concentration <- str_replace(BPlateSTD$Concentration, "B2.5", "2.5")
BPlateSTD$Concentration <- as.numeric(BPlateSTD$Concentration)

##Calculate curve and get intercept, slope, and R^2 

BPlateSTD <- BPlateSTD %>%
  group_by(Plot, Batch) %>%
  mutate(BIntercept = tidy(lm(Template_Mean ~ Concentration))[1,2], BSlope = tidy(lm(Template_Mean ~ Concentration))[2,2], BR2 <- glance(lm(Template_Mean ~ Concentration))[1,1])

##Inspect R^2 
which(BPlateSTD$r.squared < 0.985)
##Here, we do have a hard and fast rule for removing outliers or considering re-running.  
##Batch 7, R^2 = 0.983, investigate for outliers 

##Ok, now we have the MUB standard curve from the B plate and from the A plates. 


## 4) Calculate quench coefficient (slope plate A STD / slope plate B STD)
### Since all we need now from the plates are the std curves, we make the merged dataset cleaner by pairing down to those values and the batch
### Pick out batch # and intercept and slope estimates for B plates, and consolidate to only 1 entry per plate (since we used average values, don't need to keep all duplicates)  
BSTDS <- BPlateSTD[,c(3,23,24)]
BSTDS <- distinct(BSTDS)

## For A plates, keep batch, plot, std curve estimates, R^2  
ASTDS <- APlateSTD[,c(2,3, 23, 24)]
### Keep only 1 entry per plot (since we created standard curve off of averaged values)
ASTDS <- distinct(ASTDS)

### Merge plates
ABSTDS <- merge(ASTDS, BSTDS, by = "Batch", all.x = TRUE)

### Calculate Quench Coefficient 
ABSTDS$QuenchCoef <- ABSTDS$ASlope /ABSTDS$BSlope
  

## 5) Calculate net flourescence ((Assay-homogenate Control)/Quench coefficient) - substrate control

## Add Homogenate average as a column (homogenate control)
HomogenateControl <- plateMeans %>%
  group_by(TYPE, Plot, Batch) %>%
  filter(Template %in% "Hombl") %>%
  summarize(HomogenateControl = mean(PlateRead), .groups = "keep")

plateMeans <- merge(HomogenateControl, plateMeans, by = c("TYPE", "Plot", "Batch"))

## Add substrate blank average as a column (Substrate Control)
plateMeans$TemplateShort <- str_sub(plateMeans$Template, 1, 2)

SubstrateControl <- plateMeans %>%
  filter(Template %in% c("Cello_sub_blank", "NAG_sub_blank", "P_sub_blank", "BG_sub_blank")) %>%
  group_by(TYPE, Plot, Batch, Template) %>%
  summarize(SubstrateControl = mean(PlateRead), TemplateShort = TemplateShort, .groups = "keep")

SubstrateControl <- unique(SubstrateControl)
SubstrateControl <- SubstrateControl[,c(2,3,5,6)]

plateMeans <- merge(SubstrateControl, plateMeans, by = c("Plot", "Batch", "TemplateShort"))

## Take unnecessary rows (std curve calcs) out of Plate Means 
PlateMeansJustReads <- plateMeans %>%
  filter(Template %in% c("Cello_assay", "BG_assay", "NAG_assay", "P_assay"), .preserve = TRUE)

## Merge standard curves with plate means 
PlateMeansReadsCurves <- merge(PlateMeansJustReads, ABSTDS)  

## just need the averages of each enzyme type from each plate  
PlateMeansReadsCurves <- PlateMeansReadsCurves[,-c(2,3, 5, 7:10, 19, 24)]
PlateMeansReadsCurves <- distinct(PlateMeansReadsCurves)

PlateMeansReadsCurves$NetFlorescence <- ((PlateMeansReadsCurves$Template_Mean - PlateMeansReadsCurves$HomogenateControl)/PlateMeansReadsCurves$QuenchCoef) -PlateMeansReadsCurves$SubstrateControl

## 7) Calculate the emission coefficient 
### emission coefficient = fluorescence/nmol = BSlope/Assay Volume
### assay volume units** are 250uL = 0.00025 L = 0.250 mL 

PlateMeansReadsCurves$EmissionCoefficient <- PlateMeansReadsCurves$BSlope/ 0.250

## 8) Activity 
## (nmol/g/hr) = (net fluoresence x buffer vol (ml)) / (Emission coefficient x homogenate volume (ml) x time (hr) x soil mass (g))
## buffer volume = 50 ml 
## homogenate volume = 0.2 ml
PlateMeansReadsCurves$Time_Read <- strptime(PlateMeansReadsCurves$Time_Read, format = "%H:%M")
PlateMeansReadsCurves$Time_Slurry_Add <- strptime(PlateMeansReadsCurves$Time_Slurry_Add, format = "%H:%M")

PlateMeansReadsCurves$IncubationTime <- as.numeric(difftime(PlateMeansReadsCurves$Time_Read, PlateMeansReadsCurves$Time_Slurry_Add, units = "hour"))
PlateMeansReadsCurves$Activity <- (PlateMeansReadsCurves$NetFlorescence * 50)/ (PlateMeansReadsCurves$EmissionCoefficient *0.2 * PlateMeansReadsCurves$IncubationTime * PlateMeansReadsCurves$Dry_Weight_Equivalent)

###Get value out of dataframe (which is how it's save, a leftover from the broom/tidy situation)
PlateMeansReadsCurves$Activity <- PlateMeansReadsCurves$Activity[,1]

## Export as an intermediate dataset 
write.csv(PlateMeansReadsCurves, "/Users/emmalink/Documents/R/PARCE/Analysis 2020/EEA_Analysis /EEA_IntermediateDataFiles")







# 4 - Data Analysis
## You should probably do your own, experiment specific thing! But the calculations are done, you know how to do this part :)
## Note that below here, my script is not very clean.  

## Assess histogram of activities 
activityBoxplots <- ggplot(PlateMeansReadsCurves, aes(x=Template, y = Activity)) + 
  geom_boxplot()
activityBoxplots

activityHist <- ggplot(PlateMeansReadsCurves, aes(x = Activity)) + 
  geom_histogram() + 
  facet_grid(~Template)

activityHist

###Values have some really big and small outliers. Identifiy them. Consider removal of outliers in the replicates if they exist, rerunning, or other QAQC

### Put data into final format, we want columns Year, Batch, Plot, BG_Activity, Cello_Activity, NAG_Activity, P_Activity, Climate, and Crop 

FinalDataset <- PlateMeansReadsCurves[,c(1,2,3,30)]
FinalDataset <- pivot_wider(FinalDataset, names_from = Template, values_from = Activity)
FinalDataset$Year <- 2021
FinalDataset$Batch <- as.factor(FinalDataset$Batch)

## Merge with treatment key (use your own!)
FinalDataset <- merge(FinalDataset, Key_Plot_Treatment, by = "Plot")
FinalDataset$Block <- as.factor(FinalDataset$Block)

###We will need to remove 108 (index = 8) and 302 (index = 35) because they are way to high/low, something must have gone wrong
FinalDataset <- FinalDataset[-c(8,35),]

###Note that these chunks should be a function because we're calling the same thing more than twice, I am just feeling lazy and don't want to write the function lol 

hist(FinalDataset$BG_assay)
boxplot(FinalDataset$BG_assay ~ FinalDataset$Crop)
which(FinalDataset$BG_assay > 800)
#303 and 304 are outliers, both from batch 1 
boxplot(BG_assay ~ Batch, data = FinalDataset)
BG <- ggplot(FinalDataset, aes(x = Crop, y=BG_assay))+
  geom_boxplot(fill = "skyblue")+ 
  labs(x=NULL, y= NULL)+
  ggtitle("β-glucosidase Activity")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14))
BG 


hist(FinalDataset$Cello_assay)
boxplot(Cello_assay ~ Crop, data = FinalDataset)
which(FinalDataset$Cello_assay > 100)
#207, 304, 315, 401, 511 are outliers in their sections but likely not in the whole sample
boxplot(Cello_assay ~ Batch, data = FinalDataset)
Cello <- ggplot(FinalDataset, aes(x = Crop, y=Cello_assay))+
  geom_boxplot(fill = "skyblue")+ 
  labs(x=NULL, y= NULL)+
  ggtitle("Cellobiohydrolase Activity")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14))
Cello 

hist(FinalDataset$NAG_assay)
boxplot(NAG_assay ~ Crop, data = FinalDataset)
which(FinalDataset$NAG_assay >130)
#103, 211, 213 are outliers 
boxplot(NAG_assay ~ Batch, data = FinalDataset)
NAG <- ggplot(FinalDataset, aes(x = Crop, y=NAG_assay))+
  geom_boxplot(fill = "skyblue")+ 
  labs(x=NULL, y= NULL)+
  ggtitle("N-acetylglucosaminidase Activity")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14))
NAG

hist(FinalDataset$P_assay)
boxplot(P_assay ~ Crop, data = FinalDataset)
which(FinalDataset$P_assay >2000)
#306 and 401 
which(FinalDataset$P_assay < 500)
#307 is low 
boxplot(P_assay ~ Batch, data = FinalDataset)
P <- ggplot(FinalDataset, aes(x = Crop, y=P_assay))+
  geom_boxplot(fill = "skyblue")+ 
  labs(x=NULL, y= NULL)+
  ggtitle("Phosphatase Activity")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14))
P

#304 and 401 are repeat offenders...

grid.arrange(BG, Cello, NAG, P, nrow = 2, bottom = "Cropping System - 2021 Data", left = text_grob(expression(Enzyme~Activity~nmol~hr^-1~g~soil^-1), rot = 90))


## AOV models 
### Note that good stats says that if you don't see any potential in the visuals above, you don't need to run any stats tests. In fact, you shouldn't, because you will get false positives. 
PAov <- aov(P_assay ~ Crop, data = FinalDataset)
summary(PAov)
PAov_Block <- aov(P_assay ~ Crop + Block, data = FinalDataset)
summary(PAov_Block)

###Compare models 
anova(PAov, PAov_Block)
###Prefer the reduced model 

PAov_Batch <- aov(P_assay ~ Crop + Batch, data = FinalDataset)
anova(PAov, PAov_Batch)
#If we follow p values, we prefer reduced model, but it really seems like Batch adds something :/

PAov_BlockInt <- aov(P_assay ~ Crop + Block + Batch, data = FinalDataset)
anova(PAov_Batch, PAov_BlockInt)
#Prefer reduced model 


TukeyHSD(PAov_Batch)


##Save as csv 

write.csv(FinalDataset, "/Users/emmalink/Documents/R/PARCE/Analysis 2020/EEA_Analysis /EEA_IntermediateDataFiles/FinalDataset.csv")
