library("ggplot2")
library("dplyr")
library("IRdisplay")
library("readxl")
library("plyr")


#=========For Negative=============


#set the R working directory
#show files in directory
#Reading R Data Files+ naming the files
#metadata
MD <- read.csv("Metadata.csv")
head(MD)


#Quant_Table output 
FT <- read.csv("Quantification_Table.csv")
head(FT)

colnames(FT) <- sub("X","",colnames(FT))
#Selecting the columns with samples
New_FT <- select(FT, matches("mzML"))

#Removing "peak area" from column names (to match MD)
colnames(New_FT) <- gsub(".Peak.area","",colnames(New_FT))

#Changing names in metadata (to match New_FT)
MD$filename <- gsub('-','.',MD$filename)

#Finding blanks
Blank <- filter(MD, SampleType == "Blank")
Blank_data <- select(New_FT, c(Blank$filename))

#Finding samples
Sample <- filter(MD, SampleType == "SAMPLE")
Sample_data <- select(New_FT, c(Sample$filename))

#Adding average to blanks and samples and adding row.ID
Blank_data$MeanB <-apply(Blank_data,1,mean)
Sample_data$MeanS <-apply(Sample_data,1,mean)
Blank_data <- data.frame(row.ID = select(FT, row.ID), mz = select(FT, row.m.z), RT = select(FT, row.retention.time), Blank_data)
Sample_data <- data.frame(row.ID = select(FT, row.ID), mz = select(FT, row.m.z), RT = select(FT, row.retention.time), Sample_data)

#Combining averages into new file
BSaverage <- data.frame(row.ID = select(FT, row.ID), mz = select(FT, row.m.z),RT = select(FT, row.retention.time), MeanB = select(Blank_data, MeanB), MeanS = select(Sample_data, MeanS))

#Blank/Sample ratio
BSaverage$ratio <- BSaverage$MeanB/BSaverage$MeanS

#if ratio is bigger than 2 then it's 0, if not then it's 1
BSaverage <- BSaverage %>% 
  mutate(bin = if_else(BSaverage$ratio > 2, 0, 1))

#Number of features after blank removal
FeaturesRemoved <- (nrow(BSaverage)-sum(BSaverage$bin))
FeaturesRemained <- sum(BSaverage$bin)

#Filtering all samples and blank removal
Sample_data <- subset (Sample_data, select = -MeanS)
rownames(Sample_data) <- Sample_data[,1]
Sample2 <- filter(BSaverage, bin == "1")
Sample2 <- data.frame(Sample2[,-1], row.names=Sample2[,1])
Sample_clean <- Sample_data %>% filter(rownames(Sample_data) %in% rownames(Sample2))
colnames(Sample_clean) <- sub("X","",colnames(Sample_clean))
write.csv(Sample_clean,"Data__after_Cleanup.csv", row.names = FALSE)

#===============

RNA_MD <- read.csv("IBD_samples_for_16rna.csv")
RNA_data <- read.csv("16RNA_Abundance_Results.csv")
colnames(RNA_data) <- sub("X","",colnames(RNA_data))
#Selecting the columns with samples
MS_Data <- select(Sample_clean, matches("mzML"))

# Create proper mapping
name_mapping <- setNames(RNA_MD$filename, as.character(RNA_MD$RNA_sampleID))

# Get current column names and create new ones
old_names <- colnames(RNA_data)
new_names <- old_names

# Replace only numeric columns with mapped filenames
numeric_cols <- suppressWarnings(!is.na(as.numeric(old_names)))
new_names[numeric_cols] <- name_mapping[old_names[numeric_cols]]

# Create new dataframe with renamed columns
NEW_RNA_data <- RNA_data
colnames(NEW_RNA_data) <- new_names
#Selecting the columns with samples
NEW_RNA_data <- select(NEW_RNA_data, matches("mzML"))
NEW_RNA_data <- data.frame(Bacteria = select(RNA_data, Name), NEW_RNA_data)
colnames(NEW_RNA_data) <- sub("X","",colnames(NEW_RNA_data))
# Write updated data
write.csv(NEW_RNA_data, "16RNA_results_renamed.csv", row.names = FALSE)

