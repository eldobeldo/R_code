#Source file MS Analysis; output files include Distance Matrix, PCA, and allele counts

rm(list=ls())           #Clear workspace == ensure no leftover variables
graphics.off();         #Close all graphics windows

#Need to install appropriate packages; combinat and polysat 
#Run this code: install.packages(combinat), then install.packages(polysat)

require(polysat)        #Loads Polysat package

getwd()                 #Get current working directory

setwd("/Users/eldin/Desktop/R working dir/")

search()                #See current packages loaded

#Load your data into polysat; it must be in the correct format! See README and SampleFormat.txt file

MyMSdata <- read.GeneMapper("Input_data.txt")

summary(MyMSdata)       #Get summary of your data

Samples(MyMSdata)       #List your samples 

Loci(MyMSdata)          #List your Loci information

find.missing.gen(MyMSdata) #Find which had missing MS data

#Add Description to the dataset

Description(MyMSdata) <-"Thailand K13 Data" 

PopNames(MyMSdata) <- c("PopA")

#PopInfo(MyMSdata) <- rep(1:2, c(32, 365))  #Add pop info, 1 = EAST, 2 = WEST; code == repeat 1, 32 times, followed by repeat 2, 65 times

Usatnts(MyMSdata) <- c(2)                 #Set nucleotide repeats for data

#Data Anaysis; Bruvo takes mutations into account, while Lynch does not! By default == Bruvo

GeneticDistanceMatrix_Bruvo <- meandistance.matrix(MyMSdata) #Calculate genetic distances matrix between samples

#To test whether you have NA values anywhere, run is.na(YourData); It will return TRUE/FALSE values
#Sometimes the NA values don't show up, but are just empty cells. To fill them, use x[(x=="")] <- NA, to fill with NA
#Now you can replace the NA with 0s, use x[is.nan(x)] = 0

is.na(GeneticDistanceMatrix_Bruvo)

GeneticDistanceMatrix_Bruvo[(GeneticDistanceMatrix_Bruvo=="")] <- NA

GeneticDistanceMatrix_Bruvo[is.nan(GeneticDistanceMatrix_Bruvo)] = 0 

PCA <- cmdscale(GeneticDistanceMatrix_Bruvo)                 #Perform Principal Coordinate Analysis (PCA)


GeneticDistanceMatrix_Lynch <- meandistance.matrix(MyMSdata)

#Run the same code for Lynch method; fix empty cells issue

is.na(GeneticDistanceMatrix_Lynch)

GeneticDistanceMatrix_Lynch[(GeneticDistanceMatrix_Lynch=="")] <- NA

GeneticDistanceMatrix_Lynch[is.nan(GeneticDistanceMatrix_Lynch)] = 0 


PCA2 <- cmdscale(GeneticDistanceMatrix_Lynch)   

simal <- alleleDiversity(MyMSdata)  #Perform allele diversity and frequency analysis fore each loci

AlleleCounts <- simal$counts        #Save allele counts into allelecount variable

#Export Genetic Distance Matrix and PCA

#dir.create("Results")             #Create results directory in current wd()

write.table(GeneticDistanceMatrix_Bruvo, file="GeneticDistanceMatrix_BruvoMethod.txt", quote=FALSE, sep="\t", na="", row.names=TRUE, col.names=FALSE)

write.table(PCA, file="PCAresults_BruvoMethod.txt", quote=FALSE, sep="\t", na="", row.names=TRUE, col.names=FALSE)

write.table(GeneticDistanceMatrix_Lynch, file="GeneticDistanceMatrix_LynchMethod.txt", quote=FALSE, sep="\t", na="", row.names=TRUE, col.names=TRUE)

write.table(PCA2, file="PCAresults_LynchMethod.txt", quote=FALSE, sep="\t", na="", row.names=TRUE, col.names=FALSE)

write.table(AlleleCounts, file="AlleleCounts.txt", quote=FALSE, sep="\t", na="", row.names=TRUE, col.names=FALSE)



#Save current workspace 

save(MyMSdata, file="MSanalysis.RData")



