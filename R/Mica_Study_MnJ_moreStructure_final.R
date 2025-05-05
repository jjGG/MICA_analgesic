# R

# Functional Genomics Center 2024
#
# stats are also important in the real world ;)
#
# this script is developed for Michael Straessle together with JG
# set working directory to source file location
#
# @ jg@fgcz.ethz.ch
#

# libraries
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(gplots)
library(matrixStats)

# source functions
source("UsefulFunctions_MicaStudy.R")


#read in annotation file
annoFile_here <- read_tsv("../input/headerNAnnotations_2024_11_06_validatedFortxtComp_wExclud.txt")

# data
dataFileName <- "../input/data_close2final_MICA_2024-08_complete.txt"
mydat_here  <- read_tsv(dataFileName)

# work on rows
# data with excluded patients (non consent)
boolExclude <- mydat_here$ExcludeFromStudy == "x"
includeInStudy_bool <- is.na(boolExclude) # make sure that the cells are empty therefore here seen as NA
table(includeInStudy_bool)

# > table(includeInStudy_bool)
# includeInStudy_bool
# FALSE  TRUE
# 93    98

mydat_here <- mydat_here[includeInStudy_bool, ]

# work on columns
# only keep these columns where anno file == exclude from study = "x"
keepColumn_bool <- annoFile_here$useIt == "x"
annoToUse <- annoFile_here[which(annoFile_here$useIt == "x"),]
keepColumn_bool[is.na(keepColumn_bool)] <- FALSE
sum(keepColumn_bool, na.rm = T)
keep_idx <- which(annoFile_here$useIt == "x")

# data file to only relevant columns
relevantDat <- mydat_here[,keep_idx]

# check complications manually -> write out to tsv
length(grep(x = colnames(relevantDat), pattern = "txtcom"))
gg <- relevantDat[,c(1,grep(x = colnames(relevantDat), pattern = "txtcomp"))]
write_tsv(x = gg, file = "../output/output_pMica_complicationsTable_2024-11-1.tsv")

# write out some special columns
columnsOI <- c("rboschmerzzeit", "chkdafalgan", "cbodosisdafalgan", "chknovalgin", "cbodosisnovalgin", "chktramal", "cbodosistramal", "chkandere", "txtandere" ,"cbodosisandere", "radschmerzmvorop", "txtwelcheschmerzm")
resultColsOI <- matrix(nrow = nrow(relevantDat))
for (i in 1:length(columnsOI)) {
  idx <- grep(pattern = columnsOI[i], x = colnames(relevantDat))
  myColsResOfI <- relevantDat[,idx]
  resultColsOI <- cbind(resultColsOI, myColsResOfI)
}

# fix for first col
resultColsOI <- resultColsOI[, -1 ]
dim(resultColsOI)
write_tsv(resultColsOI, file = "../output/output_Analgetic_intake_20241030.tsv")
write_tsv(relevantDat[,c(grep(x = colnames(relevantDat), pattern = "txtcom") )], file = "../output/output_ShortWay_complications_20241030.tsv")



# work on dmaa
columnsOI2 <- c("dmaa")
resultColsOI2 <- matrix(nrow = nrow(relevantDat))
for (i in 1:length(columnsOI2)) {
  idx <- grep(pattern = columnsOI2[i], x = colnames(relevantDat))
  myColsResOfI <- relevantDat[,idx]
  resultColsOI2 <- cbind(resultColsOI2, myColsResOfI)
}

# fix for first col
resultColsOI2 <- resultColsOI2[, -1 ]

pre_idx <- grep(x = colnames(resultColsOI2), pattern = "pre")
colnames(resultColsOI2)[pre_idx]

lastCOI2 <- keepLastValueFromTimeResolvedResultDF(resultColsOI2[,pre_idx])

# get pre means
pre_mean <- mean(lastCOI2, na.rm = TRUE)
pre_SD <- sd(lastCOI2, na.rm = TRUE)
pre_median <- median(lastCOI2, na.rm = TRUE)
pre_min <- min(lastCOI2, na.rm = TRUE)
pre_max <- max(lastCOI2, na.rm = TRUE)

post_idx <- grep(x = colnames(resultColsOI2), pattern = "post")
colnames(resultColsOI2)[post_idx]

lastCOI3 <- keepLastValueFromTimeResolvedResultDF(resultColsOI2[,post_idx])
# get post means
post_mean <- mean(lastCOI3, na.rm = TRUE)
post_SD <- sd(lastCOI3, na.rm = TRUE)
post_median <- median(lastCOI3, na.rm = TRUE)
post_min <- min(lastCOI3, na.rm = TRUE)
post_max <- max(lastCOI3, na.rm = TRUE)

# is it different?
boxplot(lastCOI2, lastCOI3)

moreRes <- cbind(pre_mean, pre_median, pre_min, pre_max, post_mean, post_median, post_min, post_max)

# write out results
write.table(moreRes, file = "../output/output_PrePost_DMAA_241030.tsv", sep="\t", quote = FALSE)


hist(lastCOI3, main = "DMAA", xlab = "DMAA post", ylab = "Frequency")
hist_data <- hist(lastCOI3, plot = FALSE)


# plotting
pdf("../output/output_DMAA-post.pdf")
# Plot histogram
hist(lastCOI3, main="",xlab = "DMAA post", ylab = "Frequency")
# Add counts on top of the bars
for (i in 1:length(hist_data$counts)) {
  text(x = hist_data$mids[i],
       y = hist_data$counts[i],
       labels = hist_data$counts[i],
       pos = 3)
}
dev.off()


pdf("../output/output_DMAA-pre.pdf")
# Plot histogram
hist(lastCOI2, main="",xlab = "DMAA pre", ylab = "Frequency")

hist_data2 <- hist(lastCOI2, plot = FALSE)

# Add counts on top of the bars
for (i in 1:length(hist_data2$counts)) {
  text(x = hist_data2$mids[i],
       y = hist_data2$counts[i],
       labels = hist_data2$counts[i],
       pos = 3)
}
dev.off()


########
#
#  EQ5D5L
#
########

# which ones are to be further processed
keepRowsToDecode_bool <- annoFile_here$FurtherProcess == "xx"
keepRowsToDecode_bool[is.na(keepRowsToDecode_bool)] <- FALSE
table(keepRowsToDecode_bool)
# keepRowsToDecode_bool
# FALSE  TRUE
# 1534   196
annoToDecodeQuestions <- annoFile_here[keepRowsToDecode_bool, ]

#sapply(strsplit(myMat$`Pep-Phos-Prot`, split = "~"), function(x)x[4])
annoToDecodeQuestions$TimeCode <- sapply(strsplit(annoToDecodeQuestions$mainHeader, split = "_"), function(x)x[2])
unique(paste(annoToDecodeQuestions$TimeCode, annoToDecodeQuestions$weeks,annoToDecodeQuestions$question))


# decoding from Matrix with annotations
# here we need to look into!
myDecodedQs <- DecodeQuestionForOnlyAFewQWithAnnoFile(annoFile = annoToDecodeQuestions, mydat = relevantDat, keyWord = "yes")




#myDecodedQs <- DecodeQuestionMatrixWithannoFileFile(annoFile = annoToDecodeQuestions, mydat = relevantDat)
myDecodedQsNoSum <- DecodeQuestionForOnlyAFewQWithAnnoFileNoSummary(annoFile = annoToDecodeQuestions, mydat = relevantDat, keyWord = "yes")
myDecodedQsNoSum <- data.frame(myDecodedQsNoSum)

# split for shoes and EQ5DL
myShoeColumns <- myDecodedQsNoSum[,grep(x = colnames(myDecodedQsNoSum), pattern = "Schuhwerk")]

# process here with coding table to summarize all!
myEQ5D5Lcolumns <- myDecodedQsNoSum[, -grep(x = colnames(myDecodedQsNoSum), pattern = "Schuhwerk")]
dim(myEQ5D5Lcolumns)

# read in coding table
myCodingTable <- read_tsv("../input/EQ5DLcodingTable.txt")
str(myCodingTable)
relevantCoding <- cbind(myCodingTable$`5L profile`, myCodingTable$Germany)
colnames(relevantCoding) <- c("myStringCode", "myNumber")
relevantCoding <- data.frame(relevantCoding)


# how many timePoints
numTimePoint <- 7

myEQ5D5Lcodes <- matrix(nrow = nrow(myEQ5D5Lcolumns), ncol = numTimePoint)
# look up and fill matrix
for (i in 1:numTimePoint) {
  j <- i*5 - 4
  print(paste(j, "--", j+4))
  mySubEQ5D5L <- myEQ5D5Lcolumns[,j:(j+4)]
  #paste(c("hi", x), collapse = " ")
  for (jj in 1:nrow(mySubEQ5D5L)) {
    myStringWord <- paste(mySubEQ5D5L[jj,], collapse = "")
    myLookUp <- try(relevantCoding$myNumber[which(relevantCoding$myStringCode == myStringWord)])
    if (length(myLookUp) > 0) {
      myEQ5D5Lcodes[jj,i] <- myLookUp
    } else {
      myEQ5D5Lcodes[jj,i] <- NA
    }
  }
}

# label cols
colnames(myEQ5D5Lcodes) <- paste("EQ5D5Lvalues", "timepoint", 1:numTimePoint)
myEQ5D5LLastEntryInTime <- keepLastValueFromTimeResolvedResultDF(myEQ5D5Lcodes)

# we do not wanna loose the preop value
# this is only the first value if available
myEQ5D5LPreOpValue <- myEQ5D5Lcodes[,1]

########
#
#  FFI
#
########


# work on FFI
str(annoFile_here)
FFI_idx <- which(annoFile_here$group == "FFI" & annoFile_here$useIt == "x")
FFI_anno <- annoFile_here[FFI_idx, ]
FFI_anno$TimeCode <- sapply(strsplit(FFI_anno$mainHeader, split = "_"), function(x)x[2])
FFI_anno$groupingCode <- paste(FFI_anno$weeks, FFI_anno$group, FFI_anno$TimeCode)
unique(FFI_anno$groupingCode)

myFFIOverTime <- matrix(nrow = nrow(myEQ5D5Lcolumns), ncol = numTimePoint)

# look up and fill matrix
for (i in 1:length(unique(FFI_anno$groupingCode))) {
  subGrp <- unique(FFI_anno$groupingCode)[i]
  myHeaders <- FFI_anno$mainHeader[FFI_anno$groupingCode == subGrp]
  myFFIblock <- relevantDat[,myHeaders]
  # insert check if all are filled and otherwise ommit row completely
  #if (sum(is.na(as.numeric(myFFIblock[i,]))) == 0) myFFIOverTime[,i] <- rowSums(myFFIblock[i,], na.rm = TRUE)
  myFFIOverTime[,i] <- rowSums(myFFIblock, na.rm = TRUE) # keep track of NAs below
  # overwrite row sums with "NA" because not all are filled and therefore rowSum is not adquate
  boolAllComplete <- (rowSums(is.na(myFFIblock)) > 0) # new to adjust for sporadic filler
  myFFIOverTime[boolAllComplete,i] <- NA # new to adjust for sporadic filler
}

colnames(myFFIOverTime) <- paste("FFIsums", "timepoint", 1:numTimePoint)
# normalize for the number of Qs /230 * 100 = /2.3
myFFIOverTime_norm <- myFFIOverTime/2.3
myFFILastEntryInTime <- keepLastValueFromTimeResolvedResultDF(myFFIOverTime_norm)
myFFIPreOpValue <- myFFIOverTime_norm[,1]

# boxplot customized
# summarize the last two columns together and plot it together with the other columns without plotting the last two but the sum of these two
myFFIOverTime_norm <- cbind(myFFIOverTime_norm, rowMeans(myFFIOverTime_norm[,6:7], na.rm = TRUE))
colnames(myFFIOverTime_norm)[8] <- "FFIsums timepoint 6+7"
boxplot(myFFIOverTime_norm[,c(1:5,8)], main = "FFI", xlab = "timepoint", ylab = "FFI", col = "lightblue")

# install.packages("vioplot")
library(vioplot)
vioplot(myFFIOverTime_norm[,c(1:5,8)], main = "FFI", xlab = "timepoint", ylab = "FFI", col = "lightblue")


# Differenzen Prae and Post OP
myFFI_diff <- myFFILastEntryInTime - myFFIPreOpValue


# split FFI into subclasses
# look up and fill matrix
myPatternVector <- c("last 4 weeks", "complains during", "of foot complain")
descriptVector <- c("pain", "difficulty", "disability")
mySubclasses <- paste("FFIsubScales", descriptVector)
# write proper colnames in the right order
mySubClassGrps <- rep(c("last 4 weeks", "complains during", "of foot complain"), 3)
myDivisor <- c(0.9, 0.9, 0.5)

myTimesInThree <- sort(rep(c(1:numTimePoint), 3))
myColnameSubclassesFFI <- paste(rep(mySubclasses,7),myTimesInThree)


# join in categories from annoToUse to match column names
myBlocknames <- as_tibble(colnames(myFFIblock))
colnames(myBlocknames) <- "mainHeader"
myLittleFrame <- left_join(myBlocknames, annoToUse)
myMatch <- myPatternVector[1]
getIDx <- grep(x = myLittleFrame$category, pattern = myMatch)


#
myFFIResultFrame <- matrix(nrow = nrow(myFFIblock), ncol = length(myPatternVector)*numTimePoint)
colnames(myFFIResultFrame) <- myColnameSubclassesFFI
myColIterator <- 0
for (i in 1:length(unique(FFI_anno$groupingCode))) {
  myColIterator <- myColIterator + 1
  subGrp <- unique(FFI_anno$groupingCode)[i]
  myHeaders <- FFI_anno$mainHeader[FFI_anno$groupingCode == subGrp]
  myFFIblock <- relevantDat[,myHeaders]

  for (j in 1:length(myPatternVector)) {
    myBlocknames <- as_tibble(colnames(myFFIblock))
    colnames(myBlocknames) <- "mainHeader"
    myLittleFrame <- left_join(myBlocknames, annoToUse)
    myMatch <- myPatternVector[j]
    getIDx <- grep(x = myLittleFrame$category, pattern = myMatch)
    # insert check if all are filled and otherwise ommit row completely
    myBlocknames <- as_tibble(colnames(myFFIblock))
    colnames(myBlocknames) <- "mainHeader"
    myLittleFrame <- left_join(myBlocknames, annoToUse)
    if (j == 1) {
      myColIterator <- myColIterator
      myFFIblockSub <- myFFIblock[,getIDx]
      myFFIsubBlockSums <- rowSums(myFFIblockSub, na.rm = FALSE)/myDivisor[j]
      # overwrite row sums with "NA" because not all are filled and therefore rowSum is not adquate
      myFFIResultFrame[,myColIterator] <- myFFIsubBlockSums

    }
    if (j == 2)  {
      myColIterator <- myColIterator + 1
      myFFIblockSub <- myFFIblock[,getIDx]
      myFFIsubBlockSums <- rowSums(myFFIblockSub, na.rm = FALSE)/myDivisor[j]
      # overwrite row sums with "NA" because not all are filled and therefore rowSum is not adquate
      myFFIResultFrame[,myColIterator] <- myFFIsubBlockSums
    }
    if (j == 3) {
      myColIterator <- myColIterator + 1
      myFFIblockSub <- myFFIblock[,getIDx]
      myFFIsubBlockSums <- rowSums(myFFIblockSub, na.rm = FALSE)/myDivisor[j]
      # overwrite row sums with "NA" because not all are filled and therefore rowSum is not adquate
      myFFIResultFrame[,myColIterator] <- myFFIsubBlockSums
    }
  }

}

# boxplotting
# FFI subscale pain
colnames(myFFIResultFrame[,c(1,4,7,10,13,16,19)])
boxplot(myFFIResultFrame[,c(1,4,7,10,13,16,19)], main = "FFIsubscale pain", xlab = "timepoint", ylab = "FFIsub pain", col = "lightblue")

myFFIOverTime_pain_norm <- cbind(myFFIResultFrame[,c(1,4,7,10,13)], rowMeans(myFFIResultFrame[,c(16,19)], na.rm = TRUE))
colnames(myFFIOverTime_pain_norm)[6] <- "FFIsums_pain timepoint 6+7"
boxplot(myFFIOverTime_pain_norm, main = "FFIsubscale_pain", xlab = "timepoint", ylab = "FFI", col = "lightblue")

# FFI subscale difficulty
colnames(myFFIResultFrame[,c(2,5,8,11,14,17,20)])
boxplot(myFFIResultFrame[,c(2,5,8,11,14,17,20)], main = "FFIsubscale difficulty", xlab = "timepoint", ylab = "FFIsub difficulty", col = "lightblue")

myFFIOverTime_difficulty_norm <- cbind(myFFIResultFrame[,c(2,5,8,11,14)], rowMeans(myFFIResultFrame[,c(17,20)], na.rm = TRUE))
colnames(myFFIOverTime_difficulty_norm)[6] <- "FFIsums_difficulty timepoint 6+7"
boxplot(myFFIOverTime_difficulty_norm, main = "FFIsubscale_difficulty", xlab = "timepoint", ylab = "FFI", col = "lightblue")

# FFI subscale disability
colnames(myFFIResultFrame[,c(3,6,9,12,15,18,21)])
boxplot(myFFIResultFrame[,c(3,6,9,12,15,18,21)], main = "FFIsubscale disability", xlab = "timepoint", ylab = "FFIsub disability", col = "lightblue")

myFFIOverTime_disability_norm <- cbind(myFFIResultFrame[,c(3,6,9,12,15)], rowMeans(myFFIResultFrame[,c(18,21)], na.rm = TRUE))
colnames(myFFIOverTime_disability_norm)[6] <- "FFIsums_disability timepoint 6+7"
boxplot(myFFIOverTime_disability_norm, main = "FFIsubscale_disability", xlab = "timepoint", ylab = "FFI", col = "lightblue")



########
#
#  MOXFQ
#
########


# work on MOXFQ -> a bit more cumbersome since "words"..
MOXFQ_idx <- which(annoFile_here$group == "MOXFQ" & annoFile_here$useIt == "x")
MOXFQ_anno <- annoFile_here[MOXFQ_idx, ]
MOXFQ_anno$TimeCode <- sapply(strsplit(MOXFQ_anno$mainHeader, split = "_"), function(x)x[2])
MOXFQ_anno$groupingCode <- paste(MOXFQ_anno$weeks, MOXFQ_anno$group, MOXFQ_anno$TimeCode)
unique(MOXFQ_anno$groupingCode)


myMOFXOverTime <- matrix(nrow = nrow(myEQ5D5Lcolumns), ncol = numTimePoint)

# look up and fill matrix
for (i in 1:length(unique(MOXFQ_anno$groupingCode))) {
  subGrp <- unique(MOXFQ_anno$groupingCode)[i]
  myHeaders <- MOXFQ_anno$mainHeader[MOXFQ_anno$groupingCode == subGrp]
  myMOXFQblock <- relevantDat[,myHeaders]
  # recode column by column (the last column has different words)
  #myMOXFQOverTime[,i] <- rowSums(myMOXFQblock)
  myResultBlock <- matrix(nrow = nrow(myMOXFQblock), ncol = ncol(myMOXFQblock))
  for (jj in 1:ncol(myMOXFQblock)) {
    myV <- myMOXFQblock %>% pull(jj)
    #print(table(myV))
    mre <- dplyr::recode(myV, "immer" =  4,"meistens" = 3, "manchmal" = 2, "selten" = 1, "nie" = 0)
    if (length(unique(mre)) == 1) mre <-  dplyr::recode(myV, "schwer" =  4,"m\x8assig" = 3, "mild" = 2, "sehr mild" = 1, "keine" = 0)
    #print(table(mre))
    myResultBlock[,jj] <- mre
  }
  myMOFXOverTime[,i] <- rowSums(myResultBlock, na.rm = TRUE) # keep track of NAs below
  # check if there are some rows where there are few NAs (sporadic non-filler)
  boolAllComplete <- (rowSums(is.na(myResultBlock)) > 0) # new to adjust for sporadic non-filler
  myMOFXOverTime[boolAllComplete,i] <- NA # new to adjust for sporadic non-filler
}

#bak <- myMOFXOverTime
myMOFXsumDF <- data.frame(myMOFXOverTime)
colnames(myMOFXsumDF) <- paste("MOFXQsum",unique(MOXFQ_anno$groupingCode))

# boxplot over time MOXFQ
myMOFXsumDF <- cbind(myMOFXsumDF, rowMeans(myMOFXsumDF[,6:7], na.rm = TRUE))
colnames(myMOFXsumDF)[8] <- "MOXFQsums timepoint 6+7"
boxplot(myMOFXsumDF[,c(1:5,8)], main = "MOXFQ", xlab = "timepoint", ylab = "MOXFQ", col = "lightblue")

# normalize
# sum /64 * 100
myMOFXsumDF_norm <- myMOFXsumDF /64 * 100
myMOFXLastEntryInTime <- keepLastValueFromTimeResolvedResultDF(myMOFXsumDF_norm)
myMOFXPreOpValue <- myMOFXsumDF_norm[,1]
boxplot(myMOFXsumDF_norm[,c(1:5,8)], main = "MOXFQ", xlab = "timepoint", ylab = "MOXFQ", col = "lightblue")


# Differenzen Prae and Post OP
myMOFX_diff <- myMOFXLastEntryInTime - myMOFXPreOpValue

###########################################
# split MOXFQ into subclasses
# look up and fill matrix
myPatternVector <- c("painMOXFQ", "walking_n_standing", "social_interaction")
descriptVector <- c("pain", "difficulty", "disability")
mySubclasses <- paste("MOXFQsubScales", descriptVector)
# write proper colnames in the right order
mySubClassGrps <- rep(c("painMOXFQ", "walking_n_standing", "social_interaction"), 3)
myDivisor <- c(0.2, 0.28, 0.16)

myTimesInThree <- sort(rep(c(1:numTimePoint), 3))
myColnameSubclassesMOXFQ <- paste(rep(mySubclasses,7),myTimesInThree)


# join in categories from annoToUse to match column names
myBlocknames <- as_tibble(colnames(myMOXFQblock))
colnames(myBlocknames) <- "mainHeader"
myLittleFrame <- left_join(myBlocknames, annoToUse)
myMatch <- myPatternVector[1]
getIDx <- grep(x = myLittleFrame$category, pattern = myMatch)


#
myMOXFQResultFrame <- matrix(nrow = nrow(myMOXFQblock), ncol = length(myPatternVector)*numTimePoint)
colnames(myMOXFQResultFrame) <- myColnameSubclassesMOXFQ

subScaleVector <- list(c(1,11,12,15,16), c(2:8), c(9,10,13,14))
mySubScaleResultFrame <- matrix(nrow = nrow(myMOXFQblock), ncol = length(subScaleVector)*numTimePoint)

myColIterator <- 0
rm(myMOXFQblock) # remove the old block
myColIterator <- 1

for (i in 1:length(unique(MOXFQ_anno$groupingCode))) { # loop over timepoints (7)
  subGrp <- unique(MOXFQ_anno$groupingCode)[i]
  myHeaders <- MOXFQ_anno$mainHeader[MOXFQ_anno$groupingCode == subGrp]
  myMOXFQblock <- relevantDat[,myHeaders]
  rm(myResultBlock2)
  myResultBlock2 <- matrix(nrow = nrow(myMOXFQblock), ncol = ncol(myMOXFQblock))
  #mySubscaleResult <- matrix(nrow = nrow(myMOXFQblock), ncol = length(subScaleVector))
  for (jj in 1:ncol(myMOXFQblock)) {
    myV <- myMOXFQblock %>% pull(jj)
    #print(table(myV))
    mre <- dplyr::recode(myV, "immer" =  4,"meistens" = 3, "manchmal" = 2, "selten" = 1, "nie" = 0)
    if (length(unique(mre)) == 1) mre <-  dplyr::recode(myV, "schwer" =  4,"m\x8assig" = 3, "mild" = 2, "sehr mild" = 1, "keine" = 0)
    #print(table(mre))
    myResultBlock2[,jj] <- mre
  }
  #colnames(myMOXFQblock)
  #myMOXFQOverTime[,i] <- rowSums(myMOXFQblock)
  # go into the block and extract the columns (subscale)
  # subscale 1 = 1,11,12,15,16 # pain
  # subscale 2 = 2:8 # walking_n_standing
  # subscale 3 = 9,10,13,14 # social_interaction
  #
  for (jjj in 1:length(subScaleVector)) {
    mySubScaleResultFrame[,myColIterator] <- rowSums(myResultBlock2[,subScaleVector[[jjj]]])/myDivisor[jjj]
    myColIterator <- myColIterator + 1
    #mySubscaleBlock <- myResultBlock2[,subScaleVector[[jjj]]]
    #mySubscaleResult[,jjj] <- rowSums(mySubscaleBlock, na.rm = FALSE)/myDivisor[jjj]
  }

}

# finally resultDF w/ 7 * 3 subscales = 21 columns (1,2,3) for t1 and 4,5,6 for t2 and so on
dim(mySubScaleResultFrame)

# boxplot over time MOXFQ subscales
# do boxplotting for individual subscales per timepoint

# MOXFQ subscale pain
painMOXFQ <- mySubScaleResultFrame[,c(1,4,7,10,13,16,19)]
walkingNstanding <- mySubScaleResultFrame[,c(2,5,8,11,14,17,20)]
socialInteraction <- mySubScaleResultFrame[,c(3,6,9,12,15,18,21)]

# boxplot painMOXFQ
boxplot(painMOXFQ, main = "MOXFQsubscale pain", xlab = "timepoint", ylab = "MOXFQsub pain", col = "lightblue")

painMOXFQ <- cbind(mySubScaleResultFrame[,c(1,4,7,10,13)], rowMeans(mySubScaleResultFrame[,c(16,19)], na.rm = TRUE))
colnames(painMOXFQ)[6] <- "painMOXFQ timepoint 6+7"
colnames(painMOXFQ)[1] <- "painMOXFQ timepoint 1"
colnames(painMOXFQ)[2] <- "painMOXFQ timepoint 2"
colnames(painMOXFQ)[3] <- "painMOXFQ timepoint 3"
colnames(painMOXFQ)[4] <- "painMOXFQ timepoint 4"
colnames(painMOXFQ)[5] <- "painMOXFQ timepoint "
boxplot(painMOXFQ, main = "MOXFQsubscale pain", xlab = "timepoint", ylab = "MOXFQsub pain", col = "lightblue")

# boxplot walkingNstanding
boxplot(walkingNstanding, main = "MOXFQsubscale walking_n_standing", xlab = "timepoint", ylab = "MOXFQsub walking_n_standing", col = "lightblue")

walkingNstanding <- mySubScaleResultFrame[,c(2,5,8,11,14,17,20)]

walkingNstanding <- cbind(mySubScaleResultFrame[,c(2,5,8,11,14)], rowMeans(mySubScaleResultFrame[,c(17,20)], na.rm = TRUE))
colnames(walkingNstanding)[6] <- "walkingNstanding timepoint 6+7"
colnames(walkingNstanding)[1] <- "walkingNstanding timepoint 1"
colnames(walkingNstanding)[2] <- "walkingNstanding timepoint 2"
colnames(walkingNstanding)[3] <- "walkingNstanding timepoint 3"
colnames(walkingNstanding)[4] <- "walkingNstanding timepoint 4"
colnames(walkingNstanding)[5] <- "walkingNstanding timepoint 5"
boxplot(walkingNstanding, main = "MOXFQsubscale walkingNstanding", xlab = "timepoint", ylab = "MOXFQsub walkingNstanding", col = "lightblue")

# boxplot social_interaction
boxplot(socialInteraction, main = "MOXFQsubscale social_interaction", xlab = "timepoint", ylab = "MOXFQsub social_interaction", col = "lightblue")

socialInteraction <- mySubScaleResultFrame[,c(3,6,9,12,15,18,21)]

socialInteraction <- cbind(mySubScaleResultFrame[,c(3,6,9,12,15)], rowMeans(mySubScaleResultFrame[,c(18,21)], na.rm = TRUE))
colnames(socialInteraction)[6] <- "socialInteraction timepoint 6+7"
colnames(socialInteraction)[1] <- "socialInteraction timepoint 1"
colnames(socialInteraction)[2] <- "socialInteraction timepoint 2"
colnames(socialInteraction)[3] <- "socialInteraction timepoint 3"
colnames(socialInteraction)[4] <- "socialInteraction timepoint 4"
colnames(socialInteraction)[5] <- "socialInteraction timepoint 5"
boxplot(socialInteraction, main = "MOXFQsubscale socialInteraction", xlab = "timepoint", ylab = "MOXFQsub socialInteraction", col = "lightblue")

########
#
#  EQVAS
#
########


EQVAS_idx <- which(annoFile_here$group == "EQVAS" & annoFile_here$useIt == "x")
EQVAS_All <- data.frame(mydat_here[,EQVAS_idx])
myEQVASLastEntryInTime <- keepLastValueFromTimeResolvedResultDF(EQVAS_All)
myEQVASPreOpValue <- EQVAS_All[,1]

# Differenzen Prae and Post OP
myEQVAS_diff <- myEQVASLastEntryInTime - myEQVASPreOpValue

########
#
#  AOFAS
#
########
AOFAS_idx <- which(annoFile_here$group == "AOFAS" & annoFile_here$question == "sum")
AOFAS_All <- data.frame(mydat_here[,AOFAS_idx])
myAOFASLastEntryInTime <- keepLastValueFromTimeResolvedResultDF(AOFAS_All)
myAOFASPreOpValue <- AOFAS_All[,1]

# Differenzen Prae and Post OP
myAOFAS_diff <- myAOFASLastEntryInTime - myAOFASPreOpValue

# boxplot AOFAS over time
AOFAS_All <- cbind(AOFAS_All, rowMeans(AOFAS_All[,6:7], na.rm = TRUE))
colnames(AOFAS_All)[8] <- "AOFAS timepoint 6+7"
boxplot(AOFAS_All[,c(1:5,8)], main = "AOFAS", xlab = "timepoint", ylab = "AOFAS", col = "lightblue")


########
#
#  Radiolocical findings
#
########

# HV preoperative
HVpraeOp_idx <- which(annoFile_here$group == "radiological findings" & annoFile_here$question == "HV preoperative")
HVpraeOp_All <- data.frame(mydat_here[,HVpraeOp_idx])
HVpraeOpLastEntryInTime <- keepLastValueFromTimeResolvedResultDF(HVpraeOp_All)
sum(!is.na(HVpraeOpLastEntryInTime))

# HV postoperative
HVpostOp_idx <- which(annoFile_here$group == "radiological findings" & annoFile_here$question == "HV postoperative")
HVpostOp_All <- data.frame(mydat_here[,HVpostOp_idx])
HVpostOpLastEntryInTime <- keepLastValueFromTimeResolvedResultDF(HVpostOp_All)
sum(!is.na(HVpostOpLastEntryInTime))

# boxplot HV over time
HV_All <- cbind(rowMeans(HVpraeOp_All, na.rm = TRUE), HVpostOp_All[,c(2:7)])
colnames(HV_All)[1] <- "HV_preop"

HV_All <- cbind(HV_All, rowMeans(HV_All[,6:7], na.rm = TRUE))
colnames(HV_All)[8] <- "HV timepoint 6+7"

boxplot(HV_All[,c(1:5,8)], main = "HV", xlab = "timepoint", ylab = "HV", col = "lightblue")

# IMA preoperative
IMApraeOp_idx <- which(annoFile_here$group == "radiological findings" & annoFile_here$question == "IMA preoperative")
IMApraeOp_All <- data.frame(mydat_here[,IMApraeOp_idx])
IMApraeOpLastEntryInTime <- keepLastValueFromTimeResolvedResultDF(IMApraeOp_All)
sum(!is.na(IMApraeOpLastEntryInTime))

# IMA postoperative
IMApostOp_idx <- which(annoFile_here$group == "radiological findings" & annoFile_here$question == "IMA postoperative")
IMApostOp_All <- data.frame(mydat_here[,IMApostOp_idx])
IMApostOpLastEntryInTime <- keepLastValueFromTimeResolvedResultDF(IMApostOp_All)
sum(!is.na(IMApostOpLastEntryInTime))

# boxplot IMA over time
IMA_All <- cbind(rowMeans(IMApraeOp_All, na.rm = TRUE), IMApostOp_All[,c(2:7)])
colnames(IMA_All)[1] <- "IMA_preop"

IMA_All <- cbind(IMA_All, rowMeans(IMA_All[,6:7], na.rm = TRUE))
colnames(IMA_All)[8] <- "IMA timepoint 6+7"

boxplot(IMA_All[,c(1:5,8)], main = "IMA", xlab = "timepoint", ylab = "IMA", col = "lightblue")


# Differenzen Prae and Post OP
HV_diffLastEntryInTime <- HVpraeOpLastEntryInTime - HVpostOpLastEntryInTime
IMA_diffLastEntryInTime <- IMApraeOpLastEntryInTime - IMApostOpLastEntryInTime




########
#
#  Complications
#
########
# unique(annoFile_here$question)
# compl_idx <- which(annoFile_here$group == "complications" & annoFile_here$question == "another comlication")
# Compl_All <- data.frame(mydat_here[,compl_idx])
# ComplLastEntryInTime <- keepLastValueFromTimeResolvedResultDF(Compl_All)
# ComplLastEntryInTime <- dplyr::recode(ComplLastEntryInTime, "ja" = 1, "nein" = 0)
# sum(!is.na(ComplLastEntryInTime))
# all NAs .. due to <NA>


# Satisfaction
Satis_idx <- which(annoFile_here$group == "satisfaction" & annoFile_here$useIt == "x")
Satis_All <- data.frame(mydat_here[,Satis_idx])
SatisLastEntryInTime <- keepLastValueFromTimeResolvedResultDF(Satis_All)
sum(!is.na(SatisLastEntryInTime))

# recommendations
recomm_idx <- which(annoFile_here$group == "recommondation" & annoFile_here$useIt == "x")
recomm_All <- data.frame(mydat_here[,recomm_idx])
recommLastEntryInTime <- keepLastValueFromTimeResolvedResultDF(recomm_All)
sum(!is.na(recommLastEntryInTime))

# pain killers
painK_idx <- which(annoFile_here$group == "pain killers" & annoFile_here$useIt == "x" & annoFile_here$question == "Schmermitteleinnahme")
painK_All <- data.frame(mydat_here[,painK_idx])


myPainKillers <- matrix(nrow = nrow(painK_All), ncol = ncol(painK_All))
for (jj in 1:ncol(painK_All)) {
  myV <- painK_All %>% pull(jj)
  #print(table(myV))
  pk <- dplyr::recode(myV, "> 6 Wochen" =  5, "< 6 Wochen" =  4,"< 3 Wochen" = 3, "< 2 Wochen" = 2, "< 1 Woche" = 1, "0-3 Tage" = 0)
  myPainKillers[,jj] <- pk
}

# visualize painkiller
plot(1,1)
hist(myPainKillers)
#hist(myPainKillers, legend(x = "topright", legend = c("> 6 Wochen" =  5, "< 6 Wochen" =  4,"< 3 Wochen" = 3, "< 2 Wochen" = 2, "< 1 Woche" = 1, "0-3 Tage" = 0)))


painKMaxEntryInTime <- rowMaxs(myPainKillers, na.rm = TRUE)
painKMaxEntryInTime[!is.finite(painKMaxEntryInTime)] <- NA
sum(!is.na(painKMaxEntryInTime))
hist(painKMaxEntryInTime)

pdf("../output/output_analgesics.pdf")
# Plot histogram
hist(painKMaxEntryInTime, main="",xlab = "weeks", ylab = "Frequency")

hist_data3 <- hist(painKMaxEntryInTime, plot = FALSE)

# Add counts on top of the bars
for (i in 1:length(hist_data3$counts)) {
  text(x = hist_data3$mids[i],
       y = hist_data3$counts[i],
       labels = hist_data3$counts[i],
       pos = 3)
}
dev.off()

# General smoker
genSmoke_idx <- which(annoFile_here$group == "General" & annoFile_here$question == "smoking")
general_smoke <- as.vector(mydat_here[,genSmoke_idx])
general_smoke <- data.frame(general_smoke)
gsVec <- general_smoke %>% pull(1)
general_smoke <- dplyr::recode(gsVec, "Ja" =  1, "Nein" =  0)

# General other diseases
genotherDis_idx <- which(annoFile_here$group == "General" & annoFile_here$question == "other disease")
general_otherDis <- as.vector(mydat_here[,genotherDis_idx])
gsVec <- data.frame(general_otherDis) %>% pull(1)
general_otherDis <- dplyr::recode(gsVec, "ja" =  1, "nein" =  0)


# age & sex
age_idx <- which(annoFile_here$question == "age" & annoFile_here$useIt == "x")
age <- as.vector(mydat_here[,age_idx])
ageVec <- data.frame(age) %>% pull(1)

# sex
sex_idx <- which(annoFile_here$question == "sex" & annoFile_here$useIt == "x")
sex <- as.vector(mydat_here[,sex_idx])
sexVec <- data.frame(sex) %>% pull(1)
sexVec <- dplyr::recode(sexVec, "m" =  1, "f" =  0)


# find from myMOFXsumDF the last time point: clinical fu
followUpDurationClin <-  keepLastValueFromTimeResolvedResultDF(myMOFXsumDF)
dim(myMOFXsumDF)

myWeeks <- c(0, 6, 12, 26, 52, 104, 104)
Clin_FU <- vector(length = length(followUpDurationClin))

# loop over LastValueAndcheckColumn
for (i in 1:length(myWeeks)) {
  bool_hit <- myMOFXsumDF[,i] == followUpDurationClin
  Clin_FU[bool_hit] <- myWeeks[i]
}
# we have to fix that all NAs are zero, -> put NA back in
Clin_FU[Clin_FU == 0] <- NA


# find from HVpostOp_All the last time point: radiological fu
followUpDurationRad <-  keepLastValueFromTimeResolvedResultDF(HVpostOp_All)
dim(HVpostOp_All)
Rad_FU <- followUpDurationRad


pdf("../output/output_Rad_FU.pdf")
# Plot histogram
hist(Rad_FU, main="",xlab = "weeks", ylab = "Frequency")

hist_data4 <- hist(Rad_FU, plot = FALSE)

# Add counts on top of the bars
for (i in 1:length(hist_data4$counts)) {
  text(x = hist_data4$mids[i],
       y = hist_data4$counts[i],
       labels = hist_data4$counts[i],
       pos = 3)
}
dev.off()

# loop over LastValueAndcheckColumn
for (i in 1:length(myWeeks)) {
  bool_hit <- HVpostOp_All[,i] == followUpDurationRad
  bool_hit[is.na(bool_hit)] <- FALSE
  Rad_FU[bool_hit] <- myWeeks[i]
}
Rad_FU[Rad_FU == 0] <- NA
hist(Rad_FU)
which(Rad_FU == 26)

# get all with POST
table(annoFile_here$question)
posts_idx <- which(annoFile_here$question == "post" & annoFile_here$useIt == "x")
posts <- as.vector(mydat_here[,posts_idx])
postsVec <- data.frame(posts) %>% pull(1)

# get initials
initials_idx <- which(annoFile_here$question == "ID" & annoFile_here$useIt == "x")
initials <- as.vector(mydat_here[,initials_idx])
initialsVec <- data.frame(initials) %>% pull(1)



# bind all together!!

# generals
myGenerals <- cbind(initialsVec, postsVec, ageVec, sexVec, Clin_FU, Rad_FU, general_smoke, general_otherDis, painKMaxEntryInTime, recommLastEntryInTime, SatisLastEntryInTime)

# radiologicals and more
myPostOPValues <- cbind(myEQ5D5LLastEntryInTime, myFFILastEntryInTime, myMOFXLastEntryInTime, myEQVASLastEntryInTime, myAOFASLastEntryInTime)
myRadioVars <- cbind(HVpraeOpLastEntryInTime, HVpostOpLastEntryInTime, HV_diffLastEntryInTime, IMApraeOpLastEntryInTime, IMApostOpLastEntryInTime, IMA_diffLastEntryInTime)
myDiffVars <- cbind(myFFI_diff, myMOFX_diff, myEQVAS_diff, myAOFAS_diff)

# add subclasses
completeDF <- cbind(myGenerals, myPostOPValues, myRadioVars, myDiffVars, myFFIOverTime_disability_norm,myFFIOverTime_difficulty_norm, myFFIOverTime_pain_norm)

# get preOp values
myPreOpValues <- cbind(myEQ5D5LPreOpValue, myFFIPreOpValue, myMOFXPreOpValue, myEQVASPreOpValue, myAOFASPreOpValue)

# final DF
finalDF <- as.data.frame(cbind(myGenerals, myPreOpValues, myPostOPValues, myRadioVars, myDiffVars, myFFIOverTime_disability_norm,myFFIOverTime_difficulty_norm, myFFIOverTime_pain_norm, socialInteraction, walkingNstanding, painMOXFQ))
dim(finalDF)


# filtering!!
write.table(finalDF, file = "../output/output_Dataframe_finalDF_withSubclasses_20241212.txt", sep="\t", quote = FALSE)
# filter for Postoperativ
finalIncluded <- finalDF[!is.na(finalDF$initialsVec),]
write.table(finalIncluded, file = "../output/output_Dataframe_finalIncluded_allValuesBeforeMinMaxSummary.txt", sep="\t", quote = FALSE)

for (i in 3:ncol(finalIncluded)) {
  finalIncluded[,i] <- as.numeric(finalIncluded[,i])
}
str(finalIncluded)

# boxplotting HV pre-post
dim(finalDF)
colnames(finalDF)
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "HVp")]

tt <- t.test(as.numeric(ddf[,1]), as.numeric(ddf[,2]))
pV <- tt$p.value

pdf("../output/output_HV-pre-post.pdf")
boxplot(as.numeric(ddf[,1]), as.numeric(ddf[,2]))
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()


# boxplotting IMA pre-post
colnames(finalDF)
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "IMAp")]
colnames(ddf)

tt <- t.test(as.numeric(ddf[,1]), as.numeric(ddf[,2]))
pV <- tt$p.value

boxplot(as.numeric(ddf[,1]), as.numeric(ddf[,2]))
legend("topright", paste("pValue: ",round(pV,6)))


pdf("../output/output_IMA-pre-post.pdf")
boxplot(as.numeric(ddf[,1]), as.numeric(ddf[,2]))
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

# boxplotting DMAA pre-post
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "myAOFAS")]
colnames(ddf)

tt <- t.test(lastCOI2, lastCOI3)
pV <- tt$p.value

pdf("../output/output_DMAA-pre-post.pdf")
boxplot(lastCOI2, lastCOI3)
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

# boxplotting AOFAS pre-post
dim(finalDF)
colnames(finalDF)
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "myAOFAS")]
colnames(ddf)

tt <- t.test(myAOFASPreOpValue, myAOFASLastEntryInTime)
pV <- tt$p.value

pdf("../output/output_AOFAS-pre-post.pdf")
boxplot(myAOFASPreOpValue, myAOFASLastEntryInTime)
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

# boxplotting MOXFQ pre-post
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "myMOFX")]
colnames(ddf)

tt <- t.test(myMOFXPreOpValue, myMOFXLastEntryInTime)
pV <- tt$p.value

pdf("../output/output_MOXFQ-pre-post.pdf")
boxplot(myMOFXPreOpValue, myMOFXLastEntryInTime)
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

# boxplotting MOXFQsubSocial pre-post
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "social")]
colnames(ddf)

tt <- t.test(as.numeric(ddf[,1]), as.numeric(ddf[,6]))
pV <- tt$p.value

pdf("../output/output_MOXFQsubSocial_interaction-pre-post.pdf")
boxplot(as.numeric(ddf[,1]), as.numeric(ddf[,6]))
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

# boxplotting MOXFQsubwalkingNstanding pre-post
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "walking")]
colnames(ddf)

tt <- t.test(as.numeric(ddf[,1]), as.numeric(ddf[,6]))
pV <- tt$p.value

pdf("../output/output_MOXFQsubwalkingNstanding_interaction-pre-post.pdf")
boxplot(as.numeric(ddf[,1]), as.numeric(ddf[,6]))
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

# boxplotting MOXFQpain pre-post
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "painM")]
colnames(ddf)

tt <- t.test(as.numeric(ddf[,1]), as.numeric(ddf[,6]))
pV <- tt$p.value

pdf("../output/output_MOXFQsubPain_interaction-pre-post.pdf")
boxplot(as.numeric(ddf[,1]), as.numeric(ddf[,6]))
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

# boxplotting FFI pre-post
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "myFFI")]
colnames(ddf)

tt <- t.test(myFFIPreOpValue, myFFILastEntryInTime)
pV <- tt$p.value

pdf("../output/output_FFI-pre-post.pdf")
boxplot(myFFIPreOpValue, myFFILastEntryInTime)
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

# boxplotting FFIsubScales disability pre-post
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "disability")]
colnames(ddf)

tt <- t.test(as.numeric(ddf[,1]), as.numeric(ddf[,6]))
pV <- tt$p.value

pdf("../output/output_FFIsubScales disability_interaction-pre-post.pdf")
boxplot(as.numeric(ddf[,1]), as.numeric(ddf[,6]))
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

# boxplotting FFIsubScales difficulty pre-post
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "difficulty")]
colnames(ddf)

tt <- t.test(as.numeric(ddf[,1]), as.numeric(ddf[,6]))
pV <- tt$p.value

pdf("../output/output_FFIsubScales difficulty_interaction-pre-post.pdf")
boxplot(as.numeric(ddf[,1]), as.numeric(ddf[,6]))
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

# boxplotting FFIsubScales pain pre-post
dim(finalDF)
colnames(finalDF)
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "pain")]
colnames(ddf)

tt <- t.test(as.numeric(ddf[,2]), as.numeric(ddf[,7]))
pV <- tt$p.value

pdf("../output/output_FFIsubScales_pain_interaction-pre-post.pdf")
boxplot(as.numeric(ddf[,2]), as.numeric(ddf[,7]))
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

# boxplotting EQVAS pre-post
dim(finalDF)
colnames(finalDF)
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "myEQVAS")]
colnames(ddf)

tt <- t.test(myEQVASPreOpValue, myEQVASLastEntryInTime)
pV <- tt$p.value

pdf("../output/output_EQVAS-pre-post.pdf")
boxplot(myEQVASPreOpValue, myEQVASLastEntryInTime)
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

# boxplotting EQ5D5L pre-post
dim(finalDF)
colnames(finalDF)
ddf <- finalDF[,grep(x=colnames(finalDF), pattern = "myEQ5D5L")]
colnames(ddf)

tt <- t.test(myEQ5D5LPreOpValue, myEQ5D5LLastEntryInTime)
pV <- tt$p.value

pdf("../output/output_EQ5D5L-pre-post.pdf")
boxplot(myEQ5D5LPreOpValue, myEQ5D5LLastEntryInTime)
legend("topright", paste("pValue: ",round(pV,6)))
dev.off()

#
# save(finalIncluded, file = "FinalVariablesForRegression_2024-10-29.RData")

# generate summaries
colnames(finalIncluded[,3:ncol(finalIncluded)])

myMax <- apply(finalIncluded[,3:ncol(finalIncluded)],2,max,na.rm=TRUE)
myMin <- apply(finalIncluded[,3:ncol(finalIncluded)],2,min,na.rm=TRUE)
mySD <- apply(finalIncluded[,3:ncol(finalIncluded)],2,sd,na.rm=TRUE)
myMeans <- apply(finalIncluded[,3:ncol(finalIncluded)],2,mean,na.rm=TRUE)
myMedians <- apply(finalIncluded[,3:ncol(finalIncluded)],2,median,na.rm=TRUE)


myResults <- cbind(myMax, myMin, myMeans, myMedians, mySD)

# write out results
write.table(myResults, file = "../output/output_myResults_MICA_MeanSDnMore_allIncluded_wPreop_12-12-2024.txt", sep="\t", quote = FALSE)


# myAOFASPreOpValue --> is missing overall!
# final included without AOVAS preop
colnames(finalIncluded)
finalIncluded_noNAs <- finalIncluded[,-15]

colnames(finalIncluded)
finalIncluded_noPreOp <- finalIncluded[,-c(12:16)]
str(finalIncluded_noPreOp)

# general others is not numeric
colnames(finalIncluded)[7]
colnames(finalIncluded)[12:16]

# get rid of some (general others, )
finalIncluded_okForPlotting <- finalIncluded[,-c(7)]
colnames(finalIncluded_okForPlotting)



# Michi hier beginnen mit den linearen Regressionen
#load("FinalVariablesForRegression_2024-07-03.RData")

finalIncluded$ComplLastEntryInTime

pdf("../output/output_recommondation.pdf", 8,8)
# Plot histogram
hist(finalIncluded$recommLastEntryInTime, main="",xlab = "Likert scale", ylab = "Frequency")
hist_data5 <- hist(finalIncluded$recommLastEntryInTime, plot = FALSE)

# Add counts on top of the bars
for (i in 1:length(hist_data5$counts)) {
  text(x = hist_data5$mids[i],
       y = hist_data5$counts[i],
       labels = hist_data5$counts[i],
       pos = 3)
}
dev.off()

pdf("../output/output_satisfaction.pdf", 8,8)
# Plot histogram
hist(finalIncluded$SatisLastEntryInTime, main="",xlab = "Likert scale", ylab = "Frequency")
hist_data6 <- hist(finalIncluded$SatisLastEntryInTime, plot = FALSE)
# Add counts on top of the bars
for (i in 1:length(hist_data6$counts)) {
  text(x = hist_data6$mids[i],
       y = hist_data6$counts[i],
       labels = hist_data6$counts[i],
       pos = 3)
}
dev.off()



