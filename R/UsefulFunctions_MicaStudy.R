# 2025-04
# source this file to have some useful functions at hand to analyze the data

# Functionize
dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1
  return(edist)
}

# adapted function
# DecodeQuestionForOnlyAFewQWithAnnoFile(annoFile = annoToDecodeQuestions, mydat = relevantDat)
DecodeQuestionForOnlyAFewQWithAnnoFile <- function(annoFile, mydat, useTag = "xx", keyWord = "yes") {
  # introduce grouping var
  annoFile$BetterGroup <- paste(annoFile$weeks, annoFile$group, annoFile$TimeCode, sep="")
  myWalk <- unique(annoFile$BetterGroup)

  myFullResults <- matrix(nrow = nrow(mydat))
  myFullResults <- myFullResults[,-1]
  # go in find relevant categories of Qs
  for (i in 1:length(myWalk)) {
    print(myWalk[i])
    myCatTable <- annoFile[annoFile$BetterGroup == myWalk[i],]
    # subset categories with multiple questions
    myQmatrix <- matrix(nrow = nrow(mydat))
    myQmatrix <- myQmatrix[,-1]
    for (j in 1:length(unique(myCatTable$question))) {
      myQTable <- myCatTable[myCatTable$question == unique(myCatTable$question)[j],]
      # go for inidividual questions
      #print(unique(myQTable$question)[j])
      headerofinterest <- myQTable$mainHeader[myQTable$question == unique(myQTable$question)]
      mydatcols <- mydat[,headerofinterest]

      # go in each question and evalute the scores
      myQresults <- vector()
      for (k in 1:nrow(mydatcols)){
        q_bool <- as.vector(mydatcols[k,] == keyWord)
        if (sum(is.na(q_bool)) == 0 && sum(q_bool) == 1) {
          # out of stoopidity there is no grading but category
          # myQresults[k] <- as.numeric(myQTable$grading[q_bool])
          myQresults[k] <- as.numeric(myQTable$category[q_bool])
        } else {
          myQresults[k] <- NA
        }
      }
      myQmatrix <- cbind(myQmatrix, myQresults)
    }
    Qsummary <- vector()
    Qresults <- vector()
    Qsummary <- rowSums(myQmatrix, na.rm = TRUE)
    # fix zeros for NA
    Qsummary[Qsummary == 0] <- NA
    # bind to final result
    Qresults <- cbind(Qsummary, myQmatrix)
    # write proper colnames
    colnames(Qresults) <- paste(myWalk[i],unique(myCatTable$question)[j],colnames(Qresults), sep="-")

    myFullResults <- cbind(myFullResults, Qresults)
  }

  return(myFullResults)
}

# without summary
# adapted function
DecodeQuestionForOnlyAFewQWithAnnoFileNoSummary <- function(annoFile, mydat, useTag = "xx", keyWord = "yes") {
  #NumQuestions <- unique(paste(myannoFile$question,myannoFile$category))
  # scoring

  # introduce grouping var
  annoFile$BetterGroup <- paste(annoFile$weeks, annoFile$group, annoFile$TimeCode, sep="")
  #unique(annoFile$BetterGroup)
  myWalk <- unique(annoFile$BetterGroup)

  #myResultScore <- vector(length=nrow(mydat))
  myFullResults <- matrix(nrow = nrow(mydat))
  myFullResults <- myFullResults[,-1]
  # go in find relevant categories of Qs
  for (i in 1:length(myWalk)) {
    print(myWalk[i])
    myCatTable <- annoFile[annoFile$BetterGroup == myWalk[i],]
    # subset categories with multiple questions
    myQmatrix <- matrix(nrow = nrow(mydat))
    myQmatrix <- myQmatrix[,-1]
    for (j in 1:length(unique(myCatTable$question))) {
      myQTable <- myCatTable[myCatTable$question == unique(myCatTable$question)[j],]
      # go for inidividual questions
      #print(unique(myQTable$question)[j])
      headerofinterest <- myQTable$mainHeader[myQTable$question == unique(myQTable$question)]
      mydatcols <- mydat[,headerofinterest]
      # go in each question and evalute the scores
      myQresults <- vector()
      for (k in 1:nrow(mydatcols)){
        q_bool <- as.vector(mydatcols[k,] == keyWord)
        if (sum(is.na(q_bool)) == 0 && sum(q_bool) == 1) {
          # out of stoopidity there is no grading but category
          # myQresults[k] <- as.numeric(myQTable$grading[q_bool])
          myQresults[k] <- as.numeric(myQTable$category[q_bool])
        } else {
          myQresults[k] <- NA
        }
      }
      myQmatrix <- cbind(myQmatrix, myQresults)
    }
    #Qsummary <- vector()
    Qresults <- vector()
    #Qsummary <- rowSums(myQmatrix, na.rm = TRUE)
    # fix zeros for NA
    #Qsummary[Qsummary == 0] <- NA
    # bind to final result
    Qresults <- myQmatrix
    # write proper colnames
    colnames(Qresults) <- paste(myWalk[i],unique(myCatTable$question)[j],colnames(Qresults), sep="-")

    myFullResults <- cbind(myFullResults, Qresults)
  }

  colnames(myFullResults) <- unique(paste(annoFile$weeks, annoFile$question, annoFile$TimeCode, sep=""))
  return(myFullResults)
}

# keep only the last entry if not na
keepLastValueFromTimeResolvedResultDF <- function(myTimeResolvedDF) {
  resultVector <- vector(length =nrow(myTimeResolvedDF))
  resultVector <- rep(NA, nrow(myTimeResolvedDF))
  for (i in 1:ncol(myTimeResolvedDF)) {
    for (k in 1:nrow(myTimeResolvedDF)) {
      if (!is.na(myTimeResolvedDF[k,i])) resultVector[k] <- myTimeResolvedDF[k,i]
    }
  }
  return(resultVector)
}

# keep only the last entry if not na
keepLastColumnIDxFromTimeResolvedResultDF <- function(myTimeResolvedDF) {
  resultVector <- vector(length =nrow(myTimeResolvedDF))
  resultVector <- rep(NA, nrow(myTimeResolvedDF))
  for (i in 1:ncol(myTimeResolvedDF)) {
    for (k in 1:nrow(myTimeResolvedDF)) {
      if (!is.na(myTimeResolvedDF[k,i])) resultVector[k] <- i
    }
  }
  return(resultVector)
}

# panel function for pairs plot
Michipanel.cor <- function(x, y, cex.cor = 0.8, method = "pearson", ...) {
  options(warn = -1)                   # Turn of warnings (e.g. tied ranks)
  usr <- par("usr"); on.exit(par(usr)) # Saves current "usr" and resets on exit
  par(usr = c(0, 1, 0, 1))             # Set plot size to 1 x 1
  r <- cor(x, y, method = method, use = "pair")               # correlation coef
  p <- cor.test(x, y, method = method)$p.val                  # p-value
  n <- sum(complete.cases(x, y))                              # How many data pairs
  txt <- format(r, digits = 3)                                # Format r-value
  txt1 <- format(p, digits = 3)                                 # Format p-value
  txt2 <- paste0("r= ", txt, '\n', "p= ", txt1, '\n', 'n= ', n) # Make panel text
  text(0.5, 0.5, txt2, cex = cex.cor, ...)                      # Place panel text
  options(warn = 0)                                             # Reset warning
}

