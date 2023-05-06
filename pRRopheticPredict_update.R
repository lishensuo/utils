pRRopheticPredict2 <- function(testMatrix, drug, tissueType="all", batchCorrect="eb", powerTransformPhenotype=TRUE, removeLowVaryingGenes=.2, minNumSamples=10, selection=-1, printOutput=TRUE, removeLowVaringGenesFrom="homogenizeData", dataset="cgp2014")
{
  cgpTrainData <- getCGPinfo(drug, tissueType, dataset) # get the IC50 and expression data for this drug/tissueType
  
  predictedPtype <- calcPhenotype2(cgpTrainData$trainDataOrd, cgpTrainData$ic50sOrd, testMatrix, batchCorrect=batchCorrect, powerTransformPhenotype=powerTransformPhenotype, removeLowVaryingGenes=removeLowVaryingGenes, minNumSamples=minNumSamples, selection=selection, printOutput=printOutput, removeLowVaringGenesFrom=removeLowVaringGenesFrom)
  
  return(predictedPtype)
  
}

calcPhenotype2 <- function(trainingExprData, trainingPtype, testExprData, batchCorrect="eb", powerTransformPhenotype=TRUE, removeLowVaryingGenes=.2, minNumSamples=10, selection=-1, printOutput=TRUE, removeLowVaringGenesFrom="homogenizeData")
{
  
  # check if the supplied data are of the correct classes
  # if(class(testExprData)[1] != "matrix") stop("ERROR: \"testExprData\" must be a matrix.");
  # if(class(trainingExprData)[1] != "matrix") stop("ERROR: \"trainingExprData\" must be a matrix.");
  # if(class(trainingPtype) != "numeric") stop("ERROR: \"trainingPtype\" must be a numeric vector.");
  if(ncol(trainingExprData) != length(trainingPtype)) stop("The training phenotype must be of the same length as the number of columns of the training expressin matrix.");
  
  # check if an adequate number of training and test samples have been supplied.
  if((ncol(trainingExprData) < minNumSamples) || (ncol(testExprData) < minNumSamples))
  {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }

  # Get the homogenized data
  homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)
  
  # Do variable selection if specified. By default we remove 20% of least varying genes from the homogenized dataset.
  # We can also remove the intersection of the lowest 20% from both training and test sets (for the gene ids remaining in the homogenized data)
  # Otherwise, keep all genes.
  # Check batchCorrect paramter
  if(!(removeLowVaringGenesFrom %in% c("homogenizeData", "rawData")))
  {
    stop("\"removeLowVaringGenesFrom\" must be one of \"homogenizeData\", \"rawData\"")
  }
  
  keepRows <- seq(1:nrow(homData$train)) # by default we will keep everything
  if(removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1) # if the proportion of things to keep is between 0 and 1
  {
    if(removeLowVaringGenesFrom == "homogenizeData") # if we are filtering based on the homoginized data
    {
      keepRows <- pRRophetic:::doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes=removeLowVaryingGenes)
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if(printOutput) cat(paste("\n", numberGenesRemoved, "low variabilty genes filtered."));
    }
    else if(removeLowVaringGenesFrom == "rawData") # if we are filtering based on the raw data, i.e. the intersection of the things filtered from both datasets.
    {
      evaluabeGenes <- rownames(homData$test) # pull the gene names from the genes remaining in the homoginized dataset
      keepRowsTrain <- pRRophetic:::doVariableSelection(trainingExprData[evaluabeGenes, ], removeLowVaryingGenes=removeLowVaryingGenes)
      keepRowsTest <- pRRophetic:::doVariableSelection(testExprData[evaluabeGenes, ], removeLowVaryingGenes=removeLowVaryingGenes)
      keepRows <- intersect(keepRowsTrain, keepRowsTest) # only keep things that are kept (i.e. variable) in both datasets
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if(printOutput) cat(paste("\n", numberGenesRemoved, "low variabilty genes filtered."));
    }
  }
  
  # PowerTranform phenotype if specified.
  offset = 0
  if(powerTransformPhenotype)
  {
    if(min(trainingPtype) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
    {
      offset <- -min(trainingPtype) + 1
      trainingPtype <- trainingPtype + offset
    }
    
    transForm <- powerTransform(trainingPtype)[[6]]
    trainingPtype <- trainingPtype^transForm
  }
  
  # create the Ridge Regression model on our training data
  if(printOutput) cat("\nFitting Ridge Regression model... ");
  trainFrame <- data.frame(Resp=trainingPtype, t(homData$train[keepRows, ]))
  rrModel <- linearRidge(Resp ~ ., data=trainFrame)
  if(printOutput) cat("Done\n\nCalculating predicted phenotype...");
  
  # calculate the relative contribution of each gene to the prediction
  # i might report these, I don't know if there's any point.
  totBeta <- sum(abs(coef(rrModel)))
  eachBeta <- abs(coef(rrModel))
  eachContribution <- eachBeta/totBeta
  
  # predict the new phenotype for the test data.
  # if there is a single test sample, there may be problems in predicting using the predict() function for the linearRidge package
  # This "if" statement provides a workaround
  if(class(homData$test)[1] == "numeric")
  {
    n <- names(homData$test)
    homData$test <- matrix(homData$test, ncol=1)
    rownames(homData$test) <- n
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata=rbind(testFrame, testFrame))[1]
  }
  else #predict for more than one test sample
  {
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata=testFrame)
  }
  
  # if the response variable was transformed, untransform it now...
  if(powerTransformPhenotype)
  {
    preds <- preds^(1/transForm)
    preds <- preds - offset
  }
  if(printOutput) cat("Done\n\n");
  
  return(preds)
}
