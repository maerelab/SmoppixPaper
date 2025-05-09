#Estimate the missclassification rate using training and test data
estMissClass = function(dat, idTest, cellCounts, nFolds, what = "pis", ...){
    trainDat = prepPredPI(dat, -idTest, nPCs = nPCs, cellCounts = cellCounts, what = what,...)
    Nams = setdiff(colnames(trainDat$df), c("Class"))
    testDat = prepPredPI(dat, idTest, cellCounts = cellCounts, what = what, ...,
                         pis = list("nnSig" = gsub("pis\\.", "", grep(value = TRUE, invert = TRUE, "\\.\\.", grep("pis\\.", value = TRUE, Nams))),
                                    "nnPairSig" = gsub("\\.\\.", "--",gsub("pis\\.", "", grep(value = TRUE, "\\.\\.", Nams)))))
    estAccuracy(trainDat, testDat, exprNames = paste0("expr.", colnames(cellCounts)), nFolds)
}
estAccuracy = function(trainDat, testDat, exprNames, nFolds, what = "pis"){
    fitAge = glm(data = trainDat$df, Class~age, family = "binomial")
    fitAgeExpr = fitPred(trainDat$df, c("age", exprNames), nFolds = nFolds)
    fitAgeExprPI = fitPred(trainDat$df, c("age", exprNames, grep(what, names(trainDat$df), value = TRUE)), nFolds = nFolds)
    fitAgePI = fitPred(trainDat$df, c("age", grep(what, names(trainDat$df), value = TRUE)), nFolds = nFolds)
    predAge = predict(fitAge, newdata = testDat$df, type = "response")
    predAgExpr = predict(fitAgeExpr, type = "response", newx = data.matrix(testDat$df[, rownames(coef(fitAgeExpr))[-1]]))
    predAgExprPI = predict(fitAgeExprPI, type = "response", newx = data.matrix(testDat$df[, rownames(coef(fitAgeExprPI))[-1]]))
    predPI = predict(fitAgePI, type = "response", newx = data.matrix(testDat$df[, rownames(coef(fitAgePI))[-1]]))
    predMat = cbind("age" = predAge, "ageExpr" = c(predAgExpr), "ageExprPI" = c(predAgExprPI),
                    "agePI" = c(predPI))#, "agePC" = c(predPC), "agePCexpr" = predPCexpr)
    predClass = matrix(NA, nrow(predMat), ncol(predMat));predClass[predMat>0.5] = "compartimentalised"
    predClass[predMat<0.5] = "mixed"
    mc = testDat$df$Class==predClass
    return(mc)
}
