prepPredPI = function(pisRes, folds, cellCounts, impute = NULL, pis = NULL, nPCs = 2, estPCs = TRUE, trainPCAs = NULL,
                      hypFrame = pisRes$hypFrame, what = "pis", loMat, dfSpic){
    colnames(cellCounts) = paste0("expr.", colnames(cellCounts))
    if(what=="pis"){
        pisRes$hypFrame = pisRes$hypFrame[folds,]
        #Feature selection through hypothesis tests
        if(is.null(pis)){
            pisRes = addWeightFunction(pisRes)
            lmmms = fitLMMs(pisRes, fixedVars = c("age", "Class"), verbose = FALSE)
            nnRes = getResults(lmmms, "nn", "Class")
            nnSig = rownames(nnRes)[nnRes[, "pAdj"]<sigLevel]
            nnPairRes = getResults(lmmms, "nnPair", "Class")
            nnPairSig = rownames(nnPairRes)[nnPairRes[, "pAdj"]<sigLevel & !is.na(nnPairRes[, "pAdj"])]
        } else {
            nnSig = pis$nnSig;nnPairSig = pis$nnPairSig;
        }
        nnPis = vapply(pisRes$hypFrame$pimRes, FUN.VALUE = double(length(nnSig)+ length(nnPairSig)), function(x){
            c(if(length(length(nnSig))) x$pointDists$nn[nnSig],
              if(length(length(nnPairSig))) vapply(nnPairSig, FUN.VALUE = 1, function(y) {
                  tmp = getGp(gp = y, x$pointDists$nnPair)
                  if(is.null(tmp)) NA else tmp
              }))
        })
        nnPis = t(matrix(nnPis, nrow = length(c(nnSig, nnPairSig))))
        colnames(nnPis) = c(nnSig, nnPairSig)
        nnPis[is.na(nnPis)] = 0.5
        df = data.frame("Class" = factor(pisRes$hypFrame$Class, ordered = TRUE, levels = c("mixed", 'compartimentalised')),
                        "pis" = nnPis, cellCounts[folds,,drop = FALSE], "age" = pisRes$hypFrame$age)
    } else if(what == "NE"){
        dfSpiatTNBC = data.frame(prop = NA,  pisRes$hypFrame[, c("age", "Class")])
        if(is.null(pis)){
            TNBCspiatLms = lapply(seq_len(nrow(loMat)), function(j){
                dfSpiatTNBC$prop = loMat[j,]
                fit = try(lm(prop ~ age + Class, data = dfSpiatTNBC[folds,]))
                if(inherits(fit, "try-error")) 1 else summary(fit)$coef["Classmixed", "Pr(>|t|)"]
            })
            pis = rownames(loMat)[p.adjust(TNBCspiatLms, method = "BH")<sigLevel & !is.na(TNBCspiatLms)]
        }
        loMat[is.na(loMat)] = 0 #Default values
        pisMat = t(loMat[unlist(pis),folds, drop =FALSE])
        colnames(pisMat) = paste0("pis.", colnames(pisMat))
        df = data.frame("Class" = factor(dfSpiatTNBC$Class[folds], ordered = TRUE, levels = c("mixed", 'compartimentalised')),
                        pisMat, cellCounts[folds,,drop = FALSE], "age" = dfSpiatTNBC$age[folds])
    } else if(what == "spicy"){
        # dfSpicTNBC = data.frame(u = NA,  pisRes$hypFrame[, c("age", "Class")])
        # if(is.null(pis)){
        #     TNBCspicyLms = lapply(seq_len(nrow(dfSpic)), function(j){
        #         dfSpicTNBC$u = dfSpic[,j]
        #         fit = try(lm(u ~ age + Class, data = dfSpicTNBC[folds,]))
        #         if(inherits(fit, "try-error")) 1 else summary(fit)$coef["Classmixed", "Pr(>|t|)"]
        #     })
        #     pis = names(dfSpic)[p.adjust(TNBCspicyLms, method = "BH")<sigLevel]
        # }
        pisMat = dfSpic[folds, , drop =FALSE]#unlist(pis)
        #Use all variables, cause none are significant
        pisMat[is.na(pisMat)] = 0 #Default values
        if(ncol(pisMat)){
            colnames(pisMat) = paste0("pis.", colnames(pisMat))
        }
        df = data.frame("Class" = factor(pisRes$hypFrame$Class[folds], ordered = TRUE, levels = c("mixed", 'compartimentalised')),
                        pisMat, cellCounts[folds,,drop = FALSE], "age" = pisRes$hypFrame$age[folds])
    }
    list("df" = df[complete.cases(df),,drop = FALSE])
}
fitPred = function(df, xNames, nFolds){
    cv.glmnet(y = factor(df$Class), x = data.matrix(df[, xNames]), family = "binomial",
              alpha = 1, nfolds = nFolds - 1, type.measure = "class")
}
expit = function(x) {exp(x)/(1+exp(x))}
