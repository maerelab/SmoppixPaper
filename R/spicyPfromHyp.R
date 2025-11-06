spicyPfromHyp = function(hyp, pppVars = c("day", "root", "section", "drs"),
                         fixedVars = "day", ranVars = "root", imageID = "drs",
                         to = getFeatures(hyp), from = getFeatures(hyp), cores = 1 ){
    seYang = makeSEfromHyp(hyp, pppVars = pppVars)
    spicyTime = spicy(seYang, condition = fixedVars, subject = ranVars,
                      imageID = "imageID", cellType = "gene", from = from,
                      to = to, cores = cores)
    spicyTime$pairwiseAssoc = lapply(spicyTime$pairwiseAssoc, function(x){
        names(x) = hypYang[, "drs", drop = TRUE];x
    })
    pValsSpicy = lapply(selfName(names(spicyTime$weights)), function(genePair){
        dfYang = data.frame("u" = spicyTime$pairwiseAssoc[[genePair]], "weight" = spicyTime$weights[[genePair]],
                            hypYang[, c("root", "day")])
        if(sum(!is.na(dfYang$u)) < 5){
            return(NULL)
        }
        dfYang$u[is.na(dfYang$u)] = 0 #Following spicyR
        modSpicy = try(lmerTest::lmer(u ~ day + (1|root), data = dfYang,
                                      weights = weight, contrasts = list("day" = "contr.sum")),
                       silent = TRUE)
        if(is(modSpicy, "try-error")) {
            return(matrix(NA, 2, 2, dimnames = list(c("(Intercept)", "day1"), c("Estimate", "Pr(>|t|)"))))
        } else summary(modSpicy)$coef[, c("Estimate", "Pr(>|t|)")]
    })
    return(pValsSpicy)
}
