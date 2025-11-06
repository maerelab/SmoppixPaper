getLargestRange = function(owin){
    max(diff(owin$xrange), diff(owin$yrange))
}
genAltPPP = function(densScaled, what, gpRange, mu, win, scale,
                     gpRangeScale, gpRangeScaleI, p2){
    switch(what,
           "Poisson" = rpoispp(densScaled),
           "Cauchy" = rCauchy(densScaled, scale = scale, mu = mu, win = win),
           "GaussPoisson" = rGaussPoisson(densScaled, gpRange*gpRangeScale, p2 = p2, win = win),
           "MaternCluster" = rMatClust(densScaled, scale = gpRange*gpRangeScale, mu = mu, win = win),
           "MaternInhibition" = rMaternI(densScaled, r = gpRange*gpRangeScaleI, win = win),
           #Increase gpRangeScaleI
           "Thomas" = rThomas(densScaled, scale = gpRange*gpRangeScale, mu = mu, win = erosion(win, 50))
    )
}
makeSEfromHyp = function(hyp, pppVars = c("day", "root", "section", "drs")){
    df = Reduce(rbind,lapply(seq_len(nrow(hyp)), function(x){
        data.frame(coords(hyp[x, "ppp", drop = TRUE]),
                   "gene" = marks(hyp[x, "ppp", drop = TRUE]), hyp[x, pppVars])
    }))
    colnames(df)[ncol(df)] = pppVars[length(pppVars)] = "imageID"
    SpatialExperiment(assays = matrix(0, 1, nrow(df)),
                      colData = df[, c(pppVars, "gene", "x", "y")],
                      spatialCoordsNames = c("x", "y"),
                      sample_id = df[, "imageID"])
}
#Cluster alternative scenario with parent and offspring pps
genAltPPPbi = function(densScaled, what, gpRange, mu, win, scale, gpRangeScale){
    ppp = switch(what,
         "Cauchy_bi" = rCauchy(densScaled, scale = scale, mu = mu, win = win, saveparents = TRUE),
         "MaternCluster_bi" = rMatClust(densScaled, scale = gpRange*gpRangeScale,
                                     mu = mu, win = win, saveparents = TRUE),
         "Thomas_bi" = rThomas(densScaled, scale = gpRange*gpRangeScale, mu = mu,
                            win = erosion(win, 50), saveparents = TRUE)
    )
    #Split into parent and offspring
    parentsId = match(seq_along(attr(ppp, "parents")$x), attr(ppp, "parentid"))
    #First occurrence is the parent
    ppps = split.ppp(ppp, factor(seq_len(npoints(ppp)) %in% parentsId,
                                 levels = c(TRUE, FALSE), labels = c("parent", "offspring")))
}