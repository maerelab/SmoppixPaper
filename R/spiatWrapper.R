#Wrapper the spiat enrichment calculations
spiatWrapper = function(spiatList, hypFrame, rangeFac = 20, nCores = 1,
                        gridFeat = expand.grid("gene1" = getFeatures(hypFrame), "gene2" = getFeatures(hypFrame), stringsAsFactors = FALSE),
                        featsAll = getFeatures(hypFrame)){
    # Reference: overall proportions
    emptyProp = double(length(featsAll));names(emptyProp) = featsAll
    refProp = vapply(hypFrame$tabObs, FUN.VALUE = emptyProp, function(x){
        emptyProp[names(x)] = x;emptyProp
    })
    Range = median(vapply(hypFrame$ppp, FUN.VALUE = double(1),
                          function(x) max(apply(coords(x), 2, function(y) diff(range(y))))))/rangeFac
    #Roughly following vignette
    spiatMat = simplify2array(mclapply(seq_along(spiatList), mc.preschedule = FALSE, mc.cores = nCores, function(j){cat(j)
        y = spiatList[[j]]
        avP = vapply(seq_len(nrow(gridFeat)), FUN.VALUE = double(1), function(i){
            average_percentage_of_cells_within_radius(spe_object = y,
                reference_celltype = gridFeat[i, "gene1"], target_celltype = gridFeat[i, "gene2"],
                radius = Range, feature_colname = "gene")
        })
        #number_of_cells_within_radius
        avP
    }))
    idUni <- gridFeat[, "gene1"] == gridFeat[, "gene2"]
    rownames(spiatMat)[idUni] = gridFeat[idUni, "gene1"]
    rownames(spiatMat)[!idUni] = apply(gridFeat[!idUni, ], 1, paste, collapse = "--")
    # Log-ratio
    rp = refProp[gridFeat[, "gene2"],] - idUni #Subtract event itself for univariate enrichment
    loMat = log(spiatMat/t((t(rp)*100/(colSums(refProp)-1))))
    #SPIAT returns percentages + always subtract element itself in denominator
    loMat[is.infinite(loMat)] = NA
    list("spiatMat" = spiatMat, "loMat" = loMat, "Range" = Range, "gridFeat" = gridFeat)
}