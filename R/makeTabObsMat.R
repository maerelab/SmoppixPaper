makeTabObsMat = function(hypFrame, rn = rownames(hypFrame)){
    mat = matrix(0, nrow(hypFrame), length(getFeatures(hypFrame)),
                 dimnames = list(rownames(hypFrame), getFeatures(hypFrame)))
    for(to in as.character(rownames(hypFrame))){
        mat[to, names(hypFrame[[to, "tabObs"]])] = hypFrame[[to, "tabObs"]]
    }
    rownames(mat) = rn
    return(mat)
}