# A wrapper function to build a spiat object
makeSpiatCellTypes = function(df, geneName = "gene", sample_id, varNames){
    colnames(df)[match(c("x", "y"), colnames(df))] = c("Cell.X.Position", "Cell.Y.Position")
    seObj = SpatialExperiment(assays = matrix(0, 1, nrow(df)),
                               colData = df[, c(union(varNames, sample_id), geneName, "Cell.X.Position", "Cell.Y.Position")],
                               spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position"), sample_id = sample_id)
    define_celltypes(seObj, categories = NULL, category_colname = geneName)
}
makeSpiatCellTypesHF = function(hypFrame, varNames, sample_id, ...){
    lapply(seq_len(nrow(hypFrame)), function(x){
        cat(x)
        df = data.frame(coords(hypFrame[x, "ppp", drop = TRUE]), marks(hypFrame[x, "ppp", drop = TRUE]),
                        hypFrame[x, varNames])
        makeSpiatCellTypes(df, varNames = varNames, sample_id = sample_id, ...)
    })
}