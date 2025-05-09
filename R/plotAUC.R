#Plot the AUC representation
plotAUC = function(ecdfObs,ecdfNull, ...){
    pSeq = seq(0, 1, length.out = 1e4)
    plot(pSeq, y<- ecdfNull(quantile(ecdfObs, probs = pSeq)), type = "l", asp=1, ...,
         ylab = expression(G[0]~"["~G^-1 ~"(p)]"), xlab = "p", xlim = c(0,1), ylim = c(0,1))
    polygon(c(0, pSeq, 1), c(0, y, 0), col = "grey")
    abline(0,1, lty = "dotdash")
}