plotLoess = function(x, y, main, xlab = "Log10 sum of relative expressions"){
    ord = order(x)
    plot(log10(x[ord]), log10(y[ord]), xlab = xlab, ylab = "Log10 p-value",
         pch = 20, cex = 0.25, main = main, cex.main = 0.8)
    lines(x = log10(x[ord]), predict(loess(log10(y) ~ log10(x)))[ord], col = "blue")
}