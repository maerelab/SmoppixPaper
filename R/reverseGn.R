#Reverse gene pair names for spiat
reverseGn = function(x){
    c(vapply(x, FUN.VALUE = character(2), function(y){
        z = strsplit(y, split = "--")[[1]]
        c(y, paste(z[2:1], collapse = "--"))
    }))
}