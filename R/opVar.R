opVar <- function(O) {
   vv <- apply(O, 1, function(x) sum(x^2))
   vv
}