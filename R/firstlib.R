.First.lib <- function(lib, pkg) {
   library.dynam("operator", pkg, lib)
   cat("operator 1.0 loaded\n")
}

