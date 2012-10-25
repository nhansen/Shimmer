is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
is.installed("statmod")
