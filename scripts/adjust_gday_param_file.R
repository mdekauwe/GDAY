# Change various params in a G'DAY input file.
# Author: Martin De Kauwe
# Date: 12.05.2016
# Email: mdekauwe@gmail.com

if (!require("ini")){
    install.packages("ini")
    library(ini)
}

in_fname <- "base_start.cfg"
out_fname <- "test.ini"
g <- read.ini(in_fname)

# Change control
g$control$ncycle <- "true"
g$control$modeljm <- 3
g$control$print_options <- "end"

g$params$vcmax <- "55.0"
g$control$jmax <- "110.0"


write.ini(g, out_fname)
