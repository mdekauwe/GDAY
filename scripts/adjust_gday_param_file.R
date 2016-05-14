# Change various params in a G'DAY input file.
# Author: Martin De Kauwe
# Date: 12.05.2016
# Email: mdekauwe@gmail.com
#
# e.g....
#
# in_fname <- "base_start.cfg"
# out_fname <- "test.ini"
#
# replacements <- list("ncycle" = "true",
#                     "modeljm" = "3",
#                     "print_options" = "end",
#                     "jmax" = "110.0",
#                     "vcmax" = "55.0")
#
# adjust_gday_params(in_fname, out_fname, replacements)

adjust_gday_params <- function(in_fname, out_fname, replacements) {

  if (!require("ini")){
      install.packages("ini")
      library(ini)
  }

  g <- read.ini(in_fname)

  for (key in names(replacements)) {

    match_git <- key %in% names(g$git)
    match_files <- key %in% names(g$files)
    match_params <- key %in% names(g$params)
    match_state <- key %in% names(g$state)
    match_control <- key %in% names(g$control)

    if (match_git) {
      g$git[key] <- replacements[key]
    } else if (match_files) {
      g$files[key] <- replacements[key]
    } else if (match_params) {
      g$params[key] <- replacements[key]
    } else if (match_state) {
      g$state[key] <- replacements[key]
    } else if (match_control) {
      g$control[key] <- replacements[key]
    }

  }
  write.ini(g, out_fname)
}
