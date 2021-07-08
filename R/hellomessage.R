.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")

  msg <- paste0("=========================================================================================
", pkgname, " version ", version, "
Project URL: https://github.com/PrinceWang2018/SEAA
Usages: https://cran.r-project.org/web/packages/UCSCXenaTools/vignettes/USCSXenaTools.html
If you use it in published research, please cite:
Wang et al., (2021). SEAA: Splicing Efficiency Analysis and Annotation
=========================================================================================
                       --Tomorrow will be more beautiful! ^_^--")
  base::packageStartupMessage(msg)
}
