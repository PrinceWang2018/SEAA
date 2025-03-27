.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")

  msg <- paste0("=========================================================================================
", pkgname, " version ", version, "
Project URL: https://github.com/PrinceWang2018/SEAA
Reference URL: https://github.com/PrinceWang2018/SEAA_reference
If you use it in published research, please cite:
Qin, J., Huang, T., Wang, Z. et al. Bud31-mediated alternative splicing
is required for spermatogonial stem cell self-renewal and differentiation. 
Cell Death Differ (2022). https://doi.org/10.1038/s41418-022-01057-1
If you have any question, please contact: wangzixiang@sdu.edu.cn
=========================================================================================
                       --Tomorrow will be more beautiful! ^_^--")
  base::packageStartupMessage(msg)
}
