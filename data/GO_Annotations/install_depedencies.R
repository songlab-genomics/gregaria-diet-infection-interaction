# install_dependencies.R

cran_packages <- c("data.table", "dplyr")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
BiocManager::install(version = "3.21")
}

bioc_packages <- c("GO.db", "AnnotationHub", "clusterProfiler","Biostrings",
                   "rtracklayer")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}
