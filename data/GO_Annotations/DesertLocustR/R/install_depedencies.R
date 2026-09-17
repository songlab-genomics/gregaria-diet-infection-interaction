# install_dependencies.R

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
BiocManager::install(version = "3.21")
}

required_packages <- c("GO.db", "AnnotationHub", "clusterProfiler",
                       "rtracklayer", "Biostrings", "data.table", "dplyr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}
