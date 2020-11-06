install.packages('devtools')
install.packages('sets')
require('devtools')
install_version("ggplot2", version = "3.3.2", repos = "http://cran.us.r-project.org")

install_version("ggpubr", version = "0.4.0", repos = "http://cran.us.r-project.org")
install_version("gplots", version = "3.1.0", repos = "http://cran.us.r-project.org")
install_version("ggthemes", version = "4.2.0", repos = "http://cran.us.r-project.org")
install_version("Seurat", version = "3.2.2", repos = "http://cran.us.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
