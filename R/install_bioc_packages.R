install_bioc_packages <- function(packages) {
  # Check if BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("BiocManager is not installed. Installing BiocManager...\n")
    install.packages("BiocManager")
  } else {
    cat("BiocManager is already installed.\n")
  }

  # Load BiocManager
  library(BiocManager)

  # Install the specified Bioconductor packages
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(paste("Installing Bioconductor package:", pkg, "...\n"))
      BiocManager::install(pkg)
    } else {
      cat(paste("Bioconductor package", pkg, "is already installed.\n"))
    }
  }
}


Stats_COVID = read.csv("~/Downloads/Stats_COVID.csv", check.names = FALSE)
sampled_data <- Stats_COVID[sample(nrow(Stats_COVID), 50), ]
y = colnames(sampled_data)[grep("RT_min:RT_max", colnames(sampled_data))]
y = unlist(lapply(y, function(x) strsplit(x, "RT_min:RT_max_")[[1]][2]))
#y =  unlist(lapply(y, function(x) paste(strsplit(x," ")[[1]], collapse = "_")))
colnames(sampled_data)[grep("RT_min:RT_max", colnames(sampled_data))] = y
c = colnames(sampled_data)[which(colnames(sampled_data) %in% y)]
df_c = sampled_data[, c]
df_c[df_c== "NIL"] <- "0.0:0.0"
df_c_mod <- data.frame(matrix(ncol = 2*length(c), nrow = dim(sampled_data)[1]))
x <- c(paste("RT_min",c,sep="_"), paste("RT_max",c,sep="_"))
colnames(df_c_mod) <- x
#colnames(df_c_mod) <- x
for(i in 1:length(c)){
  g = colnames(df_c)[i]
  df_c_mod[,paste("RT_min", g, sep="_")] = as.numeric(as.character(unlist(lapply(df_c[,i],
                                                                               function(x) strsplit(x,":")[[1]][1]))))
  df_c_mod[,paste("RT_max",g, sep="_")] = as.numeric(as.character(unlist(lapply(df_c[,i],
                                                                               function(x) strsplit(x,":")[[1]][2]))))
}

dim(df_c_mod)

stats = data.frame(ID=rep("",50),mz = rep(0,50), RT=rep(0,50))
stats$ID = sampled_data$ID
stats$mz = sampled_data$mz
stats$RT = sampled_data$RT_median
write.table(stats, "~/Desktop/test_package_dures/test_1/Stats.txt", col.names = TRUE, row.names = FALSE, sep="\t", quote = FALSE)
write.table(df_c_mod, "~/Desktop/test_package_dures/test_1/Sample-and_feature-wise-RT-tolerance.txt", col.names = TRUE, row.names = FALSE, sep="\t", quote = FALSE)


Stats_virg = read.csv("~/Desktop/test_package_dures/test_2/COVID_ST002016_POS_Targeted_Stats.csv", check.names = FALSE)
sampled_data <- Stats_virg[sample(nrow(Stats_virg), 50), ]
stats = data.frame(ID=rep("",50),mz = rep(0,50), RT=rep(0,50))
stats$ID = sampled_data$ID
stats$mz = sampled_data$mz
stats$RT = sampled_data$RT_median
head(stats)
write.table(stats, "~/Desktop/test_package_dures/test_2/Stats.txt", col.names = TRUE, row.names = FALSE, sep="\t", quote = FALSE)

Stats_virg_qc = read.csv("~/Desktop/test_package_dures/test_2/COVID_ST002016_PooledQC_POS_5ppm_5ppm_Untargeted_Stats.csv", check.names = FALSE)
sampled_data_qc <- Stats_virg_qc[sample(nrow(Stats_virg_qc), 50), ]
stats_qc = data.frame(ID=rep("",50),mz = rep(0,50), RT=rep(0,50))
stats_qc$ID = sampled_data_qc$ID
stats_qc$mz = sampled_data_qc$mz
stats_qc$RT = sampled_data_qc$RT_median
stats_both_ms1_ms2 = rbind(stats, stats_qc)
write.table(stats_both_ms1_ms2, "~/Desktop/test_package_dures/test_2/Stats.txt", col.names = TRUE, row.names = FALSE, sep="\t", quote = FALSE)

head(sampled_data)
y = colnames(sampled_data)[grep("RT_min:RT_max", colnames(sampled_data))]
y = unlist(lapply(y, function(x) strsplit(x, "RT_min:RT_max_")[[1]][2]))
#y =  unlist(lapply(y, function(x) paste(strsplit(x," ")[[1]], collapse = "_")))
colnames(sampled_data)[grep("RT_min:RT_max", colnames(sampled_data))] = y
c = colnames(sampled_data)[which(colnames(sampled_data) %in% y)]
df_c = sampled_data[, c]
df_c[df_c== "NIL"] <- "0.0:0.0"
df_c_mod <- data.frame(matrix(ncol = 2*length(c), nrow = dim(sampled_data)[1]))
x <- c(paste(c,"RT_min",sep="_"), paste(c,"RT_max",sep="_"))
colnames(df_c_mod) <- x
#colnames(df_c_mod) <- x
for(i in 1:length(c)){
  g = colnames(df_c)[i]
  df_c_mod[,paste(g,"RT_min",sep="_")] = as.numeric(as.character(unlist(lapply(df_c[,i],
                                                                               function(x) strsplit(x,":")[[1]][1]))))
  df_c_mod[,paste(g,"RT_max",sep="_")] = as.numeric(as.character(unlist(lapply(df_c[,i],
                                                                               function(x) strsplit(x,":")[[1]][2]))))
}
df_c_mod$ID = sampled_data$ID
df_c_mod$ID = sampled_data$ID
df_c_mod$ID = sampled_data$ID


df_c_mod_qc = df_c_mod


stats_temp = read.delim("~/Desktop/test_package_dures/test_1/Sample-and_feature-wise-RT-tolerance.txt")
dt = colnames(stats_temp)
dt_m = paste("RT_min_", unlist(lapply(dt, function(x) strsplit(x, "_RT")[[1]][1]))[1:144],sep="")
dt_p = paste("RT_max_", unlist(lapply(dt, function(x) strsplit(x, "_RT")[[1]][1]))[145:288],sep="")
colnames(stats_temp) = c(dt_m,dt_p)
cleaned_strings <- gsub("X", "", colnames(stats_temp))
colnames(stats_temp) = cleaned_strings
write.table(stats_temp, "~/Desktop/test_package_dures/test_1/Sample-and_feature-wise-RT-tolerance.txt", col.names = TRUE, row.names = F, sep="\t", quote = F)
