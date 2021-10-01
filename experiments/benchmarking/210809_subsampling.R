##################################
# subsampling the 1.3 mln data set
##################################

# 1. trying to load in the data at once with TENxBrainData ####

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("scuttle")

# libraries
  library(TENxBrainData)
  # library(scuttle)
  library(scran)
  library(scater)
  
# dirs
  data.dir <- "/home/ruths/Projects/13_FlandersAI/data"

# we should have the data as a SingleCellExperiment
  tenx <- TENxBrainData()
  tenx
  tenx <- tenx[,1:1000000]
  tenx

# we can just subsample the data set to 1M cells (should help)

# 2. a rudimentary processing pipeline (OCSA)

  # first step normalization
    lsf <- librarySizeFactors(tenx)
    sizeFactors(tenx) <- lsf
    tenx <- logNormCounts(tenx)

  # feature selection
    dec.tenx <- modelGeneVarByPoisson(tenx)
    top.tenx <- getTopHVGs(dec.tenx, prop=0.1)

  # dimensionality reduction
    set.seed(10000)
    tenx <- denoisePCA(tenx, subset.row=top.tenx, technical=dec.tenx)
    ncol(reducedDim(tenx, "PCA"))
    
  # clustering: kmeans
    set.seed(500)
    clust.kmeans <- kmeans(reducedDim(tenx, "PCA"), centers = 300, nstart = 5, iter.max = 12)
    table(clust.kmeans$cluster)
    summary(clust.kmeans$size)
    colLabels(tenx) <- clust.kmeans$cluster
    
  # the computational part is finished, just save it now  
    saveRDS(tenx, file.path(data.dir,"tenx.rds"))
    
  # doing the subsetting
    prop <- table(colLabels(tenx))/dim(tenx)[2]
  
    iter <- c(1:5)  
    size <- c(3000, 6000, 12000, 25000, 50000, 75000, 100000, 125000,
              250000, 500000, 750000)
    
    set.seed(210809)

    
    for (k in iter){ 
      
      for (j in size){
        print(paste0("iter ",k))
        print(paste0("size ",j))
        subsample <- data.frame("cell" = character(), 
                           "clust" = factor())
          for (i in 1:length(prop)){
            subs.cell <- sample(colnames(tenx[,colLabels(tenx) ==i]),size = round(prop[[i]]*j), replace = F)
            subs.clust <- i
            subs <- data.frame("cell" = subs.cell, 
                                    "clust" = subs.clust)
            subsample <-rbind(subs,subsample)
        }
        print(dim(subs))
        print(dim(subsample))
        write.csv(subsample,file.path(data.dir, paste0("subsample_it_",k,"_size_",j,".csv")))
        assign(paste0("subsample_it_",k,"_size_",j),subsample)

      }
    }

  # full sample 
    fullsample <- data.frame("cell" = colnames(tenx), 
                                                 "clust" = colLabels(tenx))  
    write.csv(fullsample,file.path(data.dir,"fullsample.csv"))
  
    
  