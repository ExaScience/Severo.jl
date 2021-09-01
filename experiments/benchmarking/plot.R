library(rhdf5)
library(patchwork)
library(ggplot2)
library(reshape2)
library(cowplot)
library(dplyr)
library(tidyr)

plot.loadings <- function(X, dims=1:6, n=20) {
    top <- function(x, n) names(head(sort(abs(x)), n))

    p <- lapply(dims,
        function(i) {
            features <- top(X[,i], n=n)
            data.plot <- data.frame(x=X[features, i], n=features)
            ggplot(data.plot, aes(x=x, y=n)) +
                geom_point(col='blue') + labs(y = NULL) + theme_cowplot()
        }
    )

    wrap_plots(p)
}

loadings <- h5read("/data/thaber/210129_raw_BALPBMC.h5", "loadings")
vf <- h5read("/data/thaber/210129_raw_BALPBMC.h5", "hvf")
rownames(loadings) <- vf

loadings2 <- h5read("/data/thaber/res.h5", "loadings_hvf")
hvf <- h5read("/data/thaber/res.h5", "hvf")
rownames(loadings2) <- hvf

p <- plot.loadings(loadings, dims=1:3)
q <- plot.loadings(loadings2, dims=1:3)
wrap_plots(p,q)

if(F) {
cos.dist <- function(X) {
    sim <- X / sqrt(rowSums(X * X))
    sim <- sim %*% t(sim)
    as.dist(sim)
}

fname <- "l5_all.loom.h5"
fname <- "210129_raw_BALPBMC.h5ad.h5"
orig <- h5read(fname, "coordinates")
umap_jl <- h5read(fname, "jl_umap")
umap_R <- h5read(fname, "R_umap")
J <- sample(1:nrow(orig), 10000)

d_orig <- cos.dist(orig[J,])
d_jl <- cos.dist(umap_jl[J,])
d_R <- cos.dist(umap_R[J,])

cuts <- cut(d_orig, 50)
boxplot(d_jl~cuts,outline=FALSE,xaxt="n")
r <- cor(d_jl, d_orig, method="pearson")
text(x=par("usr")[1],y=par("usr")[4],labels=paste("r=",signif(r,2),sep=""),pos=4,xpd=NA)
}

X <- read.csv("comparison.csv", stringsAsFactors=T)
levels(X$implementation) <- list(seurat="R", scanpy="py", severo="jl")
X <- X %>% dplyr::drop_na(implementation)
ggplot(X, aes(x=implementation, y=t, fill=step)) + geom_bar(position = "stack", stat="identity") + facet_wrap(~size, scales="free_y", nrow=1)

Y <- X %>% dplyr::group_by(dataset, size, implementation) %>% dplyr::summarize(totaltime=sum(t)) %>% dplyr::ungroup()
ggplot(Y, aes(x=size, y=totaltime, fill=implementation)) + geom_bar(position="dodge", stat="identity")
