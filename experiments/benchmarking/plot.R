library(patchwork)
library(ggplot2)
library(reshape2)
library(cowplot)
library(dplyr)
library(tidyr)

imec.cols <- c("#3F98BD", "#929497", "#90298D", "#36337D", "#1582BE", "#99BDE4", "#C778AD", "#52BDC2", "#3F98BD", "#2D6C85")
cols <- c("#ff7f0e","#922b21","#76448a","#2874a6","#148f77","#1d8348","#d68910","#ba4a00", "#3F98BD")
cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

palette <- c("red4", "darkslategray3", "dodgerblue1", "darkcyan",
               "gray79", "black", "skyblue2", "dodgerblue4",
               "purple4", "maroon", "chocolate1", "bisque3", "bisque",
               "seagreen4", "lightgreen", "skyblue4", "mediumpurple3",
               "palevioletred1", "lightsalmon4", "darkgoldenrod1")

cols <- c("seurat"="#F8766D", "scanpy"="#00BA38", "severo"="#619CFF")

if(F) {
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
}

X <- read.csv("comparison_umap.csv") %>% dplyr::filter(implementation!="jl32") %>%
    dplyr::mutate(implementation = factor(implementation, levels=c("R", "py", "jl"), labels=c("seurat", "scanpy", "severo")),
        dataset = factor(dataset,
                    levels=c("1M_gz", "210129_raw_BALPBMC", "l5_all", "pbmc_68k", "3k"),
                    labels=c("Brain cells (1.3M)", "Covid-19 (500k)", "Mouse Brain Atlas (120k)", "pbmc_68k"="PBMC (68k)", "3k"="PBMC (3k)"),
        ),
    ) %>% tidyr::drop_na(implementation, dataset)

Y <- X %>% dplyr::group_by(dataset, implementation) %>% dplyr::summarize(r=first(pearson))

p <- ggplot(X) +
    geom_boxplot(aes(x=cut, ymin=ymin, lower=lower, middle=middle, upper=upper, ymax=ymax, group=cut, color=implementation), stat="identity") +
    geom_text(aes(x=-Inf, y=Inf, label=paste0("r = ", round(r,1))), color="black", data=Y, hjust=-0.5, vjust=1.5) +
    scale_color_manual(values=cols) +
    facet_grid(implementation~dataset) +
    theme(legend.position = "none") +
    xlab("Binned pairwise distance in original space") +
    ylab("Pairwise distance in embedding")
ggsave(file="comparison_umap.pdf", plot=p, width=20, height=15)

X <- read.csv("comparison.csv", stringsAsFactors=T)
X <- X %>% dplyr::filter(dataset %in% c("1M_gz", "210129_raw_BALPBMC.h5ad", "3k", "l5_all", "pbmc_68k")) %>%
    dplyr::mutate(implementation = factor(implementation, levels=c("R", "py", "jl"), labels=c("seurat", "scanpy", "severo")),
                step = factor(step,
                            levels=c("CreateSeuratObject", "FindVariableFeatures", "NormalizeData", "ScaleData", "RunPCA", "FindNeighbors", "FindClusters", "RunUMAP", "FindAllMarkers"),
#                            labels=c("filter_counts", "find_variable_features", "normalize_cells", "scale_features", "embedding", "shared_nearest_neighbours", "cluster")
                ),
                dataset = factor(dataset,
                            levels=c("1M_gz", "210129_raw_BALPBMC", "l5_all", "pbmc_68k", "3k"),
                            labels=c("Brain cells (1.3M)", "Covid-19 (500k)", "Mouse Brain Atlas (120k)", "pbmc_68k"="PBMC (68k)", "3k"="PBMC (3k)"),
                ),
    ) %>% tidyr::drop_na(implementation, dataset) %>% dplyr::select(-it, -clusters, -size)

# add speedup
Z <- X %>% dplyr::group_by(dataset, implementation) %>% dplyr::summarize(totaltime=sum(time)) %>% dplyr::ungroup(implementation) %>%
                tidyr::pivot_wider(names_from=implementation, values_from=totaltime) %>%
                dplyr::mutate(ref=severo) %>%
                tidyr::pivot_longer(c(seurat, severo, scanpy), names_to="implementation", values_to="totaltime") %>%
                dplyr::mutate(speedup=totaltime/ref) %>%
                dplyr::ungroup()

p <- ggplot(X, aes(x=implementation, y=time, fill=step)) +
    geom_bar(position = "stack", stat="identity") +
#    scale_fill_manual(values=cols) +
    geom_text(aes(x=implementation, y=totaltime, fill=NULL, label=paste0(round(speedup,1),"x")), data=Z, position=position_dodge(0.9), vjust=-1) +
    facet_wrap(~dataset, scales="free_y", nrow=1) + ylab("Time (s)") + xlab(NULL) + theme(panel.background=element_rect("#F8f8f8"))

#Y <- X %>% dplyr::group_by(dataset, implementation) %>% dplyr::summarize(totaltime=sum(time)) %>% dplyr::ungroup()
#ggplot(Y, aes(x=dataset, y=totaltime, fill=implementation)) + geom_bar(position="dodge", stat="identity")

ggsave(file="comparison_datasets.pdf", plot=p, width=20, height=15)

X <- read.csv("comparison_metrics.csv") %>%
    dplyr::mutate(implementation = factor(implementation, levels=c("R", "py", "jl"), labels=c("seurat", "scanpy", "severo")),
        dataset = factor(dataset,
                    levels=c("1M_gz", "210129_raw_BALPBMC", "l5_all", "pbmc_68k", "3k"),
                    labels=c("Brain cells (1.3M)", "Covid-19 (500k)", "Mouse Brain Atlas (120k)", "pbmc_68k"="PBMC (68k)", "3k"="PBMC (3k)"),
        ),
    ) %>% tidyr::drop_na(implementation, dataset)

Y <- X %>% dplyr::group_by(dataset, size, implementation) %>% dplyr::summarize(ari=median(ari), jaccard=median(jaccard), peakmem=median(peakmem)) %>% dplyr::ungroup()
p <- ggplot(Y, aes(x=dataset, y=peakmem, group=implementation, fill=implementation)) +
    geom_col(position="dodge") +
    scale_color_manual(values=cols) +
    ylab("Peak memory (GB)") + xlab(NULL)
p <- ggplot(Y, aes(x=size, y=peakmem, group=implementation, color=implementation)) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) +
    geom_point() +
    scale_color_manual(values=cols) +
    ylab("Peak memory (GB)") + xlab("Number of cells")

Z <- Y %>% dplyr::filter(implementation != "seurat") %>% dplyr::mutate(implementation = factor(implementation, levels=c("scanpy", "severo")))
p1 <- ggplot(Z, aes(x=dataset, y=jaccard, group=implementation, fill=implementation)) +
        geom_col(position="dodge") +
        ylab("Jaccard index") + xlab(NULL)

p2 <- ggplot(Z, aes(x=dataset, y=ari, group=implementation, fill=implementation)) +
        geom_col(position="dodge") +
        ylab("Adjusted rand index") + xlab(NULL)
combined <- p1 + p2 & theme(legend.position = "bottom")
p <- combined + plot_layout(guides = "collect")
ggsave(file="comparison_metrics.pdf", plot=p, width=20, height=15)
