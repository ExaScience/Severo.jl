library(rhdf5)
library(patchwork)
library(ggplot2)
library(reshape2)
library(cowplot)
library(dplyr)
library(tidyr)

X1 <- read.csv("nohvf.csv")
X2 <- read.csv("hvf.csv")
X <- rbind(cbind(hvf=F,X1), cbind(hvf=T,X2))

# correct for number of clusters found in FindAllMarkers
X$implementation <- as.factor(X$implementation)
X$step <- as.factor(X$step)
X[X$step == "FindAllMarkers","t"] <- X[X$step == "FindAllMarkers","t"] / X[X$step == "FindAllMarkers","clusters"] * 30

ols <- lm(t~size:implementation:step + implementation:step + I(size^2):implementation:step - 1, X)
summary(ols)

Z <- X %>% dplyr::group_by(size, implementation, step) %>%
			dplyr::summarize(mu=mean(t, na.rm=T), ci=2*sd(t, na.rm=T)) %>%
			dplyr::ungroup()
Z <- cbind(Z,predict=predict(ols, Z))

extrap <- expand.grid(implementation=levels(X$implementation), step=levels(X$step), size=c(1000000, 2000000, 5000000))
extrap$t <- predict(ols, extrap)

pd <- position_dodge(.3)
x_scale <- scale_x_continuous(breaks=c(3000, 25000, 50000, 100000, 250000, 500000, 750000, 1000000))
ggplot(Z, aes(x=size,y=mu, fill=implementation, color=implementation)) +
	geom_errorbar(aes(ymin=mu-ci, ymax=mu+ci), width=.1, position=pd) +
	geom_point(position=pd) +
	#geom_point(aes(x=size, y=t), data=extrap, shape=4, position=pd) +
	geom_line(aes(x=size, y=predict, group=implementation, color=implementation), position=pd) +
	xlab("Number of cells") +
	ylab("Time (s)") +
	x_scale +
	guides(x =  guide_axis(angle = 45)) +
	facet_wrap(~step, scales="free_y", ncol=3)

XX <- data.frame(X)
XX$size <- as.factor(XX$size)
ggplot(XX, aes(fill=step, x=size, y=t)) + geom_bar(position="fill", stat="identity") + facet_wrap(~implementation, scales="free_y", ncol=2)

Y <- X %>% dplyr::group_by(hvf, dataset, it, size, implementation) %>%
			dplyr::summarize(t=sum(t), clusters=first(clusters)) %>%
		dplyr::ungroup()

YS <- Y %>%
		dplyr::group_by(size, implementation) %>%
			dplyr::summarize(time=mean(t, na.rm=T), ci=2*sd(t, na.rm=T), clusters=median(clusters), count=sum(!is.na(t))) %>%
		dplyr::ungroup()

ols <- lm(t~size:implementation + implementation + I(size^2):implementation - 1, Y)
#ols <- lm(t~size:implementation + implementation - 1, Y)
YS <- cbind(YS, predict=predict(ols, YS))

ggplot(YS, aes(x=size,y=time, fill=implementation, color=implementation)) + geom_pointrange(aes(ymin=(time-ci), ymax=(time+ci))) + x_scale + guides(x =  guide_axis(angle = 45))

ggplot(YS, aes(x=size,y=t, fill=implementation, color=implementation)) +
	geom_point(aes(x=size, y=t), data=Y, position=pd) +
	geom_line(aes(x=size, y=predict, group=implementation, color=implementation), data=YS, position=pd) +
	xlab("Number of cells") +
	ylab("Time (s)") +
	x_scale +
	guides(x =  guide_axis(angle = 45))

ggplot(YS, aes(x=size,y=count, fill=implementation, color=implementation)) +
	geom_point(position=position_dodge(0.03)) +
	scale_y_continuous(breaks=unique(YS$count)) +
	scale_x_log10(breaks=unique(YS$size))

YY <- data.frame(Y)
YY$size <- as.factor(YY$size)
ggplot(YY) + geom_boxplot(aes(x=size, fill=implementation, y=clusters))

S <- Y %>% dplyr::group_by(hvf, dataset, it, size) %>% dplyr::select(-clusters) %>%
			 	tidyr::pivot_wider(names_from=implementation, values_from=t) %>%
					dplyr::transmute(R=R/jl_opt, py=py/jl_opt, opt=jl/jl_opt) %>%
				tidyr::pivot_longer(cols=c(R, py, opt), names_to="implementation", values_to="speedup") %>%
			dplyr::ungroup()

ggplot(S, aes(x=size,y=speedup, fill=implementation, color=implementation)) +
        geom_point(position=pd) +
        xlab("Number of cells") +
        ylab("Time (s)") +
        x_scale +
        guides(x =  guide_axis(angle = 45)) +
        scale_y_continuous(breaks=seq(0, 50, 5))

S <- X %>% dplyr::group_by(hvf, dataset, it, size, step) %>% dplyr::select(-clusters) %>%
			 	tidyr::pivot_wider(names_from=implementation, values_from=t) %>%
					dplyr::transmute(R=R/jl_opt, py=py/jl_opt, opt=jl/jl_opt) %>%
				tidyr::pivot_longer(cols=c(R, py, opt), names_to="implementation", values_to="speedup") %>%
			dplyr::ungroup()

ggplot(S, aes(x=size, y=speedup, fill=implementation, color=implementation)) +
	geom_point(position=pd) +
	xlab("Number of cells") +
	ylab("Time (s)") +
	x_scale +
	guides(x =  guide_axis(angle = 45)) +
	facet_wrap(~step, scales="free_y", ncol=3)

P <- X %>% dplyr::group_by(hvf, size, it, implementation) %>% dplyr::mutate(total=sum(t)) %>% tidyr::pivot_wider(names_from=step, values_from=t) %>% dplyr::transmute(percent=FindAllMarkers/total, total=total, markers=FindAllMarkers) %>% dplyr::ungroup()
ggplot(P, aes(x=size, y=percent, fill=implementation)) + geom_boxplot()

X1 <- read.csv("ari_nohvf.csv")
X2 <- read.csv("ari_hvf.csv")
XX <- rbind(cbind(hvf=F,X1), cbind(hvf=T,X2))
YY <- XX %>% dplyr::group_by(size, implementation) %>% dplyr::summarize(ari_mu=mean(ari), ari_se=sd(ari), ri_mu=mean(ri), ri_se=sd(ri), jaccard_mu=mean(jaccard), jaccard_se=sd((jaccard))) %>% dplyr::ungroup()
wrap_plots(
	ggplot(YY) + geom_pointrange(aes(x=size, y=ari_mu, ymin=(ari_mu-ari_se), ymax=(ari_mu+ari_se), color=implementation), size=0.5, position = position_dodge(width=.03)) + scale_x_log10(breaks=unique(YY$size)) + ylab("Adjusted Rand Index"),
	ggplot(YY) + geom_pointrange(aes(x=size, y=ri_mu, ymin=(ri_mu-ri_se), ymax=(ri_mu+ri_se), color=implementation), size=0.5, position = position_dodge(width=.03)) + scale_x_log10(breaks=unique(YY$size)) + ylab("Rand Index")
)
#ggplot(YY) + geom_pointrange(aes(x=size, y=jaccard_mu, ymin=(jaccard_mu-jaccard_se), ymax=(jaccard_mu+jaccard_se), color=implementation), size=0.5, position = position_dodge(width=.03)) + scale_x_log10(breaks=unique(YY$size))

ggplot(XX) + geom_boxplot() + geom_jitter(color="black", size=0.4, alpha=0.9) + theme_ipsum()


X1 <- read.csv("ari_nohvf.csv")
X2 <- read.csv("ari_hvf.csv")
XX <- rbind(cbind(hvf=F,X1), cbind(hvf=T,X2))
YY <- XX %>% dplyr::group_by(size, implementation) %>% dplyr::summarize(ari=median(ari), jaccard=median(jaccard)) %>% dplyr::ungroup()
ggplot(YY, aes(x=size, y=ari, group=implementation, color=implementation)) + geom_line()

#scp lynx:/data/thaber/1.3M_subsamples_processed_hvf/hvf.csv .
#scp lynx:/data/thaber/1.3M_subsamples_processed/nohvf.csv .
#scp lynx:/data/thaber/1.3M_subsamples_processed/ari_nohvf.csv .
#scp lynx:/data/thaber/1.3M_subsamples_processed_hvf/ari_hvf.csv .
