# Statistical Analysis for TnSeq data
# Combinatorial selection in environmental hosts drives the evolution of a human pathogen
#Jason M. Park1,2, Soma Ghosh1, Tamara J. O’Connor1,*
#1Department of Biological Chemistry, The Johns Hopkins University School of Medicine, Baltimore, Maryland, USA
#2Current Address: Veterinary Microbiology and Pathology, Washington State University, Pullman, Washington, USA 

# Required packages
#install.packages("gplots", dependencies = TRUE)
#install.packages("RColorBrewer", dependencies = TRUE)
#install.packages("sm", dependencies = TRUE)
#install.packages("fitdistrplus", dependencies = TRUE)
#install.packages("nls2", dependencies = TRUE)
#install.packages("gtools", dependencies = TRUE)
#install.packages("data.table", dependencies = TRUE)
#install.packages("BSDA", dependencies = TRUE)

getwd()
setwd("fitnesstables/")
allfiles <- list.files(pattern="*.csv")

##### Step 1: Import fitness tables #####

for(files in allfiles)
{
	name <- which(strsplit(files, "")[[1]]==".")
	assign(gsub(" ","",substr(files, 1, name-1)),
	read.csv(paste(files, sep=",")))
}
setwd("../")
pool_list <- list(p1r1, p1r2, p2r1, p2r2, p4r1, p4r2, p6r1, p6r2, p9r1, p9r2, p10r1, p10r2, p11r1, p11r2, p12r1, p12r2)

##### Step 2. Fit Gaussian curve to each pool replicate #####

ratio_logtrans_nozero <- lapply(pool_list, function(x) {log10(x$ratio[-which(x$ratio==0)])}) 

plot.multi.dens <- function(s)
{
	junk.x = NULL
    junk.y = NULL
    for(i in 1:length(s)) {
         junk.x = c(junk.x, density(s[[i]])$x)
         junk.y = c(junk.y, density(s[[i]])$y)
    }
    xr <- range(junk.x)
    yr <- range(junk.y)
    plot(density(s[[1]]), xlim = xr, ylim = yr, main = "")
    for(i in 1:length(s)) {
        lines(density(s[[i]]), xlim = xr, ylim = yr, col = i)
     }
}
jpeg(file="pool_distribution_b4_norm.jpg", width = 5*300, height = 5*300, res = 300, pointsize =8)
plot.multi.dens(ratio_logtrans_nozero)
legend("topright", legend=c("p1r1", "p1r2", "p2r1", "p2r2", "p4r1", "p4r2", "p6r1", "p6r2", "p9r1", "p9r2", "p10r1", "p10r2", "p11r1", "p11r2", "p12r1", "p12r2"), col=(1:length(ratio_logtrans_nozero)), lwd=2, lty=1)
dev.off()

library(nls2)
all_hists <- list()
all_fits <- list()
gauss_ratio <- list()
for(i in 1:length(pool_list))
{
	all_hists[[i]] <- hist(ratio_logtrans_nozero[[i]],breaks=seq(min(ratio_logtrans_nozero[[i]]) -0.02, max(ratio_logtrans_nozero[[i]]) +0.02,by=0.02),plot=FALSE)
	all_fits[[i]] <- nls2(y ~ k*exp(-1/2*(x-mu)^2/sigma^2),start=c(mu=0,sigma=4,k=12),data=data.frame(x=all_hists[[i]]$breaks[1:length(all_hists[[i]]$breaks)-1],y=all_hists[[i]]$counts), trace=TRUE)
	#jpeg(file=paste("fitting_gaussian_curve_pool_", i,".jpg"), width = 5*300, height = 5*300, res = 300, pointsize =8)
	#plot(all_hists[[i]]$breaks[1:length(all_hists[[i]]$breaks)-1],all_hists[[i]]$counts)
	#lines(all_hists[[i]]$breaks[1:length(all_hists[[i]]$breaks)-1],predict(all_fits[[i]]),lty=2,col="red",lwd=3)
	#dev.off()
	
##### Step 3. Normalize ratios #####
	
	gauss_ratio[[i]] <- 10^summary(all_fits[[i]])$parameters["mu","Estimate"]
	pool_list[[i]]$new_ratio <- pool_list[[i]]$ratio / gauss_ratio[[i]]
}

new_ratio_logtrans_nozero <- lapply(pool_list, function(x) {log10(x$new_ratio)}[-which(x$new_ratio==0)]) 
jpeg(file="pool_distribution_after_norm.jpg", width = 5*300, height = 5*300, res = 300, pointsize =8)
plot.multi.dens(new_ratio_logtrans_nozero)
legend("topright", legend=c("p1r1", "p1r2", "p2r1", "p2r2", "p4r1", "p4r2", "p6r1", "p6r2", "p9r1", "p9r2", "p10r1", "p10r2", "p11r1", "p11r2", "p12r1", "p12r2"), col=(1:length(new_ratio_logtrans_nozero)), lwd=2, lty=1)
dev.off()

##### Step 4. Calculate stats per gene #####
library(data.table)

setattr(pool_list, 'names' , c("p1r1", "p1r2", "p2r1", "p2r2", "p4r1", "p4r2", "p6r1", "p6r2", "p9r1", "p9r2", "p10r1", "p10r2", "p11r1", "p11r2", "p12r1", "p12r2"))  # needs to be changed acc to pool_list
all_pools <- rbindlist(pool_list, idcol = "pool")
#write.table(all_pools,"all_pools_collated.csv")
 
ins_detail <- within(all_pools, rm(chromosome,count_1,count_2,mt_freq_t1,mt_freq_t2,pop_freq_t1,pop_freq_t2,D,W))
gene_detail <- ins_detail[(grep("^lpg",ins_detail$gene))]
gene_detail_ordered <- gene_detail[order(gene_detail$gene),]
#write.table(gene_detail_ordered,"insertion_details_for_gene.csv")

library(plyr)
temp <- count(all_pools, "gene")
temp1 <- temp[order(as.numeric(gsub("lpg","",temp$gene))),]
temp2 <- rename(temp1, c("freq"="insertions"))
gene_stats <- temp2[-nrow(temp2),]

##### Step 5. Calculate MAD for each gene and a modified Z (Zi) score to identify outliers #####

dt <- data.table(all_pools)

clean_dt <- within(dt, rm(count_1,count_2,mt_freq_t1,mt_freq_t2,pop_freq_t1,pop_freq_t2,D,W)) 
temp3 <- clean_dt[,{median=median(new_ratio);dev=abs(new_ratio-median);mad=median(dev);Ziscore=0.6745*(dev/mad);list(ratio, new_ratio, median,dev,mad, Ziscore)},by=gene]
temp4 <- merge(gene_stats, temp3, "gene")
gene_ins_mad <- rename(temp4, c("V1"="ratio","V2"="new_ratio","V3"="median","V4"="dev","V5"="MAD", "V6"="Ziscore"))
#write.table(gene_ins_mad,"a_seqinf.csv")

##### Step 6. Remove outliers (Zi >3.5 && <-3.5)
# For cases, where the Zi score is NA, num and den are 0 and hence Zi can be safely put as 0. Where Zi score is INF, MAD is 0, but numerator is non-zero. In that case we will check the numerator and if num > 5, we will remove that value as an outlier.  
gene_ins_mad[is.na(gene_ins_mad)] <- 0
gene_ins_mad <- gene_ins_mad[!(gene_ins_mad$Ziscore==Inf & gene_ins_mad$dev > 5),]
gene_ins_mad[sapply(gene_ins_mad, is.infinite)] <- 0
#write.table(gene_ins_mad,"per_gene_stats.csv")
clean_set_df <- gene_ins_mad[(gene_ins_mad$Ziscore < 3.5 & gene_ins_mad$Ziscore > -3.5),]
write.table(clean_set_df,"per_gene_stats_outlier_removed.csv")

##### Step 7. Repeat calculations after removing outliers #####

new_temp <- count(clean_set_df,"gene")
new_temp1 <- new_temp[order(as.numeric(gsub("lpg","",new_temp$gene))),]
final_gene_insertions <- rename(new_temp1, c("freq"="insertions"))
fresh_dt <- data.table(clean_set_df)
gene_avg_std <- fresh_dt[,{average=mean(new_ratio);std_dev=sd(new_ratio);list(average,std_dev)}, by=gene]
final_gene_avg_std <- rename(gene_avg_std, c("V1"="average", "V2"="Sdev"))
final_gene_stats <- merge(final_gene_insertions,final_gene_avg_std,"gene")
final_gene_stats <- within(final_gene_stats, log_trans_ratio <- log10(average))

##### Step 8. Log10 convert the average ratios to examine by Gaussian curve to determine the deviation points and the 0.5 mean population ######

log_ratios <- final_gene_stats$log_trans_ratio[-which(final_gene_stats$average==0)]
#write.table(log_ratios,"new_ratio.csv")

hist_ratio <- hist(log_ratios,breaks=seq(min(log_ratios) -0.02,max(log_ratios) +0.02,by=0.02),plot=FALSE)
fit_gauss <- nls2(y ~ k*exp(-1/2*(x-mu)^2/sigma^2),start=c(mu=3,sigma=4,k=12),data=data.frame(x=hist_ratio$breaks[1:length(hist_ratio$breaks)-1],y=hist_ratio$counts), trace=TRUE)
jpeg(file="fitting gaussian to log transformed average ratios outliers removed.jpg", width = 5*300, height = 5*300, res = 300, pointsize =8)
plot(hist_ratio$breaks[1:length(hist_ratio$breaks)-1],hist_ratio$counts)
lines(hist_ratio$breaks[1:length(hist_ratio$breaks)-1],predict(fit_gauss),lty=2,col="red",lwd=3)
dev.off()
 
##### Step 9. Determine the 0.5 SD mean population #####


a <- 10^(summary(fit_gauss)$parameters["mu","Estimate"]+abs(summary(fit_gauss)$parameters["sigma","Estimate"])/2)
b <- 10^(summary(fit_gauss)$parameters["mu","Estimate"]-abs(summary(fit_gauss)$parameters["sigma","Estimate"])/2)
meanpop <- final_gene_stats[(final_gene_stats$average < a & final_gene_stats$average > b),]

#####Step 10 determine Z-score #####

final_gene_stats$zscore <- (final_gene_stats$log_trans_ratio - summary(fit_gauss)$parameters["mu","Estimate"]) /abs(summary(fit_gauss)$parameters["sigma","Estimate"])

#####Step 11. Calculate unpaired, Welch t-test, followed by Benjamini and Hochberg correction for multiple testing  #####

library(BSDA)

pop_mean = mean(meanpop$average)
pop_sd = sd(meanpop$average)

final_gene_stats$pval <- tsum.test(pop_mean,pop_sd,nrow(meanpop),final_gene_stats$average,final_gene_stats$Sdev,final_gene_stats$insertions,alternative="two.sided", var.equal=FALSE, conf.level=0.95, mu=0)$p.value

final_gene_stats$adpval <- p.adjust(final_gene_stats$pval, method="fdr", n=nrow(final_gene_stats))
write.table(file="Zscores_per_gene.csv", final_gene_stats)



