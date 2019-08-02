library(Seurat)

load('.RData')

# Downsample by UMI
Down_Sample_Matrix <- function (expr_mat) {
  min_lib_size <- min(colSums(expr_mat))
  down_sample <- function(x) {
    prob <- min_lib_size/sum(x)
    return(unlist(lapply(x, function(y) {
      print(sprintf("%s %s", y, prob))
      rbinom(1, y, prob)
    })))
  }
  down_sampled_mat <- apply(expr_mat, 2, down_sample)
  return(down_sampled_mat)
}

# idents <- seu@active.ident
# all_ages <- factor(sapply(strsplit(as.character(seu@meta.data$orig.ident), '_'), function(x) x[1]))

library(parallel)

load("data2/all_ages.RData")
load("data2/genes.RData")
load("data2/ident.Rdata")
load("data2/raw_counts.RData")
all_ages <- all_cels_by_age
idents <- ident
raw_counts <- counts
rm(all_cels_by_age)
rm(ident)
rm(counts)

cl <- makePSOCKcluster(2)
# clusterCall(cl, function() load("data2/raw_counts.RData"))
clusterExport(cl, c('raw_counts',
                    'idents',
                    'all_ages',
                    'Down_Sample_Matrix',
                    'nv.genes.labels'))

levels_1_10 <- levels(idents)[1:10]

noise_1_10 <- lapply(levels_1_10, function(l){
  cl_id <- idents == l
  ages <- all_ages[cl_id]
  # ages_freq <- table(ages)
  # cl_id <- rep(FALSE, length(cl_id))
  # cl_id[c(sample(which(ages == '2-3mo'), min(ages_freq['2-3mo'], 400)), sample(which(ages == '21-22mo'), min(ages_freq['21-22mo'], 400)))] <- TRUE
  # ages <- all_ages[cl_id]
  counts <- raw_counts[,cl_id]

  ages_freq <- table(ages)
  if(any(ages_freq < 10)) return(NULL)

  # zeros <- which(Matrix::rowSums(counts) == 0)
  # counts <- data.matrix(counts[-zeros,])
  # print(sum(is.na(counts)))

  # counts <- Down_Sample_Matrix(counts)

  do.call(rbind, lapply(levels(ages), function(a) {
    counts_age <- counts[,ages == a]
    counts_age <- counts_age[,sample(1:ncol(counts_age), min(ages_freq))]

    nv.genes.labels.present <- nv.genes.labels[nv.genes.labels %in% rownames(counts_age)]

    counts_age <- counts_age[nv.genes.labels.present,]
    counts_age <- sqrt(counts_age)

    center <- apply(counts_age, 1, mean)
    noise <- apply(counts_age, 2, function(x){
      sqrt(sum((x-center)^2))
      # sqrt(sum((x-center)^2)/(length(x) - 1))/mean(x)
    })
    data.frame(noise = noise, cluster = l, age = a)
  }))
})

noise_data_1_10 <- do.call(rbind, noise_1_10)

# ----

genes.ordered <- order(seu@assays$SCT@meta.features$sct.gmean, decreasing = TRUE)
bins <- lapply(0:9, function(x){
  start=floor(x/10*length(genes.ordered))+1
  end=floor((x+1)/10*length(genes.ordered))
  genes.ordered[start:end]
})
bins<- bins[2:9]
bins<-lapply(bins, function(bin){
  bin[order(seu@assays$SCT@meta.features$sct.variance[bin])]
})

nv.genes <- lapply(bins, function(bin){
  bin[1:(length(bin)%/%10)]
})
nv.genes <- unlist(nv.genes)

seu <- seu[nv.genes,]
nv.genes.labels <- rownames(seu@assays$SCT@counts)

# noise <- do.call(rbind, lapply(levels(seu@active.ident), function(l){
#   cnt <- table(seu@meta.data$age.ident[seu@active.ident == l])
#   cond <- all(cnt >= 10)
#   if (cond) {
#     do.call(rbind, lapply(levels(seu@meta.data$age.ident), function(a) {
#       counts <- sqrt(seu[,seu@active.ident == l & seu@meta.data$age.ident == a]@assays$RNA@counts)
#       counts <- counts[,sample(1:ncol(counts), min(cnt))]
#       center <- apply(counts, 1, mean)
#       noise <- apply(counts, 2, function(x){
#         sqrt(sum((x-center)^2))
#       })
#       data.frame(noise = noise, cluster = l, age = a)
#     }))
#   } else {
#     NULL
#   }
# }))

library(ggplot2)

ggplot(noise_data_1_10, aes(factor(cluster), noise, fill = age)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = rev(scales::hue_pal()(2))) +
  theme_bw()

library(dplyr)

noise_data_1_10 %>%
  group_by(cluster) %>%
  do({
    data.frame(p = wilcox.test(noise ~ age, .)$p.value)
  }) %>%
  mutate(padj = p.adjust(p, method = "fdr")) %>%
  write.csv("noize2.csv")

save.image('work2.RData')
