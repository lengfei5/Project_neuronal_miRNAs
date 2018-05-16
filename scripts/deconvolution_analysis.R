###########
### test know packages for deconvolution
###########
library('dtangle')
names(shen_orr_ex)
#library('dtangle.data')
#names(shen_orr_ex)
truth = shen_orr_ex$annotation$mixture

pure_samples <- lapply(1:3, function(i) {
  which(truth[, i] == 1)
})
names(pure_samples) = colnames(truth)
pure_samples

Y <- shen_orr_ex$data$log
Y[1:4,1:4]

marker_list = find_markers(Y,pure_samples,data_type="microarray-gene",marker_method='ratio')

lapply(marker_list$L,head)

q = .1
quantiles = lapply(marker_list$V,function(x)quantile(x,1-q))
K = length(pure_samples)
n_choose = sapply(1:K,function(i){max(which(marker_list$V[[i]] > quantiles[[i]]))})
n_choose

marks = marker_list$L
dc <- dtangle(Y,pure_samples,n_choose,data_type='microarray-gene',markers=marks)

phats <- dc$estimates
plot(truth,phats,xlab="Truth",ylab="Estimates",xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))

## test CellMix and DSA
require(CellMix)
x <- ExpressionMix('GSE19830', verbose=TRUE)

annotation(x)
# extract mixed samples
mix <- mixedSamples(x)
# load TIGER marker list
ml <- MarkerList('TIGER')
ml

names(ml)
basisnames(x)

ml <- ml[c('brain', 'liver', 'lung')]
summary(ml)

mlx <- convertIDs(ml, mix, verbose=TRUE)
summary(mlx)

profplot(mlx[,1:10], mix)

mlsc <- extractMarkers(mlx, expb(mix, 2), method='SCOREM', alpha=10^-12)
summary(mlsc)





