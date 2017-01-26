if (!interactive()) pdf('/home/yan/treeseq-inference/data/basic_test+aw_vs_fa.pdf',width=10, height=7)
data <- read.csv('/home/yan/treeseq-inference/data/basic_test_data.csv')

datamean <- aggregate(subset(data, select=-ARGweaver_iterations), list(data$mutation_rate), mean)
toolcols <- c('fastARG'='red','Aweaver'='green')
metrics <- c('RFrooted','RFunrooted','wRFrooted','wRFunrooted','SPRunrooted','pathunrooted')
layout(matrix(1:6,2,3))
sapply(metrics, function(m) {
    colnames = paste(names(toolcols), m, sep='_')
    matplot(data$mutation_rate, data[, colnames], type='p', pch=c(1,2), col=toolcols, main=m, 
        log='x', ylim = c(0,max(data[, colnames], na.rm=TRUE)))
    matlines(datamean$mutation_rate, datamean[, colnames], type='l', lty=1, col=toolcols)
    mtext(names(toolcols), line=seq(-1.2, by=-0.8, along.with=toolcols), adj=0.95,
        cex=0.7, col=toolcols)
})

if (!interactive()) dev.off()
