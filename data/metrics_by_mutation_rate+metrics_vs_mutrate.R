if (!interactive()) pdf('/Users/yan/Documents/Research/Wellcome/treeseq-inference/data/metrics_by_mutation_rate+metrics_vs_mutrate.pdf',width=10, height=7)
data <- read.csv('/Users/yan/Documents/Research/Wellcome/treeseq-inference/data/metrics_by_mutation_rate_data.csv')

toolcols <- c('tsibiny'='blue','tsipoly'='cyan','fastARG'='red','Aweaver'='green')
metrics <- c('RFrooted','RFunrooted','wRFrooted','wRFunrooted','SPRunrooted','pathunrooted')
error.rates <- unique(data$error_rate)
layout(matrix(1:6,2,3))
error.rates <- sort(unique(data$error_rate))
layout(matrix(1:6,2,3))
sapply(metrics, function(m) {
    colnames = paste(names(toolcols), m, sep='_')
    matplot(data$mutation_rate, data[, colnames], type='p', col=toolcols, main=paste(m, "metric"), 
        ylab='Distance between true and inferred trees', 
        xlab='mutation rate (err: dotted=0.1, dashed=0.01, solid=0.0)',
        log='x', ylim = c(0,max(data[, colnames], na.rm=TRUE)),
        pch = ifelse(data$error_rate == error.rates[1],1,ifelse(data$error_rate == error.rates[2], 2, 4)))
    d <- subset(datamean, error_rate==error.rates[1])
    matlines(d$mutation_rate, d[, colnames], lty=1, col=toolcols)
    d <- subset(datamean, error_rate==error.rates[2])
    matlines(d$mutation_rate, d[, colnames], lty=2, col=toolcols)
    d <- subset(datamean, error_rate==error.rates[3])
    matlines(d$mutation_rate, d[, colnames], type='l', lty=3, col=toolcols)

    mtext(names(toolcols), 1, line=rev(seq(-1.2, by=-0.8, along.with=toolcols)), adj=0.05,
        cex=0.7, col=toolcols)
})

if (!interactive()) dev.off()
