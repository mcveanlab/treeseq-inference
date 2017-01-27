if (!interactive()) pdf('/home/yan/treeseq-inference/data/metrics_by_mutation_rate+metrics_vs_mutrate.pdf',width=10, height=7)
data <- read.csv('/home/yan/treeseq-inference/data/metrics_by_mutation_rate_data.csv')

toolcols <- c('tsibiny'='blue','tsipoly'='cyan','fastARG'='red','Aweaver'='green')
metrics <- c('RFrooted','RFunrooted','wRFrooted','wRFunrooted','SPRunrooted','pathunrooted')
error.rates <- unique(data$error_rate)
layout(matrix(1:6,2,3))
sapply(metrics, function(m) {
    colnames = paste(names(toolcols), m, sep='_')
    d <- subset(data, error_rate==error.rates[3])
    matplot(d$mutation_rate, d[, colnames], type='l', lty=3, col=toolcols, main=m, 
        log='x', ylim = c(0,max(d[, colnames], na.rm=TRUE)))
    
    d <- subset(data, error_rate==error.rates[2])
    matlines(d$mutation_rate, d[, colnames], type='l', lty=2, col=toolcols)
        
    d <- subset(data, error_rate==error.rates[1])
    matlines(d$mutation_rate, d[, colnames], type='l', lty=1, col=toolcols)
    mtext(names(toolcols), line=seq(-1.2, by=-0.8, along.with=toolcols), adj=0.95,
        cex=0.7, col=toolcols)
})

if (!interactive()) dev.off()
