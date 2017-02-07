if (!interactive()) pdf('/Users/yan/Documents/Research/Wellcome/treeseq-inference/data/metrics_by_mutation_rate+metrics_vs_mutrate_simple.pdf',width=10, height=7)
data <- read.csv('/Users/yan/Documents/Research/Wellcome/treeseq-inference/data/metrics_by_mutation_rate_data.csv')

toolcols <- c('tsibiny'='blue','tsipoly'='cyan','fastARG'='red','Aweaver'='green','RentPls'='magenta')
metrics <- c('RFrooted','RFunrooted','SPRunrooted','pathunrooted')
makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop('alpha must be between 0 and 1')
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  return(apply(newColor, 2, .makeTransparent, alpha=alpha))
}
datamean <- aggregate(subset(data, select=-ARGweaver_iterations), list(data$mutation_rate, data$error_rate), mean, na.rm=TRUE)
error.rates <- sort(unique(data$error_rate))
layout(t(seq_along(error.rates)))
sapply(metrics, function(m) {
    colnames = paste(names(toolcols), m, sep='_')
    sapply(error.rates, function(error.rate) {
        d = subset(data, error_rate==error.rate)
        dm = subset(datamean, error_rate==error.rate)
        matplot(d$mutation_rate, d[, colnames], type='p', main=paste(m, 'metric: error', error.rate),
            col=makeTransparent(toolcols,0.1), 
            ylab='Distance between true and inferred trees',
            xlab='mutation rate',
            log='x', ylim = c(0,max(d[, colnames], na.rm=TRUE)),
            pch = ifelse(data$error_rate == error.rates[1],1,ifelse(data$error_rate == error.rates[2], 2, 4)))
        matlines(dm$mutation_rate, dm[, colnames], lty=which(error.rates==error.rate), col=toolcols)
        mtext(names(toolcols), 1, line=rev(seq(-1.2, by=-0.8, along.with=toolcols)), adj=0.05,
            cex=0.7, col=toolcols)
    })
})

if (!interactive()) dev.off()
