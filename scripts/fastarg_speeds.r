#script to analyse data from test_fastARG.py
d <- read.delim("treeseq-inference/test_files/fastARG_comparison.out")
vary.seqlen = subset(d, sample_size==1000)
vary.sampsz = subset(d, seq_length==1000000)

pdf("fastARG_sims.pdf", 7, 9)
layout(1:2)
plot(N_c_records_fastARG ~ seq_length, data=vary.seqlen, ylab="number of coal records", pch=".", col=rgb(0,0,1,0.1), las=2, cex.axis=0.5, main="Original (black) vs fastARG inference (blue)", cex.main=0.9)
mtext("(sample size = 1000)", cex=0.8,line=0.5)
points(N_c_records_orig ~ seq_length, data=vary.seqlen, pch=".", col=rgb(0,0,0,0.1))
sl.mean1 <- tapply(vary.seqlen$N_c_records_orig, list(s=vary.seqlen$seq_length), mean)
lines(as.numeric(names(sl.mean1)), sl.mean1)
sl.mean2 <- tapply(vary.seqlen$N_c_records_fastARG, list(s=vary.seqlen$seq_length), mean)
lines(as.numeric(names(sl.mean2)), sl.mean2, col=rgb(0,0,1))


plot(N_c_records_fastARG ~ sample_size, data=vary.sampsz, ylab="number of coal records", pch=".", col=rgb(0,0,1,0.1), las=2, cex.axis=0.5, main="Original (black) vs fastARG inference (blue)", cex.main=0.9)
mtext("(seq length = 1000000)", cex=0.8,line=0.5)
points(N_c_records_orig ~ sample_size, data=vary.sampsz, pch=".", col=rgb(0,0,0,0.1))
ss.mean1 <- tapply(vary.sampsz$N_c_records_orig, list(s= vary.sampsz$sample_size), mean)
lines(as.numeric(names(ss.mean1)), ss.mean1)
ss.mean2 <- tapply(vary.sampsz$N_c_records_fastARG, list(s=vary.sampsz$sample_size), mean)
lines(as.numeric(names(ss.mean2)), ss.mean2, col=rgb(0,0,1,0.1))


layout(1:2)
plot(fastARG_time ~ seq_length, data=vary.seqlen, ylab="time (seconds)", pch=".", col=rgb(0,0,1,0.1), las=2, cex.axis=0.5, main="fastARG inference time", cex.main=0.9)
mtext("(sample size = 1000)", cex=0.8,line=0.5)
sl.mean <- tapply(vary.seqlen$fastARG_time, list(s=vary.seqlen$seq_length), mean)
lines(as.numeric(names(sl.mean)), sl.mean, col=rgb(0,0,1))


plot(fastARG_time ~ sample_size, data=vary.sampsz, ylab="time (seconds)", pch=".", col=rgb(0,0,1,0.1), las=2, cex.axis=0.5, main="fastARG inference time", cex.main=0.9)
mtext("(seq length = 1000000)", cex=0.8,line=0.5)
ss.mean <- tapply(vary.sampsz$fastARG_time, list(s= vary.sampsz$sample_size), mean)
lines(as.numeric(names(ss.mean)), ss.mean, col=rgb(0,0,1))
dev.off()