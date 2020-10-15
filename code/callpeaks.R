library(Signac)
library(GenomicRanges)

fragfiles <- list.files(
  path = "fragments/", pattern = "*.sort.bed.gz$", full.names = TRUE
)
fragfiles <- unname(sapply(fragfiles, normalizePath))


peakcalls <- sapply(fragfiles, CallPeaks, macs2.path = "/home/stuartt/miniconda3/envs/biccn/bin/macs2")
peakcalls <- sapply(peakcalls, function(x) {
  x[x$score > 100]
})

peakcalls <- Reduce(f = c, x = peakcalls)
peakcalls <- reduce(x = peakcalls)
peakcalls <- subsetByOverlaps(x = peakcalls, ranges = blacklist_mm10, invert = TRUE)
peakcalls <- keepStandardChromosomes(x = peakcalls, pruning.mode = "coarse")

write.table(
  x = as.data.frame(x = peakcalls),
  file = "peaks/unified_peaks.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)