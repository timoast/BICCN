library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)


# load all peaks
all.peakfiles <- list.files(path = "peaks/", full.names = TRUE, pattern = "*.bed")
all.peakfiles <- lapply(X = all.peakfiles, FUN = read.table, sep = "\t")
# filter high scoring peaks
all.peakfiles <- lapply(X = all.peakfiles, FUN = function(x) {
  x[x$V5 > 950, ]
})
all.peaks <- lapply(
  X = all.peakfiles,
  FUN = makeGRangesFromDataFrame,
  seqnames.field = "V1",
  start.field = "V2",
  end.field = "V3"
)
all.peaks <- do.call(what = c, args = all.peaks)

# unify peaks
unified.peaks <- reduce(x = all.peaks, drop.empty.ranges = TRUE)

# remove peaks over 3 kb
unified.peaks <- unified.peaks[width(unified.peaks) < 3000]

# remove peaks on nonstandard chromosomes
seqnames.keep <- head(seqnames(BSgenome.Mmusculus.UCSC.mm10), 21)
unified.peaks <- unified.peaks[seqnames(unified.peaks) %in% seqnames.keep]

# save unified peaks
write.table(
  x = as.data.frame(x = unified.peaks),
  file = "peaks/unified_peaks.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)