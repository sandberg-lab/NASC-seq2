library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(parallel)
library(Rsamtools)
library(tictoc)

args <- commandArgs(trailingOnly=TRUE)

infile <- as.character(args[1])
outfile <- as.character(args[2])
numCPU <- as.numeric(args[3])
exclusionString <- as.character(args[4])
exclusion <- as.logical(args[5])

idxstats <- Rsamtools::idxstatsBam(infile)
idxstats <- idxstats[grep(exclusionString,idxstats$seqnames,invert=as.logical(exclusion)),]

#chunk_bcs = c('TCAGACAGAATCCGTCTTGG') ## For now, option to run this on subset of cells is not implemented

tic(msg='Conversion rates were calculated for all cells...')
taglist <- c("BC","UB","SC","TC")

rsamtools_reads <- mclapply(1:nrow(idxstats), function(x) {
  parms <- ScanBamParam(tag=taglist,
                        what="pos",
                        which = GRanges(seqnames = idxstats[x,"seqnames"], ranges = IRanges(1,idxstats[x,"seqlength"])))
  dat <- scanBam(file = infile, param = parms)
  dt <- data.table(RG = dat[[1]]$tag$BC, UB = dat[[1]]$tag$UB, SC = dat[[1]]$tag$SC, TC = dat[[1]]$tag$TC)
  return(dt)
}, mc.cores = numCPU)
rsamtools_reads <- rbindlist(rsamtools_reads, fill = TRUE, use.names = TRUE)

setDTthreads(numCPU)



rsamtools_reads[, c("gA","tC","aG","cT") := tstrsplit(SC, ';', keep = c(2,6,7,11))][
    , tC := as.numeric(substr(tC,3,nchar(tC)))][
    , aG := as.numeric(substr(aG,3,nchar(aG)))][
    , cT := as.numeric(substr(cT,3,nchar(cT)))][
    , gA := as.numeric(substr(gA,3,nchar(gA)))][
    , c('cov.a','cov.c','cov.g','cov.t') := tstrsplit(TC,';',keep=c(1,2,3,4))][
    , cov.t := as.numeric(substr(cov.t,2,nchar(cov.t)))][
    , cov.a := as.numeric(substr(cov.a,2,nchar(cov.a)))][
    , cov.c := as.numeric(substr(cov.c,2,nchar(cov.c)))][
    , cov.g := as.numeric(substr(cov.g,2,nchar(cov.g)))]

rsamtools_reads[,convRatePerCell.tC := sum(tC)/sum(cov.t),by=RG]
rsamtools_reads[,convRatePerCell.aG := sum(aG)/sum(cov.a),by=RG]
rsamtools_reads[,convRatePerCell.cT := sum(cT)/sum(cov.c),by=RG]
rsamtools_reads[,convRatePerCell.gA := sum(gA)/sum(cov.g),by=RG]

ldt <- dtplyr::lazy_dt(rsamtools_reads)

out <- ldt %>%
  select(RG,convRatePerCell.tC,convRatePerCell.aG,convRatePerCell.cT,convRatePerCell.gA) %>%
  distinct() %>%
  as.data.table()

saveRDS(rsamtools_reads,outfile)
toc()
