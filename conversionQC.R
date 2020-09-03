library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(parallel)
library(Rsamtools)
library(tictoc)

args <- commandArgs(trailingOnly=TRUE)

infile = args[1]
outfile = args[2]
numCPU = args[3]
exclusionString = args[4]
exclusion = args[5]

idxstats <- Rsamtools::idxstatsBam(infile)
idxstats <- idxstats[grep(exclusionString,idxstats$seqnames,invert=as.logical(exclusion),]

#chunk_bcs = c('TCAGACAGAATCCGTCTTGG') ## For now, option to run this on subset of cells is not implemented

tic(msg='Conversion rates were calculated for all cells...')
taglist <- c("BC","UB","SC","TC")

rsamtools_reads <- mclapply(1:nrow(idxstats), function(x) {
  parms <- ScanBamParam(tag=taglist,
                        what="pos",
                        which = GRanges(seqnames = idxstats[x,"seqnames"], ranges = IRanges(1,idxstats[x,"seqlength"])))
  dat <- scanBam(file = featfile, param = parms)
  dt <- data.table(RG = dat[[1]]$tag$BC, UB = dat[[1]]$tag$UB, SC = dat[[1]]$tag$SC, TC = dat[[1]]$tag$TC)
  return(dt)
}, mc.cores = numCPU)
rsamtools_reads <- rbindlist(rsamtools_reads, fill = TRUE, use.names = TRUE)

setDTthreads(numCPU)

rsamtools_reads[, c("tC","aG") := tstrsplit(SC, ';', keep = 6:7)][
                , tC := as.numeric(substr(tC,3,nchar(tC)))][
                , aG := as.numeric(substr(aG,3,nchar(aG)))][
                , c('cov.t','cov.a') := tstrsplit(TC,';',keep=c(4,1))][
                , cov.t := as.numeric(substr(cov.t,2,nchar(cov.t)))][
                , cov.a := as.numeric(substr(cov.a,2,nchar(cov.a)))]

rsamtools_reads[,convRatePerCell.tC := sum(tC)/sum(cov.t),by=RG]
rsamtools_reads[,convRatePerCell.aG := sum(aG)/sum(cov.a),by=RG]
saveRDS(rsamtools_reads,outfile)
toc()
