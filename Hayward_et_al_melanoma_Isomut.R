library(tidyverse)
library(BSgenome.Hsapiens.NCBI.GRCh38)

snv <- read.delim("20200629_hayward5_isomut_output/all_SNVs_after3rdround.isomut", stringsAsFactors = FALSE)

snv$Supp <- round(snv$Cov * snv$AF)

# Filter for at least 3 supporting reads, and 
# AF > .15: most of the remaining  events in normal samples are below this, 
# while the tumor mutationss show a roughly normal distribution around .3-.4
filter(snv, grepl("MELA", snv$Sample), round(Cov*AF) > 3, Cleanliness > 0.99) -> snv_tumor1
arrange(snv_tumor1, Sample, Chr, Pos) %>% 
  mutate(prevmut = Pos - lag(Pos),
         nextmut = lead(Pos) - Pos,
         prevmut = ifelse(prevmut <= 0, NA, prevmut),
         nextmut = ifelse(nextmut <= 0, NA, nextmut), muttype = paste(Ref, ">", Alt)) -> snv_tumor1


with(snv_tumor1, GRanges(seqnames = Chr, IRanges(start = Pos-1, end = Pos+1))) %>% 
  getSeq(Hsapiens, .) %>%
  as.character() -> snv_tumor1$tripmut

snv_tumor1$tripmut96 <- snv_tumor1$tripmut
snv_tumor1[which(snv_tumor1$Ref %in% c("A", "G")), "tripmut96"] <- as.character(reverseComplement(DNAStringSet(snv_tumor1[which(snv_tumor1$Ref %in% c("A", "G")), "tripmut96"])))
snv_tumor1$Alt2 <- snv_tumor1$Alt
snv_tumor1[which(snv_tumor1$Ref %in% c("A", "G")), "Alt2"] <- as.character(reverseComplement(DNAStringSet(snv_tumor1[which(snv_tumor1$Ref %in% c("A", "G")), "Alt2"])))

snv_tumor1 %>% 
  apply(1, function(x) paste0(substr(x[16], 2, 2), ">", x[17], "::", substr(x[16], 1, 1), "_", substr(x[16], 3, 3))) -> snv_tumor1$sorting


write.table(snv_tumor1, file = "hayward5_filtered.SNVs.isomut", row.names = FALSE, quote = FALSE, sep = "\t")

# plotting: triplet spectra for each sample

for (i in names(table(snv_tumor1$Sample))) {
  pdf(paste0(i, "_spectrum.pdf"), width = 12, height = 7)
  filter(snv_tumor1, Sample == i) %>% with(table(sorting)) %>% barplot(col = rep(rainbow(6), each = 16), las = 2, cex.names = .7, main = i)
  dev.off()
}

# plotting: rainfall plots

seqNames <- seqnames(Hsapiens)[c(1:23)]
seqLengths <-  seqlengths(Hsapiens)[c(1:23)]
cumulativeLengths <- Reduce(sum, seqLengths, accumulate = TRUE)
genomeLength <- sum(seqLengths)

normalizeGenomicPositions <- function(chr, pos) {
  # chr <- paste0("chr", chr)
  pos <- as.numeric(pos)
  whichChrom <- which(seqNames %in% chr)
  earlierChromsLength <- sum(seqLengths[0:(whichChrom - 1)])
  return(as.numeric(earlierChromsLength) + pos)
}

filter(snv_tumor1, Chr %in% seqNames) -> snv_tumor1
snv_tumor1$NormPos <- apply(snv_tumor1, 1, function(x) normalizeGenomicPositions(x[2], x[3]))


for (i in names(table(snv_tumor1$Sample))) {
  pdf(paste0("_", i, "_rainfall.pdf"), width = 12, height = 7)
  filter(snv_tumor1, Sample == i) %>% with(plot(x = NormPos, y = log10(nextmut), cex = .2, main = i))
  dev.off()
}

# get all close mutation pairs
which(snv_tumor1$nextmut > 1 & snv_tumor1$nextmut <= 100 & snv_tumor1$prevmut > 1000) -> ok
snv_tumor1[c(ok, ok+1),] %>% 
  arrange(Chr, Pos) %>%
  mutate(ID = paste0( gsub("_remap_RMdup.bam", "", Sample), "_", Chr, "_", substr(Pos, 1, 4))) %>%
  group_by(ID) %>%
  summarize(Mut1 = gsub(" ", "", muttype[1]),
            Mut2 = gsub(" ", "", muttype[2]),
            Dist = nextmut[1],
            Chr = Chr[1],
            Pos1 = Pos[1],
            Pos2 = Pos[2],
            Sample = Sample[1],
            Trip1 = tripmut[1],
            Trip2 = tripmut[2]) -> nears

write.table(nears, file = "hayward5_near_pairs.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# extracting pair-specific data for postfiltering using Hayward_et_al_melanoma_postfilter.py

pf <- read.delim("hayward5_near_pairs_postfilt.txt", header = TRUE, na.strings = c("NA", "NaN"))

# near pair filters:
# 1. same AF: non-significant Fisher test
# 2. at most 1 supporting reads in either positions in the normal sample
# 3. at least 41 MAPQ for bot positions in the tumor sample
# 4. at least 3 spanning reads

pf <- filter(pf, !is.na(MAPQ2_t), !is.na(MAPQ1_t))

pf$good_one <- TRUE
pf[which(pf$Spanning_read_n < 3), "good_one"] <- FALSE
pf[which(pf$Fisher.p <= 0.05), "good_one"] <- FALSE
pf[which(pf$MAPQ1_t < 41 | pf$MAPQ2_t < 41), "good_one"] <- FALSE
pf[which(pf$Supp1_n > 1 | pf$Supp2_n > 1), "good_one"] <- FALSE

pf$culprit <- "x"
for (i in 1:nrow(pf)) {
  if (pf[i,]$Spanning_read_n < 3) pf[i,]$culprit <- paste(pf[i,]$culprit, "low_spanning", sep = ",")
  if (pf[i,]$Fisher.p <= 0.05) pf[i,]$culprit <- paste(pf[i,]$culprit, "different_AF", sep = ",")
  if (pf[i,]$MAPQ1_t < 41 | pf[i,]$MAPQ2_t < 41) pf[i,]$culprit <- paste(pf[i,]$culprit, "low_mapq", sep = ",")
  if (pf[i,]$Supp1_n > 1 | pf[i,]$Supp2_n > 1) pf[i,]$culprit <- paste(pf[i,]$culprit, "present_in_normal", sep = ",")
}

















getwd()
