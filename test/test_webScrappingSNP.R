### Packages used

library(magrittr)
library(rvest)
library(stringr)


### Data loading

setwd("..")
source("code/webScrappingSNP.R")
asthma <- read.table("data/asthma.txt", header=TRUE)

all_snps <- names(asthma)[7:ncol(asthma)]
some_snps <- all_snps[c(1,    # NA gene 
                        26,   # NA 
                        27:30 # full data
                        )]

### Test functions

## load_ncbi_datanode
load_ncbi_datanode(some_snps[1])
load_ncbi_datanode(some_snps[2])
load_ncbi_datanode(some_snps[3])
datanodes <- lapply(some_snps, load_ncbi_datanode)
datanodes


## extractRow
extractRow(datanodes[[1]], "Gene")
extractRow(datanodes[[2]], "Gene")
extractRow(datanodes[[3]], "Gene")
lapply(datanodes, extractRow, row_name = "Gene")
lapply(datanodes, extractRow, row_name = "Chromosome")


## extractGen
extractGen(datanodes[[1]])
extractGen(datanodes[[2]])
extractGen(datanodes[[3]])
lapply(datanodes, extractGen)


## extractCh
extractCh(datanodes[[1]])
extractCh(datanodes[[2]])
lapply(datanodes, extractCh)


## Extract a dataframe
lapply(some_snps, function(snp_name){
  datanode <- load_ncbi_datanode(snp_name)
  chromosoma = extractCh(datanode)
  data.frame(Gene = extractGen(datanode),
             Ch = chromosoma$ch,
             Pos = chromosoma$pos)
}) %>% 
  do.call(what = rbind)
