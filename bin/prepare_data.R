### The following contains the procedures to reproduce MSEA input data.
### All files used should reside in MSEA/data

library(dplyr)

###################################################################################################
### Procedure to create 'refGene.txt.gz'.

### 1. Download 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz'
### This contains hg19 mRNA info, including exon intervals.

###################################################################################################
### Procedure to create 'protein.gbk.regions.symbol.txt' 

### 1. Download 'ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/protein/protein.gbk.gz'
### This is a flat file of all human refseq proteins, including predicted (XP_) ones

### 2. Run 'makeProteinID_list.py'. This creates 'pids.txt', a file containing all 
### established (not predicted) protein ids, based on the flat file.

### 3. Go to http://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi and upload
### 'pids.txt' then submit the search. This step has a run time of ~2 hours.

### 4. Save the resulting file locally. Use default name: hitdata.txt.

### 5. Run 'makeDomainTable.py'. This creates protein.gbk.regions.symbol.txt 
### from the flat file and hitdata.txt.

### 6. Place protein.gbk.regions.symbol.txt in MSEA/input

###################################################################################################
### Procedure to create 'tfbs.plot_ready.txt', the promoter equivalent of domains

### 1. Download 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/tfbsConsFactors.txt.gz'
### and 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/tfbsConsSites.txt.gz'

match.tfbs.to.promoter <- function(x) {
    chr <- as.character(x[1])
    mstart <- as.numeric(x[2])
    mend <- as.numeric(x[3])
    mstrand <- as.character(x[5])
    chr.match <- prom.to.test[which(chr==prom.to.test$chrom), ]
    chr.strand.match <- chr.match[which(mstrand==chr.match$strand), ]
    chr.strand.match <- chr.strand.match[match(unique(chr.strand.match$Symbol), chr.strand.match$Symbol), ]
    chr.strand.match$pstart <- chr.strand.match$start-chr.strand.match$mut_pos 
    cs.match <- chr.strand.match[which(chr.strand.match$pstart<=mstart & mend<=chr.strand.match$pstart+1000), ]
    if (length(cs.match[,1])==0) {
        symbol = 'n'
        domain.start <- -1
        domain.end <- -2
    } else {
        symbol = cs.match$Symbol
        domain.start <- mstart-cs.match$pstart
        domain.end <- domain.start+(mend-mstart)
    } 
    return(list('symbol'=symbol,'domain.start'=domain.start,'domain.end'=domain.end))
}

tf_sites.raw <- read.table('./data/tfbsConsSites.txt', sep='\t', header=F, stringsAsFactors=F)
tf_factors.raw <- read.table('./data/tfbsConsFactors.txt', sep='\t', header=F, stringsAsFactors=F)
tf_factors <- filter(tf_factors.raw, V3=='human')
tf_sites <- tf_sites.raw[ ,c(-1,-6)]
tf_sites <- filter(tf_sites, V5 %in% tf_factors$V1)    

# assign tfbs to promoter regions, calculate motif starts and ends relative to promoter start
# run time 6-9 hours
system.time(tfbs_gene_info <- apply(tf_sites, 1, match.tfbs.to.promoter))[3]

# assign returned values to respective columns, messy b/c tfbs_gene_info is list of lists of (sometimes) list
tf_sites$symbol <- lapply(tfbs_gene_info, function(u) { u[[1]]})
tf_sites$domain.start <- lapply(tfbs_gene_info, function(u) { u[[2]]})
tf_sites$domain.end <- lapply(tfbs_gene_info, function(u) { u[[3]]})

# remove all sites that didn't match up with a promoter region
tf_sites<- filter(tf_sites, symbol!='n')

# convert new columns from lists of lists to vectors of characters
tf_sites$symbol <- sapply(tf_sites$symbol, paste0, collapse=",")
tf_sites$domain.start <- sapply(tf_sites$domain.start, paste0, collapse=",")
tf_sites$domain.end <- sapply(tf_sites$domain.end, paste0, collapse=",")

# store new columnns as lists of characters to facilitate expanding any with multiple entries
a <- strsplit(tf_sites$symbol, ',')
b <- strsplit(tf_sites$domain.start, ',')
c <- strsplit(tf_sites$domain.end, ',')

# create final tfbs df with all rows that contain multiple gene entries expanded
tfbs <- data.frame(symbol = unlist(a),
                   refseq.ID = unlist(a),
                   protein.ID = NA,
                   length = 1000, 
                   domain.start = as.integer(unlist(b)),
                   domain.end = as.integer(unlist(c)),
                   domain.source = NA,
                   domain.name = rep(tf_sites$V5, sapply(a, length)),
                   domain.anno = NA,
                   domain.type = 'placehold',
                   zscore = rep(tf_sites$V8, sapply(a, length)),
                   stringsAsFactors=F)

write.table(tfbs, file='./input/tfbs.plot_ready.txt', row.names=F, quote=F, sep='\t')

###################################################################################################
### Procedure to create file used as input to msea.clust

### 1. Process annovared (http://annovar.openbioinformatics.org/en/latest/) multi-sample vcf with make_msea_data_vcf.py. This produces a more 
###    reasonably-sized file with a table-ready format for R.

### 2. Use previous step's output as input to process_cohort.R. If cohort includes non-exonic data
###    set 'prom' flag to True so promoter data is processed.