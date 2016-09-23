### Command-line input: 1. file to process: this is the output from process_cohort.R including path.
###                     2. expand flag: set to True if you want p-values calculated with variant pop freq considered.
###                        As it currently stands, expand=T will artificially deflate p-values 
###                     TODO: 1. Require more command line input to create output filename rather than deriving
###                              from input filename and path
###                           2. Implement argparse package for better command-line parsing
###                           3. Bundle output code into function, with flag to indicate desired output type

library(doParallel)
library(plyr)
library(grid)
library(dplyr)
library(jsonlite)

msea_clust <- function(fin,expand,output_str) {
    
    get_ref_gene_length <- function(x) {
        cds.start <- as.numeric(x[7])
        cds.end <- as.numeric(x[8])
        exon.starts <- as.numeric(unlist(strsplit(x[10], ',')))
        exon.ends <- as.numeric(unlist(strsplit(x[11], ',')))
        
        # get the index of first coding exon
        for (i in 1:length(exon.starts)) {
            if (cds.start >= exon.starts[i] & cds.start <= exon.ends[i]) {
                exon.starts[i] <- cds.start
                start.exon <- i
                break
            }
        }
    
        # get the index of the last coding exon
        for (i in 1:length(exon.ends)){
            if (cds.end >= exon.starts[i] & cds.end <= exon.ends[i]) {
                exon.ends[i] <- cds.end
                end.exon <- i
                break
            }
        }
        
        exon.starts <- exon.starts[start.exon:end.exon]
        exon.ends <- exon.ends[start.exon:end.exon]
        
        len <- sum(exon.ends - exon.starts)
        return(len/3) # convert to a.a. length
    }
    
    mutations <- read.table(fin, sep='\t', header=T, stringsAsFactors=F)
    ref_gene <- read.table('../input/refGene.txt.gz', sep = '\t', header=F, stringsAsFactors=F)
    domain <- read.delim('../input/protein.gbk.regions.symbol.txt', header=T, stringsAsFactors=F, sep='\t')
    tfbs <- read.table('../input/tfbs.plot_ready.txt', header=T, stringsAsFactors=F, sep='\t')
    
    refseq_length <- apply(ref_gene, 1, get_ref_gene_length)
    names(refseq_length) <- ref_gene$V2
    
    # exclude genes with < 4 mutations and genes without entries in refgene
    filter_valid <- function(var_set) {
        nMut_per_gene <- tapply(var_set$amino_acid_change, var_set$RefSeq.ID, length)
        genes_of_length <- names(which(nMut_per_gene >= 4))
        valid_genes <- genes_of_length[genes_of_length %in% names(refseq_length)]
        valid_var_set <- filter(var_set, RefSeq.ID %in% valid_genes)
        return(valid_var_set)
    }
    
    # test data
    #mutations <- mutations[mutations$Symbol=='SAMD11' | 
    #                       mutations$Symbol=='NOC2L' |
    #                       mutations$Symbol=='KLHL17' |
    #                       mutations$Symbol=='PLEKHN1', ]
    
    asbns <- filter_valid(mutations[mutations$asbns==1, ])
    ass <- filter_valid(mutations[mutations$ass==1, ])
    ansi <- filter_valid(mutations[mutations$ansi==1, ])
    ans <- filter_valid(mutations[mutations$ans==1, ])
    dns <- filter_valid(mutations[mutations$dns==1, ])
    dnsi <- filter_valid(mutations[mutations$dnsi==1, ])
    lof <- filter_valid(mutations[mutations$lof==1, ])
    miss <- filter_valid(mutations[mutations$miss==1, ])
    prom <- filter_valid(mutations[mutations$prom==1, ])
    mutations <- unique(rbind(lof,miss,asbns,ass,ansi,ans,dns,dnsi,prom))
        
    nMut_per_gene <- tapply(mutations$amino_acid_change, mutations$RefSeq.ID, length)
    genes_to_test <- names(which(nMut_per_gene >= 4))
    
    # print some numbers
    #print(paste('# input genes: ', length(nMut.per.gene), 
    #            '; eligible genes to test: ', length(genes_to_test), sep=''))
    
    if (length(genes_to_test) == 0){
        stopCluster(cl)
        stop('No eligible genes to test.')
    }
    
    #sink(file=paste0(output_str,'_lookahead.txt'))
    #for (g in 1:length(genes_to_test)) {
    #    mut <- mutations[mutations$RefSeq.ID==genes_to_test[g], c(1,2,22)]
    #    cat(paste(mut$Symbol[1],mut$RefSeq.ID[1],paste(unique(mut$vtype),collapse='\t'), sep='\t'), '\n')
    #}
    #sink()
    
    q <- foreach (i=1:length(genes_to_test)) %dopar% {
        curr_gene <- mutations[mutations$RefSeq.ID==genes_to_test[i], ]
        refseq.ID <- genes_to_test[i]
        
        #expand each row to match number of mutations (freq)
        if (expand) {
            unexpanded_gene <- curr_gene
            curr_gene <- curr_gene[rep(seq_len(nrow(curr_gene)), curr_gene$freq), ]
            row.names(curr_gene) <- NULL
        } else {
            unexpanded_gene <- curr_gene
        }
        
        result <- list()
        # parallelize this loop as well somehow?
        for (kk in 13:21) {
            curr_vset <- curr_gene[curr_gene[, kk]==1,]
            unexp_vset <- unexpanded_gene[unexpanded_gene[, kk]==1,]
            if (nrow(curr_vset) > 0) {
                ptm <- proc.time()
                symbol <- curr_vset$Symbol[1]
                if (curr_vset$prom[1]==1) {
                    gene_length <- 1000
                    gene_domains <- tfbs[tfbs$refseq.ID==symbol, ]
                } else {
                    gene_length <- unname(refseq_length[refseq.ID])
                    gene_domains <- domain[domain$refseq.ID==refseq.ID, ]
                }     
                mut_pos <- curr_vset$mut_pos
                snp <- curr_vset$snp
                freq <- curr_vset$freq
                vset <- names(curr_vset[kk])
                    
                ### es.random
                es.random <- vector(length = gene_length*10, mode = 'numeric')
                for(ii in 1:(gene_length*10)){
                    Nm <- nrow(curr_vset)
                    
                    # randomly select the same number of variants as in input data
                    mut_pos.pai <- sample(1:gene_length, Nm, replace=T)
                    inc <- 1/length(mut_pos.pai) # same as Nm
                    nMut.per.location <- table(mut_pos.pai)
                    dec <- 1/(gene_length-length(nMut.per.location))
                    inc.1 <- rep(0, gene_length)
                    inc.1[as.numeric(names(nMut.per.location))] <- inc*nMut.per.location
                    dec.1 <- rep(-dec, gene_length)
                    dec.1[as.numeric(names(nMut.per.location))] <- 0
                    ss <- inc.1 + dec.1
                    es.cum <- unname(cumsum(ss))
                    es.pai <- (max(es.cum)-min(es.cum))
                    es.random[ii] <- es.pai
                }
                
                ### es true
                inc <- 1/length(mut_pos)
                nMut.per.location <- table(mut_pos)
                dec <- 1/(gene_length-length(nMut.per.location))
                inc.1 <- rep(0, gene_length)
                inc.1[as.numeric(names(nMut.per.location))] <- inc*nMut.per.location
                dec.1 <- rep(-dec, gene_length)
                dec.1[as.numeric(names(nMut.per.location))] <- 0
                ss <- inc.1 + dec.1
                es.cum <- cumsum(ss)
                true.es.cum <- es.cum
                es.true <- (max(es.cum)-min(es.cum))
                pvalue <- sum(es.random>=es.true)/(length(es.random)+1)
                nes <- (es.true-mean(es.random))/sd(es.random)
                p1 <- sum(es.random>=es.true)/(length(es.random)+1)
                runtime <- round((proc.time() - ptm)[[3]], 4)
                print(paste('Run time of ', symbol, ' (', genes_to_test[i], '), ', 
                            vset, ': ', runtime, 's', sep = ''))
                x_coord <- 1:gene_length
                unexp_mutpos <- unexp_vset$mut_pos
                
                # structure output so that it will convert to the desired json
                result[[kk]] <- list(gene_name=symbol, 
                                     mutation_accumulation_score=c(rbind(x_coord,es.cum)), 
                                     refgene_id=refseq.ID, 
                                     vset=vset, 
                                     nes=nes,
                                     chrom=curr_vset$chrom[1], 
                                     length=gene_length, 
                                     pvalue=pvalue, 
                                     mutations=list(start=c(rbind(unexp_mutpos,unexp_vset$start)),
                                                    exonic_functions=c(rbind(unexp_mutpos,unexp_vset$exonic_func)),
                                                    frequencies=c(rbind(unexp_mutpos,unexp_vset$freq)),
                                                    mutation_types=c(rbind(unexp_mutpos,unexp_vset$snp)),
                                                    deleterious=c(rbind(unexp_mutpos,unexp_vset$deleterious))))                                     
                
            }
        }
        
        # drop any NULL elements that result from a gene not having every type of variant set
        result <- result[!sapply(result, is.null)]
    }
    
    # flatten results down to one list of lists
    q <- unlist(q, recursive=F)
    
    # extract all p-values by variant set, perform multiple testing correction, assign corrected values back in
    vsets <- vector()
    for (j in 1:length(q)) {
        vsets <- c(vsets, q[[j]]$vset)
    }
    vsets <- unique(vsets)
    for (i in 1:length(vsets)) {
        pvalues <- NULL
        vset_indices <- which(sapply(q, '[[', 4) == vsets[i])
        for (ii in 1:length(vset_indices)) {
            pvalues[[ii]] <- q[[vset_indices[ii]]]$pvalue
        }
        pvalues <- p.adjust(pvalues,method='fdr')
        for (ii in 1:length(vset_indices)) {
            q[[vset_indices[ii]]]$pvalue <- pvalues[[ii]]
        }
    }
    
    # output results to R workspace file
    save(q, file=paste0(output_str,"_Data4plot.RData"))
    
    # output results to one file as serialized json entries
    sink(file=paste0(output_str,'.json'))
    for (g in 1:length(q)) {
        cat(toJSON(q[[g]],pretty=TRUE),'\nnext_entry\n')
    }    
    sink()     
}

args <- commandArgs(trailingOnly = T)
infile <- args[1]
expand <- as.logical(casefold(args[2],upper=T))

cl <- makeCluster(8,outfile='')
registerDoParallel(cl)
DATE <- Sys.Date()
a <- strsplit(infile,'input/')[[1]][2]
b <- strsplit(a,'\\.')[[1]][1]
output_stem <- paste0(gsub('/','_',b), '_', DATE)
output_directory <- '../output'
if (expand){
    expand_str <- 'expand'
}else{
    expand_str <- 'noexpand'
}
output_str <- paste0(output_directory, '/', output_stem, '_', expand_str)

msea_clust(infile,expand,output_str)

DATE <- Sys.Date()
### WRITE SESSION INFO TO FILE #####################################################################
sink(file=paste0(output_str, "_Session_Info.txt"));

## write memory usage to file
cat('### Memory #########################################################################################\n');
print(gc());

## write process time to file
cat('\n### Time ###########################################################################################\n');
print(proc.time());

## write list of objects
cat('\n### ls() ###########################################################################################\n');
print(ls(pos = 1));

## write sessionInfo to file
cat('\n### Session Info ###################################################################################\n');
print(sessionInfo());

## close the file
sink();

stopCluster(cl)
