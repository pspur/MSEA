### Command-line input: 1. file to process: this is the output from make_msea_data_vcf.py
###                     2. promoter flag: if there is non-exonic data, set to T
###                     3. cohort name
###                     TODO: 1. Implement argparse package for better command-line parsing

library(dplyr)

process_cohort <- function(fin,prom=F,cohort='placehold') {
    
    get.mutloc.and.strand <- function(x,rg) {
        ref.gene <- rg
        chr <- as.character(x[1])
        mloc <- as.numeric(x[2])
        gene <- as.character(x[7])
        gene.match <- ref.gene[ref.gene$V13==gene, ]
        gene.chr.match <- gene.match[gene.match$V3==chr, ]
        gene.chr.match$pstart <- ifelse(gene.chr.match$V4=='+', 
                                        gene.chr.match$V5-1000, 
                                        gene.chr.match$V6)
        gene.chr.match$pend <- ifelse(gene.chr.match$V4=='+',
                                      gene.chr.match$V5,
                                      gene.chr.match$V6+1000)
        match <- filter(gene.chr.match, pstart <=mloc, mloc <= pend)
        match <- filter(match, row_number() == 1) # if multiple matches, arbitrarily use 1st match
        if (nrow(match)==0) {
            rsid <- NA
            mloc <- -1
            strand <- '0'
        } else {
            rsid <- match$V2
            mloc <- abs(mloc - match$pstart)
            strand <- match$V4
        }
        output <- list('rsid'=rsid,'mloc'=mloc,'strand'=strand)
        return(output)
    }

    if (file.exists('../input/refGene.txt.gz')) { 
        ref.gene <- read.table('../input/refGene.txt.gz', sep = '\t', header=F, stringsAsFactors=F)
        ref.gene <- ref.gene[grep('_', ref.gene$V3, invert=T), ]
    } else {
        stop('refGene file not found.')
    }

    variants <- read.table(fin, sep='\t', header=F, stringsAsFactors=F)
    colnames(variants)[1:20] <- c('chrom','pos','snpid','ref','alt','func','gene','exonic_func','aa_change','deleterious',
                                  'onekgenome','prom','asbns','ass','ansi','ans','dns','dnsi','lof','miss')

    # calc frequency field from gt fields, then drop gt fields
    gt.fields <- variants[-1:-20]
    gt.fields <- data.frame(apply(gt.fields,c(1,2),function(u){if (u=='0|0' | u=='0/0') { u <- 0 } 
                                                               else { u <- 1} }))
    variants$freq <- apply(gt.fields,1,sum)
    variants <- cbind(variants[1:20],variants[length(variants)])

    ###################################################################################################
    #process gene variants
    gvariants <- filter(variants, aa_change!='.')

    nt <- c('A','T','C','G')
    ref.snp <- gvariants$ref %in% nt
    alt.snp <- gvariants$alt %in% nt

    # build new data frame from parsing aa_change field
    var_info <- apply(gvariants, 1, function(u) {strsplit(u['aa_change'], split = ',')[[1]]})
    var_info <- lapply(var_info, function(w) {strsplit(w, split = ':')})
    gv <- data.frame(matrix((unlist(var_info)), ncol=5, byrow=T),
                     chrom = rep(gvariants$chrom, sapply(var_info, length)),
                     start = rep(gvariants$pos, sapply(var_info, length)),
                     exonic_func = rep(gvariants$exonic_func, sapply(var_info, length)),
                     deleterious = rep(gvariants$deleterious, sapply(var_info, length)),
                     freq = rep(gvariants$freq, sapply(var_info, length)),
                     snp = rep(ref.snp & alt.snp, sapply(var_info, length)),
                     stringsAsFactors=F)

    ### why comment these out?
    #mutations = mutations[which(mutations$DNA_change!=""), ]  ### remove occasional errors
    #mutations = mutations[!is.na(match(mutations$RefSeq.ID, names(refseq.length))), ]  ### remove refseq genes with no gene length

    gv$snp <- ifelse(gv$snp=='TRUE','SNP','INDEL')

    # extract mutation positions
    mut_pos <- as.numeric(unlist(lapply(gv$X5, function(u){ 
        if(grepl("_", u)){
            u1 <- substr(u, 3, regexec("_", u)[[1]][1]-1  )
            m <- regexec("[0-9]+", u1)
            a <- regmatches(u1,m)
        } else { 
            m <- regexec("[0-9]+", u)
            m <- regmatches(u,m)
        } } ) ))
    gv$mut_pos <- mut_pos

    mutations <- cbind(gv,gvariants[12:20])
    colnames(mutations)[1:5] <- c('Symbol', 'RefSeq.ID', 'Exon.index', 'DNA_change', 'amino_acid_change')

    if (prom) { 
        pvariants <- variants[variants$prom==1, ]

        # drop downstream genes from 'upstream;downstream' rows
        pvariants$gene <- as.character(lapply(pvariants$gene, function(u) { strsplit(u, split='\\\\x3b')[[1]][1]}))

        pvariants$chrom <- paste0('chr', pvariants$chrom)

        # expand rows with multiple genes, drop genotype fields
        a <- strsplit(pvariants$gene, ',')
        pvariants <- data.frame(chrom = rep(pvariants$chrom, sapply(a, length)), 
                                pos = rep(pvariants$pos, sapply(a, length)),
                                snpid = rep(pvariants$snpid, sapply(a, length)),
                                ref = rep(pvariants$ref, sapply(a, length)),
                                alt = rep(pvariants$alt, sapply(a, length)),
                                func = rep(pvariants$func, sapply(a, length)),
                                gene = unlist(a),
                                exonic_func = rep(pvariants$exonic_func, sapply(a, length)),
                                aa_change = rep(pvariants$aa_change, sapply(a, length)),
                                deleterious = rep(pvariants$deleterious, sapply(a, length)),
                                freq = rep(pvariants$freq, sapply(a, length)),
                                onekgenome = rep(pvariants$onekgenome, sapply(a, length)),
                                prom = rep(pvariants$prom, sapply(a, length)),
                                asbns = rep(pvariants$asbns, sapply(a, length)),
                                ass = rep(pvariants$ass, sapply(a, length)),
                                ansi = rep(pvariants$ansi, sapply(a, length)),
                                ans = rep(pvariants$ans, sapply(a, length)),
                                dns = rep(pvariants$dns, sapply(a, length)),
                                dnsi = rep(pvariants$dnsi, sapply(a, length)),
                                lof = rep(pvariants$lof, sapply(a, length)),
                                miss = rep(pvariants$miss, sapply(a, length)),
                                stringsAsFactors=F)

        # remove promoters for genes that don't exist in refGene
        pvariants <- filter(pvariants, gene %in% ref.gene$V13)

        # calculate mut locs referenced to gene start locations and get strand info
        system.time(prominfo <- apply(pvariants, 1, get.mutloc.and.strand, rg=ref.gene))[3]
        prominfo <- matrix(unlist(prominfo), nrow=length(prominfo), byrow=T)
        prominfo <- data.frame(prominfo, stringsAsFactors=F)
        pvariants$RefSeq.ID <- prominfo$X1
        pvariants$mut_pos <- as.integer(prominfo$X2)
        pvariants$strand <- prominfo$X3

        # remove genes from promoters that had mismatched refGene info (mut_loc = -1)
        pvariants <- filter(pvariants, mut_pos!=-1)

        nt <- c('A','T','C','G')
        ref.snp <- pvariants$ref %in% nt
        alt.snp <- pvariants$alt %in% nt
        pvariants$snp <- ref.snp & alt.snp
        pvariants$snp <- ifelse(pvariants$snp=='TRUE','SNP','INDEL')

        pv <- data.frame(Symbol = pvariants$gene,
                         RefSeq.ID = pvariants$RefSeq.ID,
                         Exon.index = rep('-', nrow(pvariants)),
                         DNA_change = rep('-', nrow(pvariants)),
                         amino_acid_change = rep('-', nrow(pvariants)),
                         chrom = pvariants$chrom,
                         start = pvariants$pos,
                         exonic_func = rep('-', nrow(pvariants)),
                         deleterious = rep('-', nrow(pvariants)),
                         freq = pvariants$freq,
                         snp = pvariants$snp,
                         mut_pos = pvariants$mut_pos,
                         stringsAsFactors=F)

        pvariants <- cbind(pv, pvariants[13:21])
        mutations <- rbind(mutations,pvariants)
    }
    write.table(mutations, file=paste0('../input/',cohort,'/vartype_db.plot_ready.txt'), quote=F, sep='\t')
}

args <- commandArgs(trailingOnly = T)
infile <- args[1]
prom <- as.logical(casefold(args[2],upper=T))
#cohort <- strsplit(tail(strsplit(infile,'/')[[1]], n=1),'\\.')[[1]][1]
cohort <- args[3]
process_cohort(infile,prom,cohort)

