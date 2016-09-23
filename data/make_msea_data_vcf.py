# Creates variant file including variant set info from annovared .vcf for use by MSEA function
# python make_msea_input.py [input file] [output prefix]
# TODO: modify for use with additional predictors

from __future__ import print_function
import sys
import re

def get_freq(samples):
    freq = 0
    discard = ('0/0','0|0')
    for sam in samples:
        s = sam.split(':')[0]
        if (s not in discard):
            freq += 1
    return(freq)

def assign_to_file(infile,outprefix):
    nt = ['A','C','T','G']
    #func_fields = [] # testing
    #snp_func_fields = [] # testing
    with open(infile,'r') as fin, open(outprefix + '_vartype_db.txt','w') as fout:
        #fout.write('prom\tasbns\tass\tansi\tans\tdns\tdnsi\tchr\tstart\tend\tref\talt\tfunc\tgene\tgenedetail\texonicfunc\taachange\tonekg_all_idx\n')
        
        for line in fin:
            if (not line.startswith('#')):
                line = line.strip()
                cols = line.split('\t')
                vartype = [('isPROM',0),('isASBNS',0),('isASS',0),('isANSI',0),
                           ('isANS',0),('isDNS',0),('isDNSI',0),('isLoF',0),('isMISS',0)]
                func_refgene_idx = re.search('Func.refGene=(.*?);',cols[7]).groups()[0]
                gene_refgene_idx = re.search('Gene.refGene=(.*?);',cols[7]).groups()[0]
                snp_func_idx = re.search('ExonicFunc.refGene=(.*?);',cols[7]).groups()[0]
                esp6500_idx = re.search('esp6500siv2_all=(.*?);',cols[7]).groups()[0]
                onekg_all_idx = re.search('1000g2015aug_all=(.*?);',cols[7]).groups()[0]
                gatk_filter = cols[6]
                ref = cols[3]
                alt = cols[4]
                gt_fields = cols[9: ]
                gt_fields = [gt.split(':')[0] for gt in gt_fields]

                # testing
                #if func_refgene_idx not in func_fields:
                #    func_fields.append(func_refgene_idx)
                #if snp_func_idx not in snp_func_fields:
                #    snp_func_fields.append(snp_func_idx)
                
                if ((gatk_filter == 'PASS') and
                   ( './.' not in gt_fields) and
                   ( 'upstream' in func_refgene_idx)):
                    vartype[0] = ('isPROM',1) # promoter
                    vt = '\t'.join([str(tuple[1]) for tuple in vartype])
                    outcols = '\t'.join(cols[0:5])
                    fout.write('{0}\t{1}\t{2}\t{3}\t.\t.\t{4}\t{5}\t{6}\n'
                               .format(outcols,func_refgene_idx,gene_refgene_idx,
                               snp_func_idx,onekg_all_idx,vt,'\t'.join(gt_fields)))
                
                #elif ((gatk_filter == 'PASS') and
                #     ( './.' not in gt_fields) and
                #     ( 'splicing' in func_refgene_idx) and
                #     ( 'ncRNA_splicing' not in func_refgene_idx)): 
                #    vartype[7] = ('isLoF',1) # splicing
                #    vt = '\t'.join([str(tuple[1]) for tuple in vartype])
                #    outcols = '\t'.join(cols[0:5])
                #    fout.write('{0}\t{1}\t{2}\t{3}\t.\t.\t{4}\t{5}\t{6}\n'
                #               .format(outcols,func_refgene_idx,gene_refgene_idx,
                #               snp_func_idx,onekg_all_idx,vt,'\t'.join(gt_fields)))

                elif ((gatk_filter == 'PASS') and
                     ( './.' not in gt_fields)):
                   #( onekg_all_idx not in '.' and float(onekg_all_idx) <= af_filter) and
                   #( esp6500_idx not in '.' and float(esp6500_idx) <= af_filter)):
               
                    result = []
                    aa_change_idx = re.search('AAChange.refGene=(.*?);',cols[7]).groups()[0].split(',')
                    for entry in aa_change_idx:
                        if (len(entry.split(':')) == 5):
                            tmp = '{0}\t{1}\t{2}\t{3}\t{4}'.format('\t'.join(cols[0:5]),func_refgene_idx,gene_refgene_idx,snp_func_idx,entry) 
                            result.append(tmp)

                    #freq = get_freq(gt_fields)
                    sift_idx = re.search('SIFT_pred=(.*?);',cols[7]).groups()[0]
                    polyphen2_idx = re.search('Polyphen2_HDIV_pred=(.*?);',cols[7]).groups()[0]
                    lrt_idx = re.search('LRT_pred=(.*?);',cols[7]).groups()[0]
                    fathmm_idx = re.search('FATHMM_pred=(.*?);',cols[7]).groups()[0]
                    mut_taster_idx = re.search('MutationTaster_pred=(.*?);',cols[7]).groups()[0]
                    mut_assessor_idx = re.search('MutationAssessor_pred=(.*?);',cols[7]).groups()[0]
                    metaSVM_idx = re.search('MetaSVM_pred=(.*?);',cols[7]).groups()[0]
                    metaLR_idx = re.search('MetaLR_pred=(.*?);',cols[7]).groups()[0]
                
                    is_deleterious = 'No'
                    if (sift_idx == 'D' or
                        polyphen2_idx in 'D|B' or
                        mut_taster_idx in 'A|D' or
                        fathmm_idx == 'D' or
                        metaSVM_idx == 'D' or
                        metaLR_idx == 'D' or
                        mut_assessor_idx in 'H|M' or
                        lrt_idx == 'D'):
                
                        is_deleterious = 'Yes'

                    if (sift_idx == '.' and
                        polyphen2_idx == '.' and
                        mut_taster_idx == '.' and
                        fathmm_idx == '.' and
                        metaSVM_idx == '.' and
                        metaLR_idx == '.' and
                        mut_assessor_idx == '.' and
                        lrt_idx == '.'): 
 
                        is_deleterious = 'Unknown'
                
                    for r in result:
                        if (snp_func_idx == 'synonymous_SNV' or is_deleterious == 'No'):
                            vartype[1] = ('isASBNS',1) # all silent plus benign nonsilent SNV

                            if (snp_func_idx == 'synonymous_SNV'):
                                vartype[2] = ('isASS',1) # all silent SNV

                        if ('nonsynonymous_SNV' in snp_func_idx or
                            'stop' in snp_func_idx or
                            'frameshift' in snp_func_idx):
                            vartype[3] = ('isANSI',1) # all non-silent SNV and Indel

                            if (snp_func_idx.startswith('non')):
                                vartype[8] = ('isMISS',1) # missense

                            if ((ref in nt) and (alt in nt)):
                                vartype[4] = ('isANS',1) # all non-silent SNV

                                if (is_deleterious == 'Yes'):
                                   vartype[5] = ('isDNS',1) # deleterious non-silent SNV

                            if (is_deleterious == 'Yes'):
                                vartype[6] = ('isDNSI',1) # deleterious non-silent SNV and Indel

                        if (snp_func_idx.startswith('frame') or
                            'stopgain' in snp_func_idx):
                            vartype[7] = ('isLoF',1) # loss of function

                        vt = '\t'.join([str(tuple[1]) for tuple in vartype])
                        fout.write('{0}\t{1}\t{2}\t{3}\t{4}\n'
                                   .format(r,is_deleterious,onekg_all_idx,vt,'\t'.join(gt_fields)))    
    
    # testing
    #for i in func_fields:
    #    print(i)

    #print()
    #for i in snp_func_fields:
    #    print(i)

if __name__ == '__main__':
    #af = 0.1
    filein = sys.argv[1]
    outfile_prefix = sys.argv[2]
    assign_to_file(filein,outfile_prefix)
