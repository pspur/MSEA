from __future__ import print_function
import re

# Returns dict created from provided .gbk flat file.
# Dict structure: protein.ID:[symbol,refseq.ID,aa.len]
def createPIDdict():
  genes = {}
  attributes = []
  with open('./protein.gbk','r') as fin:
    for line in fin:
      line = line.strip()
      if line.startswith('LOCUS'):
        cols = line.split()
        protein_ID = cols[1]
        aa_length = cols[2]
      elif line.startswith('DBSOURCE'):
        refseq_ID = line.split()[-1].split('.')[0]
      elif 'gene=' in line:
        symbol = re.search('"([-\w\.]+)"',line).group(1)
        if 'NP' in protein_ID:
          genes[protein_ID] = [symbol,refseq_ID,aa_length]
  return(genes)

# Takes dict created by createPIDdict.
# Returns nothing, outputs protein domain file.
def parseHitData(genes):
  with open('./hitdata.txt','r') as fin, \
         open('../input/protein.gbk.regions.symbol.txt','w') as fout:
    fout.write('symbol\trefseq.ID\tprotein.ID\taa.length\tdomain.start\t'
               'domain.end\tdomain.source\tdomain.name\tdomain.anno\tdomain.type\n')
    for line in fin:
      line = line.strip()
      if line.startswith('Q#'):
        cols = line.split('\t')
        protein_ID = cols[0].split(' - ')[-1]
        dom_type = cols[1]
        dom_start = cols[3]
        dom_end = cols[4]
        dom_anno = cols[7]
        dom_name = cols[8]
        if protein_ID in genes and (not dom_type.startswith('super')):
          symbol = genes[protein_ID][0]
          refseq_ID = genes[protein_ID][1]
          aa_length = genes[protein_ID][2]
          fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format
            (symbol,refseq_ID,protein_ID,aa_length,dom_start,dom_end,
             'db_xref',dom_name,dom_anno,dom_type))
          
if __name__ == '__main__':
  pids = createPIDdict()
  df = parseHitData(pids)
