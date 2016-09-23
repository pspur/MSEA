# Creates file containing all NP_ protein ids, each id on a new line.
# To be used as input to ncbi conserved domain batch search.
# Everything hard-coded for now.

with open('./protein.gbk','r') as fin, open('./pids.txt','w') as fout:
  for line in fin:
    line = line.strip()
    if line.startswith('LOCUS'):
      cols = line.split()
      protein_ID = cols[1]
      if 'NP' in protein_ID:
        fout.write('{0}\n'.format(protein_ID))

