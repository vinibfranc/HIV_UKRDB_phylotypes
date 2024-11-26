import sys
from Bio import SeqIO

FastaFile = open(sys.argv[1], 'r')
FastaDroppedFile = open(sys.argv[2], 'w')
drop_cutoff = float(sys.argv[3])

for seqs in SeqIO.parse(FastaFile, 'fasta'):
    name = seqs.id
    seq = seqs.seq
    #print(seq)
    seqLen = len(seqs)
    #print(seqLen)
    gap_count = 0
    for z in range(seqLen):
        if seq[z]=='-': #before was ?
            gap_count += 1
    #print(seqLen-gap_count)
    if (seqLen-gap_count) < drop_cutoff:
        #print(seqLen-gap_count)
        print(' %s was removed.' % name)
    else:
        SeqIO.write(seqs, FastaDroppedFile, 'fasta')

FastaFile.close()
FastaDroppedFile.close()
