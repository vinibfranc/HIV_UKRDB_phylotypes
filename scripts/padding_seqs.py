from Bio import SeqIO
import sys

file_in = sys.argv[1]
file_out = sys.argv[2]

sequences = [s for s in SeqIO.parse(file_in, 'fasta')]
max_len = max([len(s.seq) for s in sequences])
GAPs = "-"
for seq in sequences:
    padding = GAPs*(max_len - len(seq.seq)) # creating the padding string
    seq.seq += padding

SeqIO.write(sequences, file_out, 'fasta')
