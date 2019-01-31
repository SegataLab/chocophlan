from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ifn', type=str)
parser.add_argument('ofn', type=str)

args = parser.parse_args()
readsize=150
ifn = args.ifn
ofn = args.ofn
with open(ofn, "w") as output_handle:
    for record in SeqIO.parse(ifn, "fasta"):
        start = 0
        tlen = len(record.seq)//readsize
        for i in range(1,tlen+1):
            subseq = record.seq[start:readsize*i]
            new_record = SeqIO.SeqRecord(seq=subseq, id='{}.{}/{}'.format(record.id,i, tlen), description  = '')
            start=readsize*i
            SeqIO.write(new_record, output_handle, "fasta")