from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def align_sequences(seq1, seq2):
    """
    Perform pairwise sequence alignment using Biopython's pairwise2 module.
    """
    from Bio import pairwise2
    alignments = pairwise2.align.globalxx(seq1, seq2)
    return alignments[0]

def trim_sequence(sequence, start, end):
    """
    Trim the sequence to the specified start and end positions.
    """
    return sequence[start:end]

def count_unique_sequences(sequences):
    """
    Count the frequency of unique sequences in a list of sequences.
    """
    sequence_counts = {}
    for seq in sequences:
        if seq in sequence_counts:
            sequence_counts[seq] += 1
        else:
            sequence_counts[seq] = 1
    return sequence_counts

#Read reference sequence

#interest_seq = "AGGGGTATCTCTCGATAAAAGAGAGGCTGAAGCTTGGTGCNNNNNNNNNGGACAACCTTGTTGGTAATAG"
reference_seq = "AGCTTCAGCCTCTCTTTTATCGAGAGATACC"
#lengh 70 bp

# Read fastq file
fastq_records = SeqIO.parse("8FoAgo2_S8_L001_R1_001.fastq", "fastq")


# Quality control, alignment, trimming, and counting
aligned_sequences = []
counter=0
for record in fastq_records:
    counter+=1
    alignment = align_sequences(reference_seq, record.seq)
    #print(alignment)
    trimmed_seq = trim_sequence(record.seq, len(reference_seq)+35, len(reference_seq)+71)
    aligned_sequences.append(str(trimmed_seq))

    if counter %1000 ==0:
        print (counter, "trim seq is ",trimmed_seq)

# Count unique sequences
sequence_counts = count_unique_sequences(aligned_sequences)
#print(sequence_counts)


# Convert counts to DataFrame
df = pd.DataFrame.from_dict(sequence_counts, orient='index', columns=['Frequency'])
df.index.name = 'Sequence'
df.reset_index(inplace=True)

# Sort DataFrame by frequency
df.sort_values(by='Frequency', ascending=False, inplace=True)


df.to_csv('results8FoAgo2.csv', index=False)


#Plot
top200=df.head(200)
translationlist=[]

for i in top200["Sequence"]:
    seq=Seq(i)
    translation=seq.translate()
    translationlist.append(str(translation))
top200["translation"]=translationlist
print(top200)

top200.to_csv('top200_7FoAgo1.csv', index=False)
yo=top200.head(30)
fig, ax = plt.subplots()
ax.barh(yo["translation"], yo["Frequency"])
