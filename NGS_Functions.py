# Functions for handling/analyzing NGS data

import matplotlib.pyplot as plt

# Convert the ASCII-encoded quality using "Phred+33"
def Q_to_Phred33(Q):
    # chr converts the ASCII character to its corresponding integer
    return chr(Q+33)

# Convert the quality integer to its ASCII-encoded character
def Phred33_to_Q(n):
    # ord reverse-converts integer to ASCII character
    return ord(n)-33

# Read the genome sequence from a .fa file
def FASTA_read(file):
    genome = ''
    with open(file, 'r') as f:
        for line in f:
            # The first line of FASTA files begin with ">"
            # The second line is the start of the sequence
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

# Read the sequences and qualities of a FASTQ file
def FASTQ_read(file):
    sequences = []
    qualities = []
    with open(file) as f:
        while True:
            # Skip name line
            f.readline().rstrip()
            # Read sequence line
            seq = f.readline().rstrip()
            # Skip placeholder line
            f.readline().rstrip()
            # Read qualities line
            qual = f.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

# Create histogram of qualities
def hist(qualities):
    h = [0] * 50
    for qual in qualities:
        for phred in qual:
            q = Phred33_to_Q(phred)
            h[q] += 1
    return h

def find_GC_bypos(reads):
    gc_bases = [0] * 100
    total_bases = [0] * 100
    for read in reads:
        for nuc in range(len(read)):
            if read[nuc] == 'C' or read[nuc] == 'G':
                gc_bases[nuc] += 1
                total_bases[nuc] += 1
    for nuc in range(len(gc_bases)):
        if total_bases[nuc] > 0:
            gc_bases[nuc] /= float(total_bases[nuc])
    return gc_bases

# Naive exact matching algorithm for matching artificial reads
def naive(pattern, genome):
    # Indices where pattern matches against genome
    occurrences = []
    for nuc in range(len(genome) - len(pattern) + 1):
        match = True
        for j in range(len(pattern)):
            if not genome[nuc+j] == pattern[j]:
                # If there is a mismatch
                match = False
                break
        if match:
            occurrences.append(nuc)
    # Offset of every match of pattern within genome
    return occurrences