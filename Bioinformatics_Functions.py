# .py file containing various Bioinformatics functions

import random
import itertools
from itertools import product
from collections import defaultdict
import matplotlib.pyplot as plt
import math

# Count the number of nucleotides
def nucleotide_count(seq):
    nuc_counts = {'A':seq.count('A'), 'C':seq.count('C'), 'G':seq.count('G'), 'T':seq.count('T')}
    return nuc_counts

# Obtain the reverse complement of a DNA sequence
def reverse_complement(seq):
    # Assign complementary letters
    comp_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    # Replace each nucleotide with its complementary letter
    comp_seq = ''.join(comp_dict[n] for n in seq)
    # Reverse the entire string
    rev_comp_seq = comp_seq[::-1]
    return rev_comp_seq

# Randomly generate a DNA sequence of length n
def random_sequence(n):
    nucleotides = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(nucleotides) for _ in range(n))

# Generate a list of random DNA strings to test DNAMatrix class functions
def random_dna_list(num_sequences, sequence_length):
    sequences = []
    for _ in range(num_sequences):
        sequence = ''.join(random.choice('ATCG') for _ in range(sequence_length))
        sequences.append(sequence)
    return sequences

# Count the occurrences of a given pattern in a sequence
def pattern_count(seq, pattern):
    # Store 1 if pattern is encountered, and sum total
    return sum([1 for nuc in range(len(seq) - len(pattern)) if seq[nuc:nuc + len(pattern)] == pattern])

# Find locations of a k-mer pattern within a DNA sequence
def find_pattern_indices(seq, pattern):
    occurrences = []
    k = len(pattern)
    # For whole sequence - k
    for i in range(len(seq) - k + 1):
        # Checking for our pattern using a scanning window of size k 
        if seq[i:i + k] == pattern:
            # If our pattern is found, append the index to occurences list
            occurrences.append(i)
    occurrences_nocommas = ' '.join(map(str, occurrences))
    return occurrences_nocommas

# More efficient HammingDistance function using list comprehension
def hamming_distance_V2(p, q):
    if len(p) != len(q):
        raise ValueError("Strings must have the same length")
    return sum([1 for nuc in range(len(p)) if p[nuc] != q[nuc]])

# Count the number of approximate pattern occurrences with d or fewer mismatches
def approx_pattern_count(seq, pattern, d):
    positions = []
    for nuc in range(len(seq) - len(pattern) + 1):
        # Sliding window
        if hamming_distance_V2(pattern, seq[nuc:nuc + len(pattern)]) <= int(d):
            positions.append(nuc)
    return len(positions)

# Find all approximate occurrences of a pattern in a string (d or fewer mismatches from a given pattern)
def approx_pattern_occurrences(pattern, seq, d):
    pattern_len = len(pattern)
    seq_len = len(seq)
    positions = []
    # Sliding window
    for i in range(len(seq) - len(pattern) + 1):
        if hamming_distance_V2(pattern, seq[i:i + len(pattern)]) <= int(d):
            positions.append(i)
    return ' '.join(map(str, positions))

# Create a frequency map of k-mer pattern occurrences as a dictionary
def freq_map(seq, k):
    freqMap = {}
    # For the whole sequence - k
    for i in range(len(seq) - int(k) + 1):
        # Scanning window of letters (of length k)
        pattern = seq[i:i + int(k)]
        # For each new k-length segment, it sets its value to 1 in dictionary
        if pattern not in freqMap:
            freqMap[pattern] = 1
        # If k-length segment is already in dictionary, it increments its value (initially 1)
        freqMap[pattern] += 1
    return freqMap

# Retrieve the most frequent patterns in frequency table
def most_freq_patterns(seq, k):
    freqPatterns = []
    # Get freqMap dictionary
    freqMap = freq_map(seq, int(k))
    # For each pattern (key), count (value) item in dictionary
    for pattern, count in freqMap.items():
        # Find patterns with max count values
        if count == max(freqMap.values()):
            # Append freqPatterns list with those patterns
            freqPatterns.append(pattern)
    return freqPatterns

# Generate all theoretically possible neighbors
def neighbors(pattern, d):
    nucleotides = 'ACGT'
    neighborhood = set()
    # Positions in the pattern that could be mismatched: (0,0), (0,1), etc. if d=2
    for mismatched_indices in product(range(len(pattern)), repeat=d):
    # Result if pattern length=3 and d=2:
    # (0,0), (0,1), (0,2)
    # (1,0), (1,1), (1,2)
    # (2,0), (2,1), (2,2)

        # Nucleotide substitutions that will be applied: (A,A), (A,C), etc. if d=2
        for mismatches in product(nucleotides, repeat=d):
        # Result if pattern length=3 and d=2:
        # (A,A), (A,C), (A,G), (A,T)
        # (C,A), (C,C), ...
        # (G,A), ...
        # (T,A), ...

            neighbor = list(pattern)
            # Combine elements from mismatched_indices and mismatches as tuples: (0, A), (0, C), etc.
            for index, mismatch in zip(mismatched_indices, mismatches):
                # Create mismatch substitutions
                neighbor[index] = mismatch
            # All theoretically possible mismatches
            neighborhood.add(''.join(neighbor))
    return list(neighborhood)
    
# Generates all possible neighbors of a pattern with d or fewer mismatches
def neighbors_V2(pattern, d):
    nucleotides = 'ACGT'
    
    # Case: d = 0
    if d == 0:
        return {pattern}
    # Case: pattern is 1 letter
    if len(pattern) == 1:
        return set(nucleotides)

    # Case: d > 0 and pattern > 1 letter
    neighborhood = set()
    # Starting with possible mismatches at the first letter, keeping the suffix of the pattern the same
    # Then, use recursion to iterate over all possible mismatches at each remaining letter index
    suffix_neighbors = neighbors_V2(pattern[1:], d)

    for neighbor in suffix_neighbors:
        # If the Hamming distace of the suffix is < d, the neighbor differs by at most d mismatches
        if hamming_distance_V2(pattern[1:], neighbor) < int(d):
            for nucleotide in nucleotides:
                # Concatenate each of the four nucleotides to the suffix of the pattern
                neighborhood.add(nucleotide + neighbor)
        else:
            neighborhood.add(pattern[0] + neighbor)
    return list(neighborhood)

# Find the most frequently occurring mismatched patterns
    # This function returns the most frequent neighbor(s) AND a freqMap of all neighbors
def most_freq_mismatched_patterns(seq, k, d):
    neighbor_freqMap = {}
    most_freq_neighbors = []
    # Find existing neighbors
    for nuc in range(len(seq) - int(k) + 1):
        # Sliding window of pattern length
        pattern = seq[nuc:nuc + int(k)]
        # Function call to generate neighborhood of all possible neighbors of a pattern with up to d mismatches
        neighborhood = neighbors(pattern, int(d))
        for neighbor in neighborhood:
            # Set all patterns = 0 in freqMap
            if neighbor not in neighbor_freqMap:
                neighbor_freqMap[neighbor] = 0
            # If pattern is found, increment in freqMap
            else:
                neighbor_freqMap[neighbor] += 1
                
    # Find most occurring k-mer with d or fewer mismatches and add it to patterns
    most_freq_neighbors = [neighbor for neighbor, freq in neighbor_freqMap.items() if freq == max(neighbor_freqMap.values())]
    # Most frequent neighbor(s), all neighbors
    return ' '.join(most_freq_neighbors), neighbor_freqMap

# Finding the most frequent mismatched k-mers and obtaining their reverse compliments
    # This function returns the most frequent mismatched k-mer and its reverse compliment
def most_freq_mismatched_patterns_revcomp(seq, k, d):
    patterns = []
    freqMap = {}
    for nuc in range(len(seq) - int(k) + 1):
        pattern = seq[nuc:nuc + int(k)]
        # Generate the pattern's reverse compliment
        pattern_rc = reverse_complement(pattern)
        neighborhood = neighbors(pattern, int(d))
        # Generate a neighborhood of reverse compliments
        neighborhood_rc = [reverse_complement(n) for n in neighborhood]
        
        for neighbor in neighborhood + neighborhood_rc:
            if neighbor not in freqMap:
                freqMap[neighbor] = 1
            else:
                freqMap[neighbor] += 1

    max_freq = max(freqMap.values())
    patterns = [pattern for pattern, freq in freqMap.items() if freq == max_freq]
    return ' '.join(patterns)

# The difference between the total number of occurrences of G and C in the first i nucleotides of Genome
def skew(genome):
    # Initialize cumulative sum
    skew_vals = [0]
    for nuc in genome:
        if nuc == 'C':
            # if C, subtract one from previous value
            skew_vals.append(skew_vals[-1] - 1)
        # if G, add one to previous value
        elif nuc == 'G':
            skew_vals.append(skew_vals[-1] + 1)
        # if A or T, simply retain the previous value
        else:
            skew_vals.append(skew_vals[-1])
    return skew_vals

# Find the minimum skew position(s), where the reverse half-strand ends and the forward half-strand begins
def min_skew_positions(genome):
    skew_vals = skew(genome)
    min_skew_val = min(skew_vals)
    # Obtain all minimum skew value indices
        # enumerate() returns an index, value pair (i, skew) in the iterable (skew_val list)
    positions = [i for i, skew in enumerate(skew_vals) if skew == min_skew_val]
    return positions

# Create skew plot from skew()
def plot_skew(genome):
    # y-axis
    skew_vals = skew(genome)
    # x-axis
    positions = range(len(genome) + 1)    
    plt.plot(positions, skew_vals, color='black')
    plt.xlabel('Position (i)')
    plt.ylabel('Skew')
    plt.title('Skew Plot')
    plt.grid(True)
    plt.show()
    
# Find (k,d)-motifs in a given collection of DNA sequences
    # The motifs/neighbors must be present in all sequences/strings
    # This function uses a list of strings (DNA sequences)
def motif_enumeration(Dna: list, k: int, d: int) -> str:
    patterns = set()
    # For each DNA sequence - k substring
    for nuc in range(len(Dna[0]) - int(k) + 1):
        # Sliding window of length k to find all patterns in each sequence
        pattern = Dna[0][nuc:nuc + int(k)]
        # Get all possible neighbors of each pattern in sliding window
        pattern_neighbors = neighbors_V2(pattern, d)
        motifs_found = set()
        # For each neighbor in the set of all possible neighbors for all patterns across all sequences
        for neighbor in pattern_neighbors:
            is_motif = True
            # Check each DNA sequence for the current neighbor
            for seq in Dna:
                motif_found = False
                # Check current DNA sequence substring
                for nuc in range(len(seq) - int(k) + 1):
                    # If the number of mismatches is d or fewer in a sliding window of patterns
                    if hamming_distance_V2(neighbor, seq[nuc:nuc + int(k)]) <= int(d):
                        # Motif found and move to the next substring in the current sequence
                        motif_found = True
                        # Break and go to next DNA sequence
                        break
                if not motif_found:
                    is_motif = False
                    break
            # After all sequences are parsed, if still true, neighbor was found in all sequences, so add that motif to "motifs_found" set
            if is_motif:
                motifs_found.add(neighbor)
        patterns.update(motifs_found)
    return ' '.join(list(set(patterns)))

# Find patterns that form (L,t)-clumps i.e., occurring t times within an L-nucleotide stretch
def find_clumps(seq, k, L, t):
    patterns = []
    for i in range(len(seq) - int(L) + 1):
        # Scanning window of length L
        window = seq[i:i + int(L)]
        # Get freqTable dictionaries of k-mers within L-length windows of sequence
        freqMap = freq_map(window, k)
        # For each item in dictionary
        for j in freqMap:
            # If value of key (counts) is t or higher, add to Patterns list
            if freqMap[j] >= int(t) and j not in patterns:
                patterns.append(j)
    patterns_nocommas = ' '.join(map(str, patterns))
    return patterns_nocommas

# Faster/more efficient code for longer genomes
def find_clumps_V2(seq, k, L, t):
    patterns = set()
    # Make one freqMap dictionary, default key values = 0
    freqMap = defaultdict(int)
    # Add each pattern to freqMap for first window
    for nuc in range(int(L)):
        pattern = seq[nuc:nuc + int(k)]
        freqMap[pattern] += 1
        if freqMap[pattern] == t:
            patterns.add(pattern)
    # Slide the window and update the frequency table incrementally
    for nuc in range(1, len(seq) - int(L) + 1):
        old_pattern = seq[nuc - 1:nuc - 1 + int(k)]
        # Decrement the k-mer that exited to avoid double counts in the current window
        freqMap[old_pattern] -= 1
        new_pattern = seq[nuc + int(L) - int(k):nuc + int(L)]
        # Increment the new k-mer that entered window
        freqMap[new_pattern] += 1
        if freqMap[new_pattern] == t:
            patterns.add(new_pattern)
    return list(patterns)

# Take a theoretical k-mer pattern: "How close is this pattern to a matching k-mer (present in all 
# DNA sequences)?"
# Let's say the pattern "AAA" is not very close to what exists, giving a very high total Hamming distance
# But maybe "CGT" has a very low total Hamming distance, i.e., a plausible k-mer with minimal mismatches 
# between sequences (say 2 or fewer mismatches in each sequence)
def sum_hamming_distances(pattern, Dna):
    # Cumulative sum of hamming distances
    distance_sum = 0
    k = len(pattern)
    # For each sequence in Dna
    for seq in Dna:
        min_distance = float('inf')
        # Substring sequence - k
        for nuc in range(len(seq) - int(k) + 1):
            # Scanning window of k-mers
            kmer = seq[nuc:nuc + int(k)]
            # Hamming distance between current k-mer and all theoretically possible k-mer patterns
            distance = hamming_distance_V2(pattern, kmer)
            # Create a running minimum Hamming distance for current DNA sequence
                # If current distance < current minimum, update minimum
            if distance < min_distance:
                min_distance = distance
        # Update cumulative sum of minimum distances for each sequence
        distance_sum += min_distance
    return distance_sum

# Find all k-mer Patterns ("median string") that minimize d(Pattern, Dna)
# This function finds the k-mers that minimize the total Hamming distance across a set of DNA sequences
# This function checks the Hamming distance between all k-mers in each DNA sequence and all possible k-mers, and finds the k-mer that minimizes the hamming distance sum for all sequences
# Thus, pattern candidates are determined, corresponding to 
def median_string(Dna, k):
    min_distance = float('inf')
    medians = []
    # For all possible variations of k nucleotides (patterns)
        # Iterates over 4 nucleotides and creates all theoretically possible k-length pattern combos
    for pattern in [''.join(p) for p in itertools.product('ACGT', repeat = int(k))]:
        # Update min_distance generated by Sum_Distances(pattern, Dna) using each pattern
        if sum_hamming_distances(pattern, Dna) < min_distance:
            min_distance = sum_hamming_distances(pattern, Dna)
            medians = [pattern]
        elif sum_hamming_distances(pattern, Dna) == min_distance:
            medians.append(pattern)
    return medians

# Find profile-most probable k-mer given a DNA sequence, integer k, and Profile matrix
def profile_most_probable(seq, k, profile):
    max_prob = -1
    most_probable_kmer = ""
    for i in range(len(seq) - int(k) + 1):
        # For each kmer of scanning window size k
        kmer = seq[i:i + int(k)]
        prob = 1
        for nuc in range(int(k)):
            # If current nucleotide == A, C, G or T
            if kmer[nuc] == 'A':
                # Set probability = corresponding probability in Profile
                # As we iterate through to k-th nucleotide in kmer, use that column in Profile
                prob *= profile[0][nuc]
            elif kmer[nuc] == 'C':
                prob *= profile[1][nuc]
            elif kmer[nuc] == 'G':
                prob *= profile[2][nuc]
            elif kmer[nuc] == 'T':
                prob *= profile[3][nuc]
            # kmer now has an associated probability
        if prob > max_prob:
            # Update the maximum probability found and assign the corresponding kmer to moxt_probable_kmer
            max_prob = prob
            most_probable_kmer = kmer
    return most_probable_kmer

# Build a profile matrix Profile for a given k-mer
    # Takes a list of motifs as input
def build_profile_matrix(motifs):
    # Length of each kmer/motif or row in the matrix
    k = len(motifs[0])
    # Initialize profile matrix with 4 rows (A,C,G,T) of lists containing k zeroes
    profile = [[0] * k for _ in range(4)]
    num_rows = len(motifs)
    # For each column
    for nuc in range(k):
        # Initialize dictionary of nucleotide counts
        count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        # For each row
        for motif in motifs:
            # Increment nucleotide counts in dictionary
                # motif[nuc] accesses the element of motif at the nuc-th column
                # count[motif[nuc]] is like count['A'] where motif[nuc] becomes the key whose value gets incremented
            count[motif[nuc]] += 1
        # For each row of Profile matrix
        for j in range(4):
            # Get value of j-th key (nucleotide) in {count} and divide by total number of rows in Motifs
                # the Profile rows are consistently 'A,C,G,T' so row 0 corresponds to A, etc.
            profile[j][nuc] = count['ACGT'[j]] / num_rows
    return profile

# Build a profile matrix Profile for a given k-mer
    # Takes a list of motifs as input
    # V2: pseudocounts using Laplace's Rule of Succession
def build_profile_matrix_V2(motifs):
    # Length of each kmer/motif or row in the matrix
    k = len(motifs[0])
    # Initialize profile matrix with 4 rows (A,C,G,T) of lists containing k zeroes
    profile = [[0] * k for _ in range(4)]
    num_rows = len(motifs)
    # For each column
    for nuc in range(k):
        # Initialize dictionary of nucleotide counts
        # LAPLACE'S RULE OF SUCCESSION: This is the only line we need to change
        count = {'A': 1, 'C': 1, 'G': 1, 'T': 1}
        # For each row
        for motif in motifs:
            # Increment nucleotide counts in dictionary
                # motif[nuc] accesses the element of motif at the nuc-th column
                # count[motif[nuc]] is like count['A'] where motif[nuc] becomes the key whose value gets incremented
            count[motif[nuc]] += 1
        # For each row of Profile matrix
        for j in range(4):
            # Get value of j-th key (nucleotide) in {count} and divide by total number of rows in Motifs
                # the Profile rows are consistently 'A,C,G,T' so row 0 corresponds to A, etc.
            profile[j][nuc] = count['ACGT'[j]] / num_rows
    return profile

# Score a motif matrix, where each column is assigned a score and all columns are totaled
# A column's score corresponds to the number of least conserved/most mutated nucleotides across all rows
    # The lower the score, the better / most conserved
def score(motifs):
    consensus = ""
    # Number of matrix columns
    num_cols = len(motifs[0])
    # Number of matrix rows
    num_rows = len(motifs)
    # For each column
    for nuc in range(num_cols):
        # Initialize dictionary of nucleotide counts, which will update for every column in Motifs
        count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        # For each row
        for motif in motifs:
            # Increment each nucleotide in count
            count[motif[nuc]] += 1
        # Build the consensus motif using the most-occurring nucleotide in current column for all rows
        consensus += max(count, key=count.get)
    # Calculate score
    score = 0
    # For each row
    for motif in motifs:
        # For each column
        for nuc in range(num_cols):
            # If the nucleotide differs from the corresponding consensus necleotide, increment score
            if motif[nuc] != consensus[nuc]:
                score += 1
    return score

# Construct motif matrices using Profile, and compare their scores to find the one with the lowest score
    # 't' is number of strings in Dna
def greedy_motif_search(Dna, k, t):
    best_motifs = [string[:int(k)] for string in Dna]
    # Specify/Stay in the first Dna string: building each motif matrix from Motif_1
    for i in range(len(Dna[0]) - int(k) + 1):
        # List of k-mers for each motif matrix, where the first in the list is from the first string Dna_1
        # At first, [kmers] has only Motif_1 from Dna_1, which iterates as a scanning window creating each motif matrix
        kmers = [Dna[0][i:i + int(k)]]
        # For each string in Dna
        for string in range(1, int(t)):
            # Construct Profile using preceding k-mers
            profile = build_profile_matrix(kmers)
            # Append the profile-most probable kmer in the next string
            kmers.append(profile_most_probable(Dna[string], k, profile))
        if score(kmers) < score(best_motifs):
            best_motifs = kmers
    return ' '.join(best_motifs)

# Calculate the entropy of a list of numbers rounded to 2 decimal points
def calculate_entropy(numbers):
    ans = 0
    for num in numbers:
        if num == 0:
            # Creating rule that log_2(0) = 0
            ans -= 0
        else:
            ans -= num * math.log2(num)
    return round(ans, 2)

# Find the GC content (%) from a FASTA file
def GC_content(file):
    seq_IDs = []
    seqs = []
    seq_lengths = []
    GC_contents = []
    
    current_seq = ""
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            # If a header is reached
            if line.startswith(">"):
                # If empty, "False", doesn't evaluate, ID is appended (first line)
                # If data in current_seq, "True", evaluates and appends, then resets
                if current_seq:
                    seqs.append(current_seq)
                    current_seq = ""
                seq_IDs.append(line)
            # If a line is reached that isn't a header
            else:
                # Concatenate current sequence variable
                current_seq += line
        # Last line
        if current_seq:
            seqs.append(current_seq)
                
    for seq in seqs:
        seq_lengths.append(len(seq))
        GCs = 0
        for nuc in seq:
            if nuc == 'G' or nuc == 'C':
                GCs += 1
        GC_contents.append(GCs)
        
    GC_percentages = []
    for seq_len, GC_content in zip(seq_lengths, GC_contents):
        GC_percentages.append((GC_content / seq_len) * 100)
        
    ans = []
    max_GCs = 0
    for ID, GC_percentage in zip(seq_IDs, GC_percentages):
        if GC_percentage > max_GCs:
            max_GCs = round(GC_percentage, 6)
            max_GC_seqID = ID
    ans.append(max_GCs)
    ans.append(max_GC_seqID)
    
    return seqs, seq_IDs, seq_lengths, GC_contents, GC_percentages, ans