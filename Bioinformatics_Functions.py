# .py file containing various Bioinformatics functions
import random
import itertools
from itertools import product
from collections import defaultdict
import matplotlib.pyplot as plt
import math


# Count the number of nucleotides
def nucleotide_count(seq):
    
    # Initialize dictionary
    nuc_counts = {'A':seq.count('A'), 'C':seq.count('C'), 'G':seq.count('G'), 'T':seq.count('T')}
    
    return nuc_counts


# Obtain the reverse complement of a DNA sequence
def reverse_complement(seq):
    
    # Assign complementary letters
    comp_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    
    # Replace each nucleotide with its complement
    comp_seq = ''.join(comp_dict[nuc] for nuc in seq)
    
    # Reverse the string
    return comp_seq[::-1]


# Randomly generate a DNA sequence of length n
def random_sequence(n):
    
    return ''.join(random.choice('ACGT') for _ in range(n))


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
    return sum([1 for nuc in range(len(seq) - len(pattern) + 1) if seq[nuc:nuc + len(pattern)] == pattern])


# Find locations of a k-mer pattern within a DNA sequence
def find_pattern_indices(seq, pattern):
    
    # Initialize list of indices
    indices = []
    # Length of pattern
    k = len(pattern)
    
    # For whole sequence - k
    for nuc in range(len(seq) - k + 1):
        
        # Checking for our pattern using a window of size k 
        if seq[nuc:nuc + k] == pattern:
            # If our pattern is found, append the index to indices list
            indices.append(nuc)
            
    return ' '.join(map(str, indices))


# Find indices of point mutations between two sequences
def hamming_distance(seq1, seq2):
    
    # Ensure lengths of sequences are equal
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be the same length")
        
    return sum([1 for nuc in range(len(seq1)) if seq1[nuc] != seq2[nuc]])


# Count the number of approximate pattern occurrences with d or fewer mismatches
def approx_pattern_count(seq, pattern, d):
    
    # Initialize list of indices
    positions = []
    
    # Sliding window
    for nuc in range(len(seq) - len(pattern) + 1):
        if hamming_distance(pattern, seq[nuc:nuc + len(pattern)]) <= int(d):
            positions.append(nuc)
            
    return len(positions)


# Find all approximate occurrences of a pattern in a string (d or fewer mismatches from a given pattern)
def approx_pattern_occurrences(pattern, seq, d):
    
    # Initialize list of indices
    indices = []
    k = len(pattern)
    
    # Sliding window (+ 1 because range() goes from 0 to n - 1)
    for i in range(len(seq) - k):
        if hamming_distance(pattern, seq[i:i + k]) <= int(d):
            indices.append(i)
    
    # Return indices list as a string separated by spaces
    return ' '.join(map(str, indices))


# Create a frequency map of k-mer pattern occurrences as a dictionary
def freq_map(seq, k):
    
    # Initialize final freqMap dict
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


# Retrieve the most frequent patterns in frequency map
def most_freq_patterns(seq, k):
    
    # Initialize final list of most frequenct patterns
    freqPatterns = []
    # Get freqMap dictionary
    freqMap = freq_map(seq, int(k))
    
    # For each pattern, count item in freqMap
    for pattern, count in freqMap.items():
        # Find patterns with max count values
        if count == max(freqMap.values()):
            # Append freqPatterns list with those patterns
            freqPatterns.append(pattern)
            
    freqPatterns.sort()
    return freqPatterns


# Find patterns that form (L,t)-clumps i.e., occurring t times within an L-nucleotide stretch
    # This function may not work properly
def find_clumps(seq, k, L, t):
    
    # Final k-mers that form (L,t)-clumps
    clump_kmers = []
    
    # For  each nuc in sequence - L
    for nuc in range(len(seq) - int(L)):
        # Window size = L
        window = seq[nuc:nuc + int(L)]
        # Make freq map of all kmers for each window L
        freqMap = freq_map(window, k)
        
        for kmer in freqMap:
            if freqMap[kmer] == int(t) and kmer not in clump_kmers:
                clump_kmers.append(kmer)
                
    return ' '.join(clump_kmers)


# Faster/more efficient code for longer genomes
    # Uses one freqMap and adds k-mers as it goes, instead of making several maps
def find_clumps_V2(seq, k, L, t):
    
    # Final k-mers that form (L,t)-clumps
        # Set to avoid duplicates
    clump_kmers = set()
    
    # Make one single freqMap with default key values = 0
    freqMap = defaultdict(int)
    
    # Add each k-mer to freqMap for first window L
    # For each nuc in first window L
    for nuc in range(int(L)):
        # k-mer window
        kmer = seq[nuc:nuc + int(k)]
        # Increment k-mer in freqMap
        freqMap[kmer] += 1
        
        # Add k-mer if = t occurrences
        if freqMap[kmer] == t:
            clump_kmers.add(kmer)
            
    # Slide the window one nuc and update the frequency table incrementally
    for nuc in range(1, len(seq) - int(L) + 1):
        # k-mer that exited to the left
        old_kmer = seq[nuc - 1:nuc - 1 + int(k)]
        # Decrement the k-mer that exited to avoid double counts in the current window
        freqMap[old_kmer] -= 1
        # New k-mer at the far right of the current window
        new_kmer = seq[nuc + int(L) - int(k):nuc + int(L)]
        # Increment the new k-mer that entered window from the right
        freqMap[new_kmer] += 1
        
        if freqMap[new_kmer] == t:
            clump_kmers.add(new_kmer)
            
    # Get answer as a list
    result = list(clump_kmers)
    # Sort list alphabetically
    result.sort()
    
    return result


# The running difference between occurrences of G and C in genome
def skew(genome):
    
    # Initialize cumulative sum
    skew_vals = [0]
    
    for nuc in genome:
        if nuc == 'C':
            # Subtract one from last value
            skew_vals.append(skew_vals[-1] - 1)
        elif nuc == 'G':
            # Add one to last value
            skew_vals.append(skew_vals[-1] + 1)
        else:
            # Duplicate the last value
            skew_vals.append(skew_vals[-1])
            
    return skew_vals


# Find the minimum skew position(s), where the reverse half-strand ends and the forward half-strand begins
def min_skew_positions(seq):
    
    # Get skew(seq)
    skew_vals = skew(seq)
    # Find minimum value(s) once for higher efficiency
    min_skew_val = min(skew_vals)
    
    # Obtain all minimum skew value indices
    return [idx for idx, skew in enumerate(skew_vals) if skew == min_skew_val]


# Create skew plot from skew()
def plot_skew(genome):
    
    # y-axis
    y = skew(genome)
    # x-axis
        # + 1 accounts for extra 0 at the beginning of skew_vals
    x = range(len(genome) + 1)    
    
    plt.plot(x, y, color='black')
    plt.xlabel('Position (i)')
    plt.ylabel('Skew')
    plt.title('Skew Plot')
    plt.grid(True)
    
    plt.show()

    
# Generate all theoretically possible neighbors for a given pattern
def neighbors(pattern, d):
    """
    Use product() to find all possible mismatched indices and all possible nucleotide swaps.
    Then, use zip() to convert each mismatch into a neighbor.
    """
    
    # Case: d = 0
    if d == 0:
        return {pattern}
    
    # Case: pattern is 1 letter
    if len(pattern) == 1:
        return set('ACGT')

    # Case: d > 0 and pattern > 1 letter
    # Initialize neighborhood set (avoid duplicates)
    neighborhood = set()
    
    # All combinations of indices in pattern that could be mismatched
    for mismatched_indices in product(range(len(pattern)), repeat=d):
    # Result if pattern length=3 and d=2:
    # (0,0), (0,1), (0,2)
    # (1,0), (1,1), (1,2)
    # (2,0), (2,1), (2,2)

        # Nucleotide substitutions that will be applied: (A,A), (A,C), etc. if d=2
        for mismatches in product('ACGT', repeat=d):
        # Result if pattern length=3 and d=2:
        # (A,A), (A,C), (A,G), (A,T)
        # (C,A), (C,C), ...
        # (G,A), ...
        # (T,A), ...

            # Convert pattern to a list so it is mutable
            neighbor = list(pattern)
            
            # Combine elements from mismatched_indices and mismatches as tuples: (0, A), (0, C), etc.
            for idx, mismatch in zip(mismatched_indices, mismatches):
                # Create mismatch substitutions
                neighbor[idx] = mismatch
            # All possible mismatches
            neighborhood.add(''.join(neighbor))
            
    neighborhood = list(neighborhood)
    return sorted(neighborhood)
    
# Generates all possible neighbors of a pattern with d or fewer mismatches
    # More efficient approach
def neighbors_V2(pattern, d):
    """
    Use recursion to generate all possible neighbors with d or fewer mismatches for the suffix of pattern (i.e. pattern[1:])
    Check hamming dist. between suffix neighbors and pattern. If < d, introduce 4 new neighbors by changing first nucleotide
        If = d, keep original nucleotide the same and add the neighbor
    """
    
    # Case: d = 0
    if d == 0:
        return {pattern}
    
    # Case: pattern is 1 nucleotide
    if len(pattern) == 1:
        return set('ACGT')

    # Case: d > 0 and pattern > 1 nucleotide
    # Initialize neighborhood set (avoid duplicates)
    neighborhood = set()
    
    # Use recursion to iterate over all possible mismatches for the suffix (excluding the first nucleotide)
        # Generate all possible neighbors of the suffix with d or fewer mismatches
            # Example: pattern = 'ACGT'
            # When recursion reaches the last nucleotide ('T'), it calls itself with a pattern = 1 nucleotide and generates A,C,G,T
            # Then it calls itself with pattern = 'GT', and checks hamming_dist('T', suffix neighbors ('A,C,G,T'))
                # This generates suffix neighbors ('GA', 'GC', etc.)
            # Then it calls iself with pattern = 'CGT' and generates ('CAT', 'CAA', etc.)
    suffix_neighbors = neighbors_V2(pattern[1:], d)

    for neighbor in suffix_neighbors:
        
        # If the Hamming distace between the pattern and the suffix neighbor is < d, the neighbor differs by at most d mismatches
        if hamming_distance(pattern[1:], neighbor) < int(d):
            
            for nuc in 'ACGT':
                # Concatenate each of the four nucleotides to the suffix of the pattern, handling mismatches for the first nucleotide
                neighborhood.add(nuc + neighbor)
                
        # If Hamming dist. between pattern and suffix neighbor = d, it is already a neighbor using the original first nucleotide
            # the Hamming dist. cannot be > d because of the nature of the recursive call
        else:
            neighborhood.add(pattern[0] + neighbor)
            
    # Get neighborhood as a list
    result = list(neighborhood)
    # Sort alphabetically
    result.sort()
            
    return result

# Find the most frequently occurring patterns in a sequence with d or fewer mismatches
def most_freq_mismatched_patterns(seq, k, d):
    """
    Use neighbors() to generate all possible neighbors of every k-mer in sequence
    Then, find most frequent neighbors
    """
    
    # Initialize freqMap dict of neighbors
    neighbor_freqMap = {}
    # Initialize final list of most frequent neighbors
    most_freq_neighbors = []
    
    # Find all possible neighbors of every k-mer in sequence
    for nuc in range(len(seq) - int(k) + 1):
        # Sliding window of k-mers
        kmer = seq[nuc:nuc + int(k)]
        # Get neighborhood of all possible neighbors of current k-mer with up to d mismatches
            # Current k-mer acts as the given pattern for neighbors()
        neighborhood = neighbors(kmer, int(d))
        
        # Convert neighbors to a freqMap dictionary
        for neighbor in neighborhood:
            # Set all patterns = 0 in freqMap
            if neighbor not in neighbor_freqMap:
                neighbor_freqMap[neighbor] = 0
            # If pattern is found, increment in freqMap
            else:
                neighbor_freqMap[neighbor] += 1
                
    # Sort neighbor_freqMap alphabetically
    temp = list(neighbor_freqMap.keys())
    temp.sort()
    neighbor_freqMap = {neighbor:neighbor_freqMap[neighbor] for neighbor in temp}
    
    # Find most frequent k-mer(s) with d or fewer mismatches
    most_freq_neighbors = [neighbor for neighbor, count in neighbor_freqMap.items() if count == max(neighbor_freqMap.values())]
    # Most frequent neighbor(s) as a space-separated string AND neighbor_freqMap
    return ' '.join(most_freq_neighbors), neighbor_freqMap


# Finding the most frequent mismatched k-mers and obtaining their reverse complements
    # This function returns the most frequent mismatched k-mer(s) and the reverse complement(s)
def most_freq_mismatched_patterns_revcomp(seq, k, d):
    
    patterns = []
    freqMap = {}
    
    for nuc in range(len(seq) - int(k) + 1):
        pattern = seq[nuc:nuc + int(k)]
        # Generate the pattern's reverse complement
        pattern_rc = reverse_complement(pattern)
        neighborhood = neighbors(pattern, int(d))
        # Generate a neighborhood of reverse complements
        neighborhood_rc = [reverse_complement(n) for n in neighborhood]
        
        for neighbor in neighborhood + neighborhood_rc:
            
            if neighbor not in freqMap:
                freqMap[neighbor] = 1
            else:
                freqMap[neighbor] += 1

    max_freq = max(freqMap.values())
    patterns = [pattern for pattern, freq in freqMap.items() if freq == max_freq]
    
    return ' '.join(patterns)

    
# Find (k,d)-motifs in a given collection of DNA sequences
    # (k,d)-motifs are k-mers that appear in ALL sequences with d or fewer mismatches
    # This function takes a list of strings (DNA sequences)
def motif_enumeration(Dna: list, k: int, d: int) -> str:
    """
    Generate all possible neighbors for every k-mer in the first Dna sequence.
    Then, for each neighbor, check the hamming dist. between neighbor and every k-mer in each sequence.
    If neighbor has a hamming dist. <= d for a k-mer in every sequence, it is a (k,d)-motif.
    Do this for all neighbors.
    """
    
    # Initialize final set of (k,d)-motifs (avoid duplicates)
    motifs = set()
    
    # For first DNA sequence
    for nuc in range(len(Dna[0]) - int(k) + 1):
        # Sliding window k-mer
        kmer = Dna[0][nuc:nuc + int(k)]
        # Get all possible neighbors of each k-mer in first sequence
        kmer_neighbors = neighbors_V2(kmer, d)
        # Initialize "motifs_found" set for (k,d)-motifs found in all sequences
        motifs_found = set()
        
        # For each neighbor in the set of all possible neighbors for all k-mers first sequence
        for neighbor in kmer_neighbors:
            is_motif = True
            
            # Check DNA sequences for (k,d)-motifs with the current neighbor
            for seq in Dna:
                motif_found = False
                
                # For first DNA sequence
                for nuc in range(len(seq) - int(k) + 1):
                    # Check mismatches between each k-mer window and first neighbor
                    # If hamming_distance <= d
                    if hamming_distance(neighbor, seq[nuc:nuc + int(k)]) <= int(d):
                        # Motif found
                        motif_found = True
                        # Break and check the next DNA sequence
                        break
                
                # If no motif was found in the current sequence for the current neighbor
                    # motif_found is still false, move to next neighbor
                if motif_found == False:
                    is_motif = False
                    break
                    
            # If still true, that neighbor was found in all sequences, so add neighbor to "motifs_found"
            if is_motif:
                motifs_found.add(neighbor)
                
        motifs.update(motifs_found)

    # Remove duplicates
    motifs = list(set(motifs))
    # Sort alphabetically
    motifs.sort()
    
    return ' '.join(motifs)

# Take a theoretical pattern: "How close is this pattern to a matching k-mer (present in all DNA sequences)?"
# Let's say the pattern "AAA" is not very close to what exists, giving a very high total Hamming distance
# But maybe "CGT" has a very low total Hamming distance
    # i.e., a plausible k-mer with minimal mismatches between sequences (say 2 or fewer mismatches in each sequence)
def sum_hamming_distances(pattern, Dna):
    """
    Takes a pattern, and a list of Dna sequences
    Finds k-mers in each sequence with the minimum hamming dist. between the given pattern
    Returns a cumulative sum of the total minimum hamming dist. for all sequences
    """
    
    # Cumulative sum of hamming distances
    total_cum_sum = 0
    k = len(pattern)
    
    # For each sequence in Dna
    for seq in Dna:
        # Initialize minimum distance to infinity
        min_distance = float('inf')
        
        # For each nuc in sequence
        for nuc in range(len(seq) - int(k) + 1):
            # Scanning window k-mer
            kmer = seq[nuc:nuc + int(k)]
            # Hamming distance between current k-mer and the given pattern
            distance = hamming_distance(pattern, kmer)
            
            # Create a running minimum Hamming distance for current DNA sequence
                # If current hamming dist. < current minimum, update minimum
            if distance < min_distance:
                min_distance = distance
                
        # Update cumulative sum of minimum distances for each sequence
        tot_cum_sum += min_distance
        
    return tot_cum_sum


# Find all k-mers ("median strings") that minimize sum_hamming_distance
def median_string(Dna, k):
    """
    Takes a list of DNA sequences, and an integer k
    Uses sum_hamming_distances() to find the k-mer(s) with the smallest total minimum hamming dist. in Dna
    For each k-mer in Dna, it becomes "pattern" in sum_hamming_distances
    This function returns the k-mer(s) than minimize the total from sum_hamming_distances
    """
    
    # Initialize minimum distance to infinity
    min_distance = float('inf')
    # Initialize final list of median strings
    medians = []
    
    # Generate a list of all possible k-mers
    kmers = [kmer for kmer in itertools.product('ACGT', repeat = int(k))]
    
    # For all k-mers
    for kmer in kmers:
        
        # If sum_hamming_distances() < min_distance
        if sum_hamming_distances(kmer, Dna) < min_distance:
            # Update min_distance
            min_distance = sum_hamming_distances(kmer, Dna)
            # Add k-mer to medians
            medians = [kmer]
            
        # If sum_hamming_distances() = min_distance
        elif sum_hamming_distances(kmer, Dna) == min_distance:
            # Add k-mer to medians
            medians.append(kmer)
            
    return medians


# Find profile-most probable k-mer given a DNA sequence, integer k, and Profile matrix
def profile_most_probable(seq, k, profile):
    """
    Takes a given DNA sequence, integer k, and a profile matrix
    Finds the k-mer with the highest probability from the profile matrix (i.e. the profile-most-probable k-mer)
    Rows in Profile are: A, C, G, T
    """
    
    # Initialize max probability of -1
    max_prob = -1
    # Initialize empty most probable k-mer
    most_probable_kmer = ""
    
    # For each nucleotide in seq - k
    for nuc in range(len(seq) - int(k) + 1):
        # Scanning window k-mer
        kmer = seq[nuc:nuc + int(k)]
        # Initialize probability = 1
        prob = 1
        
        # For each nucleotide in current k-mer
        for nuc in range(int(k)):
            
            # If current nucleotide == A, C, G or T
            if kmer[nuc] == 'A':
                # Set probability = corresponding nucleotide (row) probability in Profile
                    # As we iterate through to k-th nucleotide in kmer, use that column in Profile
                prob *= profile[0][nuc]
                
            elif kmer[nuc] == 'C':
                prob *= profile[1][nuc]
                
            elif kmer[nuc] == 'G':
                prob *= profile[2][nuc]
                
            elif kmer[nuc] == 'T':
                prob *= profile[3][nuc]
            # k-mer now has an associated probability, given the Profile matrix
            
        # If k-mer's probability > max_prob    
        if prob > max_prob:
            # Update the maximum probability found and assign the corresponding kmer to moxt_probable_kmer
            max_prob = prob
            most_probable_kmer = kmer
            
    return most_probable_kmer


# Build a profile matrix Profile for a given k-mer
    # Takes a list of motifs as input
def build_profile_matrix(motifs):
    """
    Takes a list of motifs (motif matrix)
    Gets the counts for every nucleotide in motif matrix
    Converts the nucleotide counts to frequencies/probabilities
    Returns profile matrix
    """
    
    # Length of each row/motif in the matrix
    k = len(motifs[0])
    # Initialize profile matrix with 4 rows (A,C,G,T), containing k zeroes each (number of nucleotides)
    profile = [[0] * k for _ in range(4)]
    # Number of rows in motif matrix
    num_rows = len(motifs)
    
    # For each column/nucleotide in motif matrix
    for nuc in range(k):
        # Initialize dictionary of nucleotide counts
        count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        
        # For each row/motif
        for motif in motifs:
            # Increment nucleotide counts in dictionary
            count[motif[nuc]] += 1
            
        # For each row of Profile matrix
        for j in range(4):
            # Convert the nucleotide counts to frequencies/probabilities
                # Gets value of j-th key (nucleotide) in {count} and divides by total number of rows in motif matrix
            profile[j][nuc] = count['ACGT'[j]] / num_rows
            
    return profile


# Build a profile matrix Profile for a given k-mer
    # This version uses Laplace's Rule of Succession to eliminate zeroes from the profile matrix (pseudocounts)
def build_profile_matrix_V2(motifs):
    """
    Takes a list of motifs (motif matrix)
    Gets the counts for every nucleotide in motif matrix
    Converts the nucleotide counts to frequencies/probabilities
    Returns profile matrix
    """
    
    # Length of each row/motif in the matrix
    k = len(motifs[0])
    # Initialize profile matrix with 4 rows (A,C,G,T), containing k zeroes each (number of nucleotides)
    profile = [[0] * k for _ in range(4)]
    # Number of rows in motif matrix
    num_rows = len(motifs)
    
    # For each column/nucleotide in motif matrix
    for nuc in range(k):
        # Initialize dictionary of nucleotide counts
            # LAPLACE'S RULE OF SUCCESSION: all counts start from 1
        count = {'A': 1, 'C': 1, 'G': 1, 'T': 1}
        
        # For each row/motif in motif matrix
        for motif in motifs:
            # Increment nucleotide counts in dictionary
            count[motif[nuc]] += 1
            
        # For each row of Profile matrix
        for row in range(4):
            # Convert the nucleotide counts to frequencies/probabilities
                # Gets value of row-th key (nucleotide) in {count} and divides by total number of rows in motif matrix
            profile[row][nuc] = count['ACGT'[row]] / num_rows
            
    return profile


# Score a motif matrix, where each column is assigned a score and all columns are totaled
    # If all motifs are highly conserved, the score will be low
# A column's score corresponds to the conservation of nucleotides across all rows for each column
    # The lower the score, the better / most conserved
        # e.g., if every motif contains a 'G' for the first column/nucleotide, score = 0, maximum conservation
        # If 8/10 motifs contain a 'G', and two contain a different letter, score = 2
def score(motifs):
    """    
    Takes a given motif matrix
    Creates a dictionary of nucleotide counts for each column across all motifs
        If the motifs have length 10 nucleotides each, we will end up with 10 dictionaries of nucleotide counts for each column
    Returns a total score using a consensus string
        Consensus string is a string containing the most conserved nucleotide for each column
    """
    
    # Initialize consensus string (most frequent nucleotide for each column)
    consensus = ""
    # Number of columns/nucleotides in motif matrix
    num_cols = len(motifs[0])
    # Number of rows/motifs in motif matrix
    num_rows = len(motifs)
    
    # For each column/nucleotide in motif matrix
    for nuc in range(num_cols):
        # Initialize dictionary of nucleotide counts
            # Sets all nucleotide counts = 0
        count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        
        # For each row/motif in motif matrix
        for motif in motifs:
            # Increment each nucleotide in count
            count[motif[nuc]] += 1
        # Find the most frequent nucleotide in column dictionary add it to the consensus string
        consensus += max(count, key=count.get)
        # Move to the next column and create a new dictionary | repeat
    
    # Now we have a consensus string with the most conserved nucleotide at each position for all motifs
    # Initialize score
    score = 0
    
    # For each row/motif
    for motif in motifs:
        
        # For each column/nucleotide
        for nuc in range(num_cols):
            
            # If the nucleotide differs from the nucleotide in the consensus string, increment score
            if motif[nuc] != consensus[nuc]:
                score += 1
    
    # Return total score for the given motif matrix
    return score


# Given a collection of DNA sequences, find the best motifs
    # i.e., motifs present in ALL sequences with the lowest score/fewest mismatches compared to a consensus string/motif
    # The consensus string represents the most common/conserved nucleotide at each position
        # Each nucleotide position across all sequences will have a "most common/conserved" nucleotide. That is the consensus string
def greedy_motif_search(Dna, k, t):
    """
    Takes a list of DNA sequences, an integer k, and an integer t (number of sequences in Dna)
    Initializes "best_motifs" using the first k-mer from each sequence
    Starts by building a profile matrix with the first k-mer in the first DNA sequence
        Creates a motif matrix by finding the profile-most-probable k-mer in all remaining sequences
        Calculates the score of this "motif matrix"
        Repeats with second k-mer in first sequence. Compares the scores continuously
    Returns a list of best motifs (lowest score)
    """
    
    # Initialize final motif matrix using the first k-mer from each sequence
    best_motifs = [seq[:int(k)] for seq in Dna]
    
    # For nucleotide in first sequence - k
    for nuc in range(len(Dna[0]) - int(k) + 1):
        # Scanning window k-mer        
        # Initializes a list of k-mers starting with the first k-mer in the first sequence
        kmers = [Dna[0][nuc:nuc + int(k)]]
        
        # For each remaining DNA sequence
        for seq in range(1, int(t)):
            # Build a profile matrix using a single k-mer
            profile = build_profile_matrix_V2(kmers)
            # Find the profile-most-probable k-mer in all remaining sequences and add them to kmers
            kmers.append(profile_most_probable(Dna[seq], k, profile))
            # We now have a "motif matrix"
            
        # Score the motif matrix
            # If current score is < best_motifs
        if score(kmers) < score(best_motifs):
            # Update best_motifs
            best_motifs = kmers
            
        # Build a new profile matrix using the next k-mer in the first sequence | repeat
            
    return ' '.join(best_motifs)

# Calculate the entropy of a list of numbers rounded to 2 decimal points
def calculate_entropy(numbers):
    
    # Initialize final answer = 0
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


def Needleman_Wunsch_global_alignment(seq1, seq2):
    """
    Needleman-Wunsch algorithm for globally aligning two sequences
    Uses score and traceback matrices
    Works when seq2 is something like 'AC-GCTACCCT-AAAGT-'
        Where each '-' represents >= 1 deletion from seq1
    """
    
    # Get the lengths of each sequence
    len1 = len(sequence_one)
    len2 = len(sequence_two)

    # Initialize score and traceback matrices with zeroes
    score = [[0] * (len2 + 1) for i in range(len1 + 1)]
    traceback = [[0] * (len2 + 1) for i in range(len1 + 1)]

    # Boundary conditions: one extra row and column at the beginning of the matrix
    # For all columns in the first row (this should say all rows in the first column, oops)
    for i in range(1, len1 + 1):
        # Fill in the first column with boundary condition (-1, -2, -3, etc.)
        score[i][0] = score[i - 1][0] - 1
        # Boundary condition for first column of traceback (left = 2)
        traceback[i][0] = 2
        
    # First rows in the first column (this should say for all columns in the first row, oops)
    for j in range(1, len2 + 1):
        # Fill in the first row with boundary condition
        score[0][j] = score[0][j - 1] - 1
        # Boundary condition for first row of traceback (up = 3)
        traceback[0][j] = 3

    # Fill score matrix and update traceback (left to right, top to bottom)
    # For each row
    for i in range(1, len1 + 1):
        
        # For each column
        for j in range(1, len2 + 1):
            
            # If there's a match (starting at seq1[0] and seq2[0])
            if sequence_one[i - 1] == sequence_two[j - 1]:
                # Diagonal + 1
                score[i][j] = score[i - 1][j - 1] + 1
                # Update traceback with 'diagonal' (=1)
                traceback[i][j] = 1
                
            # If there's a gap/mismatch
            else:
                # If above >= left
                if score[i - 1][j] >= score[i][j - 1]:
                    # Above - 1
                    score[i][j] = score[i - 1][j] - 1
                    # Update traceback with 'up' (=2)
                    traceback[i][j] = 2
                    
                # If left >= above
                else:
                    # Left - 1
                    score[i][j] = score[i][j - 1] - 1
                    # Update traceback with 'left' (=3)
                    traceback[i][j] = 3

    # Use traceback to find alignment and missing nucleotides
    seq1_aligned = []
    seq2_aligned = []
    
    # Start from bottom right corner of the matrix (end of both sequences)
    nuc1 = len1
    nuc2 = len2

    # Until both sequences are finished
    while nuc1 > 0 or nuc2 > 0:
        
        # If diagonal (either a match or a '-')
        if nuc1 > 0 and nuc2 > 0 and traceback[nuc1][nuc2] == 1:
            # Add current nucleotide from both sequences
            seq1_aligned.append(sequence_one[nuc1 - 1])
            seq2_aligned.append(sequence_two[nuc2 - 1])
            # Move up and to the left
            nuc1 -= 1
            nuc2 -= 1
            
        # If above/up (gap in seq2)
        elif nuc1 > 0 and (nuc2 == 0 or traceback[nuc1][nuc2] == 2):
            # Add nucleotide from seq1
            seq1_aligned.append(sequence_one[nuc1 - 1])
            # Add '-' from seq2
            seq2_aligned.append('-')
            # Move up
            nuc1 -= 1
            
        # If left (gap in seq1)
        else:
            # Add '-' to seq1
            seq1_aligned.append('-')
            # Add nucleotide from seq2
            seq2_aligned.append(sequence_two[nuc2 - 1])
            # Move left
            nuc2 -= 1

    # Reverse the aligned sequences
    seq1_aligned.reverse()
    seq2_aligned.reverse()

    # Append list with seq1 nucleotide if corresponding seq2 nucleotide == '-'
    missing_nucs = [nuc1 for nuc1, nuc2 in zip(seq1_aligned, seq2_aligned) if nuc2 == '-']
    # Concatenate into a string
    result = ''.join(missing_nucs)
    
    # Remove '-'s
    return result.replace("-", "")
