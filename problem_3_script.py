#!/Library/Frameworks/Python.framework/Versions/3.6/bin/python3
# problem_3_script.py
# Corban Swain, March 2018

# importing required packages
import numpy as np



## Loading in Data ##
# protein sequences
prots = [['G', 'I', 'G', 'V', 'D', 'D', 'V', 'I', 'K', 'L'],
         ['Q', 'L', 'G', 'P', 'K', 'D', 'R', 'L', 'S', 'I'],
         ['G', 'L', 'G', 'V', 'F', 'D', 'A', 'I', 'S', 'L'],
         ['Q', 'L', 'G', 'I', 'N', 'D', 'Q', 'I', 'S', 'M'],
         ['G', 'I', 'G', 'K', 'Q', 'N', 'E', 'V', 'S', 'L']]

# convert to numpy array for easier slicing
prots = np.array(prots)

# number of sequences
M = len(prots)

# amino acid (aa) backgroud frequencies
q_dict = {'A' : 0.073, 'C' : 0.025, 'D' : 0.050, 'E' : 0.061,
          'F' : 0.042, 'G' : 0.072, 'H' : 0.023, 'I' : 0.053,
	  'K' : 0.064, 'L' : 0.089, 'M' : 0.023, 'N' : 0.043,
	  'P' : 0.052, 'Q' : 0.040, 'R' : 0.052, 'S' : 0.073,
	  'T' : 0.056, 'V' : 0.063, 'W' : 0.013, 'Y' : 0.033}

# list of amino acids
aas = q_dict.keys()



## Part 1 ##
# function for counting the frequency of aas at a given position
# in a multiple sequence alignment (MSA)
def aa_freqs(seqs, pos):
    # number of sequences
    M = len(seqs)
    
    # constant to save a bit of computation time
    freq_per_count = 1.0 / M

    # initialize empty dictionary
    freq_dict = {}
    
    # looping through amino acids in each position
    for aa in seqs[:, pos]:
        # get the number of counts for that aa in the
        # dictionary, zero if not present
        current_count = freq_dict.get(aa, 0.0)

        # increase the count by one
        freq_dict[aa] = current_count + freq_per_count

    return freq_dict

# positions to analyze (zero indexed)
positions = [0, 4, 9]

# list of dictionaries to keep track of counts
freqs_out = []

# looping through queried positions
for pos in positions:
    # calculating the frequenciec for each aa at 'pos'
    # in the MSA
    freqs_out.append(aa_freqs(prots, pos))

# function for pretty printing (pp) the dictionary lists
# used throughout this script
def pp_dict_list(dict_list):
    for iPos, d in enumerate(dict_list):
        print('  Position %d:' % (positions[iPos] + 1))
        for key, val in d.items():
            print('    {:s} - {:6.4f}'.format(key, val))

# printing the calulated output
print('\nPart 1: AA Frequencies')
pp_dict_list(freqs_out)           



## Parts 2 and 3 ##
# useful math functions
fact = np.math.factorial
ln = np.math.log
exp = np.math.exp

# functions for calulating the probabaluty of observing a given
# amino acid at a specific frequency withing a set of aligned
# sequences. Formulated according to Ranganathan Lab - Note 109:
# http://systems.swmed.edu/rr_lab/Note109_files/Note109_v3.html
def p_observed(aa, freq, M):
    counts = round(M * freq)
    noncounts = M - counts
    q = q_dict[aa]
    temp_1 = fact(M) / (fact(counts) * fact(noncounts))
    temp_2 = q ** counts
    temp_3 = (1 - q) ** noncounts
    return temp_1 * temp_2 * temp_3

# relative entropy calulation
def rel_entropy(aa, freq, M):
    nonfreq = 1 - freq
    q = q_dict[aa]
    temp_1 = freq * ln(freq / q)
    temp_2 = nonfreq * ln(nonfreq / (1 - q))
    return temp_1 + temp_2

# uses Sterling's approximation    
def p_observed_approx(aa, freq, M):
    D = rel_entropy(aa, freq)
    return exp(-M * D)

# initializing ouput caches
p_obsv_out = []
entropy_out = []

# looping through each position
for obsv_freqs in freqs_out:
    # adding a new empty dictionary to each output cache
    p_obsv_out.append({})
    entropy_out.append({})

    # looping through each aa present
    for aa, freq in obsv_freqs.items():
        # calulating probability of observing each aa it its
        # present frequency
        p_obsv_out[-1][aa] = p_observed(aa, freq, M)

        # calulating the relative entropy assosiated with
        # observing each aa at its present frequency
        entropy_out[-1][aa] = rel_entropy(aa, freq, M)
        
# printing the calulated output
print('\nPart 2: Probabilities of Observation')        
pp_dict_list(p_obsv_out)
print('\nPart 3: Relative Entropy for Observation')        
pp_dict_list(entropy_out)



## Part 4 ##
print('\nPart 4: Meaning of Relative Entropy')        
print('''
  A relative entropy of zero represents the case when an amino acid
  at a given position in the multiple sequence alighment appears with
  frequency equivalent to the backgroud frequency for that amino acid.
''')


## Part 5 ##
# function for calulating the statistical coupling energy
def coupling_E(seqs, aa1, pos1, aa2, pos2):
    # number of protein sequences
    M = len(seqs)
    
    # frequency of aa2 at pos2 in the MSA
    freq = aa_freqs(seqs, pos2).get(aa2, 0.0)
    
    # subset of the sequences with aa1 at pos1
    enriched_seqs = [s for s in seqs if s[pos1] == aa1]
    enriched_seqs = np.array(enriched_seqs)

    # frequency of aa2 at pos2 in the subset
    coupled_freq = aa_freqs(enriched_seqs, pos2).get(aa2, 0.0)

    # calculating the statistical coupling energy
    temp1 = ln(p_observed(aa2, freq, M))
    temp2 = ln(p_observed(aa2, coupled_freq, len(enriched_seqs)))    
    return -1.0 / M * (temp1 - temp2)

# performing the calculation for the spcified inputs
ddG = coupling_E(prots, 'G', 0, 'L', 9)

# printing the calculated output
print('\nPart 5: Statistical Coupling Energy')        
print('  \Delta\DeltaG = %f' % ddG)
