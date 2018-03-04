#!/Library/Frameworks/Python.framework/Versions/3.6/bin/python3

import numpy as np
import pprint as pp

## Loading in Data ##
# protein sequences
prots = [['G', 'I', 'G', 'V', 'D', 'D', 'V', 'I', 'K', 'L'],
         ['Q', 'L', 'G', 'P', 'K', 'D', 'R', 'L', 'S', 'I'],
         ['G', 'L', 'G', 'V', 'F', 'D', 'A', 'I', 'S', 'L'],
         ['Q', 'L', 'G', 'I', 'N', 'D', 'Q', 'I', 'S', 'M'],
         ['G', 'I', 'G', 'K', 'Q', 'N', 'E', 'V', 'S', 'L']]
prots = np.array(prots)

M = len(prots)

# amino acid backgroud frequencies
q_dict = {'A' : 0.073, 'C' : 0.025, 'D' : 0.050, 'E' : 0.061,
          'F' : 0.042, 'G' : 0.072, 'H' : 0.023, 'I' : 0.053,
	  'K' : 0.064, 'L' : 0.089, 'M' : 0.023, 'N' : 0.043,
	  'P' : 0.052, 'Q' : 0.040, 'R' : 0.052, 'S' : 0.073,
	  'T' : 0.056, 'V' : 0.063, 'W' : 0.013, 'Y' : 0.033}

# list of amino acids
aas = aa_freq_dict.keys()

## Part 1 ##
# positions to analyze (zero indexed)
positions = [0, 4, 9]

# list of dictionaries to keep track of counts
freqs_out = []

# constant to save a bit of computation time
freq_per_count = 1.0 / M
             
# looping through queried positions
for pos in positions:
    # temporary dictionary to store counts of each aa
    # in the MSA
    temp_dict = {}

    # looping through amino acids in each position
    for aa in prots[:, pos]:
        # get the number of counts for that aa in the
        # dictionary, zero if not present
        current_count = temp_dict.get(aa, 0.0)

        # increase the count by one
        temp_dict[aa] = current_count + freq_per_count

    # appending the dict to the ouput cache
    freqs_out.append(temp_dict)

# Printing the calulated output
print('\nPart 1: AA Frequencies')        
pp.pprint(freqs_out)

## Part 2 and 3 ##
# useful math functions
fact = np.math.factorial
ln = np.math.log
exp = np.math.exp

# precalulate to save comp. time
M_fact = fact(M)

# functions for calulating the probabaluty of observing a given
# amino acid at a specific frequency withing a set of aligned
# sequences. Formulated according to Ranganathan Lab - Note 109:
# http://systems.swmed.edu/rr_lab/Note109_files/Note109_v3.html
def p_observed(aa, freq):
    counts = round(M * freq)
    noncounts = M - counts
    q = q_dict[aa]
    
    temp_1 = M_fact / (fact(counts) * fact(noncounts))
    temp_2 = q ** counts
    temp_3 = (1 - q) ** noncounts
    return temp_1 * temp_2 * temp_3

def rel_entropy(aa, freq):
    nonfreq = 1 - freq
    q = q_dict[aa]
    
    temp_1 = freq * ln(freq / q)
    temp_2 = nonfreq * ln(nonfreq / (1 - q))
    return temp_1 + temp_2

# uses sterlings approximation    
def p_observed_approx(aa, freq):
    D = rel_entropy(aa, freq)
    return exp(-M * D)

# looping through each position
p_obsv_out = []
entropy_out = []
for obsv_freqs in freqs_out:
    p_obsv_out.append({})
    entropy_out.append({})
    for aa, freq in obsv_freqs.items():
        p_obsv_out[-1][aa] = p_observed(aa, freq)
        entropy_out[-1][aa] = rel_entropy(aa, freq)
        

    

# Printing the calulated output
print('\nPart 2: Probabilities of Observation')        
pp.pprint(p_obsv_out)
print('\nPart 3: Relative Entropy for Observation')        
pp.pprint(entropy_out)
