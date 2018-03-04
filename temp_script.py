#!/Library/Frameworks/Python.framework/Versions/3.6/bin/python3

## Loading in Data ##
# protein sequences
prots = [['G', 'I', 'G', 'V', 'D', 'D', 'V', 'I', 'K', 'L'],
         ['Q', 'L', 'G', 'P', 'K', 'D', 'R', 'L', 'S', 'I'],
         ['G', 'L', 'G', 'V', 'F', 'D', 'A', 'I', 'S', 'L'],
         ['Q', 'L', 'G', 'I', 'N', 'D', 'Q', 'I', 'S', 'M'],
         ['G', 'I', 'G', 'K', 'Q', 'N', 'E', 'V', 'S', 'L']]

num_seqs = # FIXME

# amino acid backgroud frequencies
aa_freq_dict = {'A' : 0.073, 'C' : 0.025, 'D' : 0.050, 'E' : 0.061,
                'F' : 0.042, 'G' : 0.072, 'H' : 0.023, 'I' : 0.053,
	        'K' : 0.064, 'L' : 0.089, 'M' : 0.023, 'N' : 0.043,
	        'P' : 0.052, 'Q' : 0.040, 'R' : 0.052, 'S' : 0.073,
	        'T' : 0.056, 'V' : 0.063, 'W' : 0.013, 'Y' : 0.033}

# list of amino acids
aa_list = aa_freq_dict.keys()

## Problem 1 ##
# positions to analyze (zero indexed)
positions = [0, 4, 9]

# list of dictionaries to keep track of counts
counts_out = []

# looping through queried positions
for i, iPos in enumerate(positions):
    # appending a new dict for each position
    counts_out.append({})

    # looping through each protein sequence
    for seq in prots:
        # extracting the amino acid we care about
        aa = seq[iPos]

        # increase or set the index to one
        current_count = counts_out[i].get(aa, 0)
        counts_out[i][aa] = current_count + 1
            

    

