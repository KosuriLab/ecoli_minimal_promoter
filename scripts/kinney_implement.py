import pandas as pd

filename = '../ref/kinney_RNAP_matrix.txt'

# list of dictionaries, index 0 corresponds to position -41 and last list index 
# corresponds to position -1 each entry is a dictionary where keys are nucleotides 
# and values are the probabilities for that nucleotide at that position
energy_matrix = []
# first row is position -41, last row is position -1
with open(filename) as infile:
    nts = ['A', 'C', 'G', 'T']
    for line in infile:
        if not line.startswith('#'):
            probs = map(float, line.split())
            energy_matrix.append(dict(zip(nts, probs)))  

# read in library information and expression
data = pd.read_table('rlp5Min_SplitVariants.txt', sep=' ')

# grab -41:-1 equivalent sequence for each variant
# there are 29 bp of initial transcribed region which is at the end of the sequence, 
# preceded by the TSS
# all sequences are 150bp
data['rnap_site'] = [variant[-70:-30] for variant in data.variant]
data['rnap_site_no_spacer'] = [site[:11] + site[-12:] for site in data.rnap_site]
energy_matrix_no_spacer = energy_matrix[:11] + energy_matrix[-12:]

def score_site(mat, site):
    # first letter in string corresponds to position -41 of promoter, first index of energy matrix list
    score_total = 0.
    for i in range(len(site)):
        nt = site[i]
        score_total += mat[i][nt]
    return score_total

data['rnap_site_score'] = [score_site(energy_matrix, site) for site in data.rnap_site]
data['rnap_site_no_spacer_score'] = [score_site(energy_matrix_no_spacer, site) for site in data.rnap_site_no_spacer]

# save to text so we can do the modeling in R
data.to_csv('../processed_data/min_rnap_scores.txt', sep='\t', index=False)