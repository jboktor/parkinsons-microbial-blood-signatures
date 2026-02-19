# stitchr_wrapper.py

from Stitchr import stitchrfunctions as fxn
from Stitchr import stitchr as st

def run_stitchr(v, j, c, cdr3, chain):
    # Initialise the necessary data
    species='HUMAN'
    tcr_dat, functionality, partial = fxn.get_imgt_data(chain, st.gene_types, species)
    codons = fxn.get_optimal_codons('', species)
    
    # Define the TCR bits dictionary
    tcr_bits = {
        'v': v,
        'j': j,
        'cdr3': cdr3,
        'l': v,
        'c': c,
        'skip_c_checks': False,
        'species': species,
        'seamless': False,
        '5_prime_seq': '', 
        '3_prime_seq': '', 
        'name': 'TCR'
    }
    
    
    # Run stitchr on the rearrangement
    stitched = st.stitch(tcr_bits, tcr_dat, functionality, partial, codons, 3, '')
    
    # Return only the amino acid sequence
    return stitched if stitched[-1] == 0 else None
