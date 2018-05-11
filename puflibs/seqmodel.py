# Import modules
import scipy
import numpy as np
import pandas as pd
import sys
import os
import ipdb
import itertools

kb = 0.0019872041 # kcal/(mol deg K)
kelvin = 273.15 # Kelvin at 0 deg C

def load_params(model_param_basename='annotations/RNAmap/qMotif_20180302_'):
    """Load the model parameters based on the basename"""
    base_params     = pd.read_csv(model_param_basename + 'term1.csv', index_col=0) #.stack().values
    flip_params     = pd.read_csv(model_param_basename + 'term2_single.csv', index_col=0) #.stack().values
    dflip_params    = pd.read_csv(model_param_basename + 'term2_double.csv', index_col=0, squeeze=True) #.stack().values, 10])
    coupling_params = pd.read_csv(model_param_basename + 'term3.csv', index_col=0, squeeze=True) #.stack().values
    
    return flip_params, base_params, coupling_params, dflip_params
        

def get_ddG_conversion(temperature):
    return -(temperature+kelvin)*kb


def perfect_match_ddG_coupling_terms(sequence, base_penalties_exp, coupling_terms_exp, first_coupling, second_coupling):
    # Inputs:
    # first_coupling--set to True if first (7G) coupling should be included if the conditions are met (this is included so that this can be set to false
    # when flips occur in the coupling region, which will likely prevent coupling)
    # second_coupling--set to True if second (7C) coupling should be included if the necessary conditions are met
    # base_penalties--base penalties in partition function space (exp(-ddG_base/kT))
    # coupling_terms--coupling penalties in partition function space (exp(-ddG_base/kT))
    # sequence register that the mutational penalty contribution is being computed for
    # Output: ddG transformed into partition function space for passed sequence

    # Initialize to a penalty of 0 kcal: 10^(0) = 1; initialize to one in partition function space
    ddG = 1
    # Iterate through base penalties
    for i, base in enumerate(sequence):
        if i==8 and sequence[7]!='A':
            # exception at position 8--nothing happens
            continue
        ddG = ddG*base_penalties_exp.loc[i, base]

    # Apply coupling corrections
    if first_coupling:
        if ((sequence[4] == 'U' or sequence[4] == 'C') and
            sequence[5] == 'A' and
            sequence[6] == 'G' and
            sequence[7] != 'A'):
            ddG = ddG*coupling_terms_exp.loc['c1']
    if second_coupling:
        if sequence[6] == 'C' and sequence[7] != 'A':
            ddG = ddG*coupling_terms_exp.loc['c2']

    
    return ddG

def compute_ensemble_ddG_set(single_dG_values, temperature):
    """Same as below but better starting with an array. Also assumes inputs are in dG,
    not in 'partition function space'"""
    ddG_conversion_factor = get_ddG_conversion(temperature)
    return ddG_conversion_factor*np.log(np.exp(single_dG_values/ddG_conversion_factor).sum(axis=1))
    


def compute_ensemble_ddG(single_dG_values, temperature, needs_exponentiating=False):
    # sums the individual contributions ot the partition function to get the compute partition 
    # function and then converts that into a ddG for the ensemble
    # Inputs:
    # single_dG_values--a list of the single contributions to the partition function from all possible registers
    # Outputs:
    # final_ddG--final ddG of the ensemble
    ddG_conversion_factor = get_ddG_conversion(temperature)

    if needs_exponentiating:
        single_dG_values = np.exp(single_dG_values/ddG_conversion_factor).copy()

    # Sum the logged ddG values to compute the partition function        
    partition_function = np.sum(single_dG_values)

    # Convert the partition function to the ensemble free energy
    
    final_ddG = ddG_conversion_factor*np.log(partition_function)

    return final_ddG

def get_coupling_bool_term1(flip_pos):
    # oonly apply the first coupling term if there is no flip at position 4, 5
    if flip_pos==4 or flip_pos==5:
        return False
    else:
        return True
    
def get_coupling_bool_term2(flip_pos):
    # oonly apply the second coupling term if there is no flip at position 5
    if flip_pos==5:
        return False
    else:
        return True

def get_length9_registers(sequence, base_penalties_exp, coupling_params_exp):
    """for a sequence string, find the set of energies of each 1 nt register."""
    
    seq_length = 9
    registers = {}
    for i in range(len(sequence)-seq_length+1):
        ddG = perfect_match_ddG_coupling_terms(sequence[i:i+seq_length], base_penalties_exp, coupling_params_exp, True, True)
        registers['noflip_%d'%i] = ddG
    return pd.Series(registers)


    

def additive_PUF_flip_model(passed_sequence, flip_params, base_penalties, coupling_params, double_flip_params, temperature, return_ensemble=False):
    # Inputs
    # passed sequence--sequence to compute the affinity for
    # flip_params--list of single flip param penalties (listed as 3/4A, 3/4C, 3/4G, 3/4U, 4/5A, ...)
    # base_penalties--list of single mutation penalties (listed as 1A, 1C, 1G, 1U, 2A, ...)
    # double_flip_params--list of double flip penalties (listed as 3 double flip, 4 double flip, ...)
    # coupling_params--list of coupling adjustments(listed as 7G adjustment followed by 7C adjustment)
    # Outputs
    # ddG predicted for the pased sequence in kcal/mol

    # Compute conversion factor at given Temperature (-T*k*(factor to convert exp to base 10))
    ddG_conversion_factor = get_ddG_conversion(temperature)

    # convert penalties from ddG space to "partition function" space (should be ordered A, C, G, T)
    flip_params_exp = np.exp(flip_params/ddG_conversion_factor)
    double_flip_params_exp = np.exp(double_flip_params/ddG_conversion_factor)
    base_penalties_exp = np.exp(base_penalties/ddG_conversion_factor)
    coupling_params_exp = np.exp(coupling_params/ddG_conversion_factor)

    possible_flip_positions = flip_params_exp.index.tolist()
    possible_double_flip_pos = double_flip_params_exp.index.tolist()
    
    # Convert sequence to a list for easy indexing
    sequence = list(passed_sequence)

    # initialize a list to store affinity for each possible register (note that this will store not true ddG values, but ddG values in "partition function" space)
    single_ddG_values = []
    registers = []
    # compute the non flip ddG values
    for i in range(len(sequence)-8):
        single_ddG_values.append(perfect_match_ddG_coupling_terms(sequence[i:i+9], base_penalties_exp, coupling_params_exp, True, True))
        registers.append('noflip_%d'%i)
    
    # compute the 1 nt flip ddG values
    for i in range(len(sequence)-9):
        current_sequence = sequence[i:i+10]
        for flip_pos in possible_flip_positions:
            seq_not_flipped = current_sequence[:flip_pos]+current_sequence[flip_pos+1:]
            flip_base = current_sequence[flip_pos]

            dG = (flip_params_exp.loc[flip_pos, flip_base]* # this is the penalty of flipping the residue
                  perfect_match_ddG_coupling_terms(seq_not_flipped, base_penalties_exp, coupling_params_exp,
                                                   get_coupling_bool_term1(flip_pos),
                                                   get_coupling_bool_term2(flip_pos)))
            single_ddG_values.append(dG)
            registers.append('flip_%d;pos_%d'%(i, flip_pos))
    
    # double flips
    for i in range(len(sequence)-10):
        current_sequence = sequence[i:i+11]
        
        # 2x1nt flips
        for flip_pos1, flip_pos2 in itertools.combinations(possible_flip_positions, 2):
            seq_not_flipped = current_sequence[:flip_pos1]+current_sequence[flip_pos1+1:flip_pos2] + current_sequence[flip_pos2+1:]
            flip_base1 = current_sequence[flip_pos1]
            flip_base2 = current_sequence[flip_pos2]

            dG = (flip_params_exp.loc[flip_pos1, flip_base1]* # this is the penalty of flipping the residue 1
                  flip_params_exp.loc[flip_pos2, flip_base2]* # this is the penalty of flipping the residue 2
                  perfect_match_ddG_coupling_terms(seq_not_flipped, base_penalties_exp, coupling_params_exp,
                                                   get_coupling_bool_term1(flip_pos1) and get_coupling_bool_term1(flip_pos2),
                                                   get_coupling_bool_term2(flip_pos1) and get_coupling_bool_term2(flip_pos2)))
            single_ddG_values.append(dG)
            registers.append('doubleflip_%d;pos_%d;pos_%d'%(i, flip_pos1, flip_pos2))
        
        # 1x2nt flips
        for flip_pos in possible_double_flip_pos:

            seq_not_flipped = current_sequence[:flip_pos]+current_sequence[flip_pos+2:]
            dG = (double_flip_params_exp.loc[flip_pos]*
                  perfect_match_ddG_coupling_terms(seq_not_flipped, base_penalties_exp, coupling_params_exp,
                                                   get_coupling_bool_term1(flip_pos),
                                                   get_coupling_bool_term2(flip_pos)))
            single_ddG_values.append(dG)
            registers.append('doubleflip_%d;pos_%d'%(i, flip_pos))

    ddG = compute_ensemble_ddG(single_ddG_values, temperature)
    
    if return_ensemble:
        return ddG, pd.Series(ddG_conversion_factor*np.log(single_ddG_values), index=registers)
    else:
        return ddG

def interpret_col_names(col_names):
    """ given the 'register' annotations above, group them more meaningfully"""
    annotations = pd.Series(1, index=col_names)
    annotations.loc[[idx for idx in col_names if idx.find('noflip')==0 and idx!='noflip_0']] = 2
    annotations.loc[[idx for idx in col_names if idx.find('flip')==0]] = 4
    annotations.loc[[idx for idx in col_names if idx.find('doubleflip')==0]] = 8

    return annotations

def flag_ensemble(ddG_ensemble_vec, cutoff=1):
    """Given a set of measurements, find things within cutoff and generate flag with annotations"""
    close_enough_vec = [idx for idx, val in ddG_ensemble_vec.iteritems() if val < cutoff]
    annotations = interpret_col_names(ddG_ensemble_vec.index.tolist())
    flag =  annotations.loc[close_enough_vec].unique().sum()
    return flag
    
def determine_seq_occupancy(seq, temperature):
    """Break a seq up into 11 nt chunks and return occupancy relative to consensus site for each seq."""
    
