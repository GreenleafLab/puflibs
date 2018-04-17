import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

        
RT=0.592090122
class Single():
    
    def __init__(self, single_base_ddG=None, s=None):
        
        if single_base_ddG is None:
            single_base_ddG = pd.read_table('/lab/sarah/pufProject/K562_encode/annotations/RNAmap/single_base_penalties.dat', index_col=0)
        self.single_base_ddG = single_base_ddG
        self.single_base_Kd = np.exp(single_base_ddG/RT)
        self.pwm = find_pwm(self.single_base_Kd)
    
    
    def get_ddG(self, s=None):
        """Using the single base penalties, calculate the ddG"""
        # chieck input
        s = self._check_sequence(s)
        return self._get_ddG(s, self.single_base_ddG)
    
    def _check_sequence(self, s=None):
        """Need to be sure sequence is good"""
        # chieck input
        if s is None:
            if self.s is None:
                raise KeyError("Need to define input sequence s")
            s = self.s
            
        # check length of s
        if len(s) != 9:
            raise ValueError("input sequence must be length 9 bases")
        
        # check bases
        if not all([c in ['A', 'G', 'C', 'U'] for c in s]):
            if all([c in ['A', 'G', 'C', 'U', 'T'] for c in s]):
                s = s.replace('T', 'U')
                warnings.warn("T's detected: replacing with U's", UserWarning)
            else:
                raise ValueError("input sequence must not have any bases but standard nucleotides")
        return s
    
    
    def _get_ddG(self, s, single_base_ddG):    
        # calculate
        total_ddG = 0.
        for i, c in enumerate(s):
            # additive element to ddG is the single base penalty
            val = single_base_ddG.loc[i, c]
            
            # exception at position 8 if position 7 is G and position 5 != A
            if i == 7:
                if s[6] == 'G' and s[4] != 'A':
                    val = 0
            
            # exception at position 9 if position 7 is G and position 5 != A, or if position 8 does not equal A
            if i == 8:
                if s[7] != 'A' or (s[6] == 'G' and s[4] != 'A'):
                    val = 0
            
            # append to total energy
            total_ddG += val
        return total_ddG
        

            
    def find_log_odds_score(self, s):
        """For a given sequence, determine the log odds score given the PWM."""
        # chieck input
        s = self._check_sequence(s)        
        pwm = self.pwm
        return log_odds_score(s, pwm)

            
def log_odds_score(s, pwm):
    """For a given sequence, determine the log odds score given the PWM."""
    pseudocount = 0.001
    scores = []
    for i, c in enumerate(s):
        scores.append(np.log((pwm.loc[i, c] + pseudocount)/0.25))    
        
    return np.sum(scores)
        
def find_pwm(single_base_Kd):        
    """Using single base penalties, determine a model for the PWM"""
    
    pwm = {}
    for i in range(9):
        rel_kds = single_base_Kd.loc[i]
        rel_occ = 1./rel_kds
        rel_occ_norm = rel_occ/rel_occ.sum()
        pwm[i] = rel_occ_norm
    
    pwm = pd.concat(pwm).unstack()
    return pwm
                
def insert_base_in_pwm(position, pwm, num_insertions=1, rel_kds=pd.Series([1]*4, index=['A', 'C', 'G', 'U'])):
    """insert a base in the pwm at position"""
    rel_occ = 1./rel_kds
    rel_occ_norm = rel_occ/rel_occ.sum()    
    pwm_new = pd.concat([pwm.iloc[:position],
                         pd.DataFrame([rel_occ_norm]*num_insertions),
                         pwm.iloc[position:]], ignore_index=True)
    return pwm_new
        
def save_pwm(pwm, filename, odds_threshold, name=None):
    """Save a PWM in the HOMER format"""
    consensus_seq = ''.join([row.idxmax() for i, row in pwm.iterrows()])
    eps = 0.001
    if name is None:
        name = consensus_seq
    with open(filename, 'w') as f:
        f.write('>%s\n'%('\t'.join([consensus_seq, name, '%.2f'%odds_threshold])))
        for i, row in pwm.iterrows():
            f.write('%s\n'%('\t'.join(['%.3f'%(val+eps) for val in row])))
            
        