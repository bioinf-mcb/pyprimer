import os
from fuzzysearch import find_near_matches
from fuzzywuzzy import fuzz
import numpy as np
import warnings
import tables
from numba import njit

class TOOLS:
    @staticmethod
    def match_fuzzily(pattern,
                      sequence,
                      deletions,
                      insertions=0,
                      substitutions=2):
        if pattern in sequence:
            start = sequence.index(pattern)
            return (start, pattern)
        else:
            result = find_near_matches(pattern.encode(),
                                       sequence.encode(),
                                       max_substitutions=substitutions,
                                       max_insertions=insertions,
                                       max_deletions=deletions)
            if len(result) >= 1:
                return result
            else:
                return None
            
    @staticmethod
    def calculate_PPC(F_primer, F_match, R_primer, R_match):
        f_ratio_f = fuzz.ratio(F_primer, F_match)
        f_ratio_r = fuzz.ratio(R_primer, R_match)
        return TOOLS._calculate_PPC(F_primer, F_match, R_primer, R_match, f_ratio_f, f_ratio_r)
    
    @staticmethod
    @njit(cache=True)
    def _calculate_PPC(F_primer, F_match, R_primer, R_match, f_ratio_f, f_ratio_r):
        Fl = float(len(F_primer))
        Fm = np.round((f_ratio_f / 100) * Fl)
        Rl = float(len(R_primer))
        Rm = np.round((f_ratio_r / 100) * Rl)
        f_r_array = np.array([Fm, Rm])
        sigma_m = np.std(f_r_array)
        mi_m = np.mean(f_r_array)
        if mi_m == 0:
            PPC = 0
            return PPC
        CV_m = sigma_m / mi_m
        PPC = (Fm/Fl) * (Rm/Rl) * (1-CV_m)
        if PPC == np.nan:
            PPC = 0
        return PPC