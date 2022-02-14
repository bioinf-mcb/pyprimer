import os
from fuzzysearch import find_near_matches
from fuzzywuzzy import fuzz
import numpy as np
import warnings
import tables
from numba import njit
import string
import random

class TOOLS:
    @staticmethod
    def match_fuzzily(pattern,
                      sequence,
                      deletions = None,
                      insertions = None,
                      substitutions = None,
                      distance = None):
        result = find_near_matches(pattern,
                                    sequence,
                                    max_substitutions=substitutions,
                                    max_insertions=insertions,
                                    max_deletions=deletions,
                                    max_l_dist = distance)
        if len(result) >= 1:
            return result
        else:
            return None
    @staticmethod
    def extract_template(sequence, start_, end_):
        end = len(sequence) - start_
        start = len(sequence) - end_
        result = sequence[start:end]
        return result
    
    @staticmethod
    def check_correctness(res_bind, res_stability, gibbs_threshold, temp_c):
        if res_bind.structure_found:
            if res_bind.tm > temp_c:
                if res_bind.dg < gibbs_threshold:
                    if res_stability.structure_found:
                        if res_stability.tm > temp_c:
                            if res_stability.dg < gibbs_threshold:
                                verdict = True
                            else:
                                verdict = False
                        else:
                            verdict = False
                    else:
                        verdict = False
                else:
                    verdict = False
            else:
                verdict = False
        else:
            verdict = False
        return verdict

    @staticmethod
    def id_generator(chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(8))

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