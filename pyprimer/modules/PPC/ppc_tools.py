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
                      deletions=0,
                      insertions=0,
                      substitutions=2):
        """Match strings fuzzily usign Levenshtein distance

        Arguments:
            pattern {str} -- Searched pattern
            sequence {str} -- The sequence queried for the pattern

        Keyword Arguments:
            deletions {int} -- How many deletions can differentiate the patterns (default: {0})
            insertions {int} -- How many deletions can insertions the patterns (default: {0})
            substitutions {int} -- How many deletions can substitutions the patterns (default: {2})

        Returns:
            tuple((int, str)) -- Position and value of the closest matching pattern. 
            (None, "") if the pattern was not found.
        """
        if pattern in sequence:
            start = sequence.index(pattern)
            return (start, pattern)
        else:
            result = find_near_matches(pattern,
                                       sequence,
                                       max_substitutions=substitutions,
                                       max_insertions=insertions,
                                       max_deletions=deletions)
            return result
            # if len(result) > 1:
            #     dist_list = np.ndarray(shape = (len(result),))
            #     for i in range(len(result)):
            #         dist_list[i] = result[i].dist
            #     idx_max = np.where(dist_list == np.max(dist_list))


            # if len(result) != 0:
            #     idx_max = np.argmax([i.dist for i in result])
            #     return result[idx_max].start, result[idx_max].matched
            # else:
            #     return None, ""

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
    

class _MemSaver:
    def __init__(self, tempdir, fname, col_list):
        self.tempdir = tempdir
        self.fname = fname
        self.col_list = col_list
        os.makedirs(self.tempdir, exist_ok=True)
    
    def save_group(self, group_df, key):
        # extra weird fix for this problem:
        # https://stackoverflow.com/questions/60677863/python-pandas-append-dataframe-with-array-content-to-hdf-file
        # https://mlog.club/article/5501050
        for i in range(len(self.col_list)-1):
            group_df[self.col_list[i]
                        ] = group_df[self.col_list[i]].astype(str)

        group_df[self.col_list[len(
            self.col_list)-1]] = group_df[self.col_list[len(self.col_list)-1]].astype(float)
        warnings.filterwarnings(
            'ignore', category=tables.NaturalNameWarning)

        group_df.to_hdf(path_or_buf=os.path.join(self.tempdir, self.fname),
                        key=key,
                        mode="a",
                        format="table",
                        data_columns=True)
        del group_df
