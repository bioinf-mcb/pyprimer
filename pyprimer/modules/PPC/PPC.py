# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os
from numba import jit, prange
from tqdm import trange
from tqdm import tqdm
import h5py

import operator
import dask.dataframe as dd
import dask.multiprocessing
import dask.threaded
import sys
import time
import warnings
from pyprimer.modules.PPC.ppc_tools import TOOLS, _MemSaver

from pyprimer.utils.sequence import Essentials

class PPC(object):
    COL_LIST = ["F Primer Name",
                "F Primer Version",
                "R Primer Name",
                "R Primer Version",
                "Sequence Header",
                "Amplicon Sense",
                "Amplicon Sense Length",
                "Amplicon Sense Start",
                "Amplicon Sense End",
                "PPC"]

    SUMMARY_COL_LIST = ["Primer Group",
                        "F Version",
                        "R Version",
                        "Mean PPC",
                        "Sequences matched(%)"]

    def __init__(self, primer_df, 
                        sequence_df,
                        memsave=False,
                        tempdir="./tmp/",
                        fname="PCRBenchmark.h5"):
        self.primers = primer_df
        self.sequences = sequence_df
        self.memsave = memsave
        if memsave:
            self._saver = _MemSaver(tempdir, fname, self.COL_LIST)

    def analyse_primers(self,
                        deletions=0,
                        insertions=0,
                        substitutions=2,
                        nCores=2) -> pd.DataFrame:
        """Analyse effectiveness of the PCR primers against sequences.

        Args:
            deletions (int, optional): Number of deletions allowed in the target sequence. Defaults to 0.
            insertions (int, optional): Number of insertions allowed in the target sequence. Defaults to 0.
            substitutions (int, optional): Number of substitutions allowed in the target sequence. Defaults to 2.
            nCores (int, optional): Number of cores for concurrent processing. Defaults to 2.

        Returns:
            pd.DataFrame: Info about the group of primers with structure specified by SUMMARY_COL_LIST const.
        """
        unique_groups = self.primers["ID"].unique()

        bench_df = pd.DataFrame(columns=self.COL_LIST)
        summary = pd.DataFrame(columns=self.SUMMARY_COL_LIST)

        filter_forward = self.primers["Type"] == "F"
        filter_reverse = self.primers["Type"] == "R"
        filter_probe = self.primers["Type"] == "R"

        with tqdm(unique_groups) as pbar:
            for group in pbar:
                filter_group = self.primers["ID"] == group
                Fs = self.primers.loc[filter_group & filter_forward].values
                Rs = self.primers.loc[filter_group & filter_reverse].values
                Ps = self.primers.loc[filter_group & filter_probe].values

                f_versions = Fs[:, 5]
                f_names = Fs[:, 2]
                r_versions = Rs[:, 5]
                r_names = Rs[:, 2]
                p_versions = Ps[:, 5]
                p_names = Ps[:, 2]

                group_stats = self._calculate_group_summary(self.sequences, group, 
                                    f_versions, r_versions, p_versions,
                                    f_names, r_names, p_names,
                                    deletions, insertions, substitutions, nCores)
                summary = summary.append(group_stats)
        
        if self.memsave:
            print("Extended benchmark results were written to {}".format(
                os.path.join(self._saver.tempdir, "PCRBenchmark.h5")))
        
        return summary

    def get_primer_metrics(self, F, R, P=None, deletions=0, insertions=0, substitutions=2, nCores=2):
        """Calculate Design Points for a primer or group of primers

        Args:
            F (str or np.array): Forward primer(s)
            R (str or np.array): Reverse primer(s)
            P (str or np.array, optional): Probe(s)
            deletions (int): Number of deletions allowed by the fuzzy search. Defaults to 0.
            insertions (int): Number of insertions allowed by the fuzzy search. Defaults to 0.
            substitutions (int): Number of substitutions allowed by the fuzzy search. Defaults to 2.
            nCores (int, optional): Number of cores for concurrent processing. Defaults to 2.
        """
        GC_GOAL = 60
        NORM_FACTOR = GC_GOAL**2/50

        if type(F)==str:
            F = np.array([F])
        if type(R)==str:
            R = np.array([R])
        if P is not None and type(P)==str:
            P = np.array([P])

        if type(F)==list:
            F = np.array(F)
        if type(R)==list:
            R = np.array(R)
        if P is not None and type(P)==list:
            P = np.array(P) 

        if type(F)==np.ndarray and F.ndim==2:
            F_new = []
            for i in F:
                F_new.append("".join(i))
            F = np.asarray(F_new)
        if type(R)==np.ndarray and R.ndim==2:
            R_new = []
            for i in R:
                R_new.append("".join(i))
            R = np.asarray(R_new)
        
        stats = self._calculate_group_summary(self.sequences, "DP", F, R, P, deletions=0, insertions=0, substitutions=2, nCores=nCores, permute=False)

        gc = []
        mean_ppc = []
        matched = []
        for f, r in zip(F, R):
            filt = (stats["F Version"]==f) & (stats["R Version"]==r)
            gc.append(
                (Essentials.GCcontent(f) + Essentials.GCcontent(r) - GC_GOAL*2)**2 
                / NORM_FACTOR
                )
            mean_ppc.append((stats[filt]["Mean PPC"].values*100)[0])
            matched.append(stats[filt]["Sequences matched(%)"].values[0])

        gc = np.asarray(gc).round(3)
        mean_ppc = np.asarray(mean_ppc).round(3)
        matched = np.asarray(matched).round(3)

        return mean_ppc + matched - gc

    def _calculate_group_summary(self, sequences, group_name,         
        f_versions, r_versions, p_versions,
        f_names=None, r_names=None, p_names=None, 
        deletions=0, insertions=0, substitutions=2, nCores=2,
        permute = True):
        
        def __calculate(x):
            # print(x, 
            #     f_versions, r_versions, p_versions, 
            #     f_names, r_names, p_names,
            #     deletions, insertions, substitutions, permute)
            return self._calculate_stats(x, 
                f_versions, r_versions, p_versions, 
                f_names, r_names, p_names,
                deletions, insertions, substitutions, permute)

        if f_names is None:
            f_names = np.array([""]*len(f_versions))
        if r_names is None:
            r_names = np.array([""]*len(r_versions))
        if p_names is None:
            p_names = np.array([""]*len(p_versions))

        dsequences = dd.from_pandas(self.sequences, npartitions=nCores)
        # self.bar = tqdm(total=len(self.sequences))
        df_series = dsequences.map_partitions(
            lambda df: df.apply(__calculate, axis=1), meta=('df', None)
            ).compute(scheduler='processes')
        
        print("Concatenating")
        group_df = pd.concat(df_series.tolist())
        print("Concatenated")
        if self.memsave:
            self._saver.save_group(group_df, group_name)

        print("Saved (or not)")

        v_stats = self._craft_summary(group_df, group_name)
        print("Crafted summary")

        group_stats = pd.DataFrame(v_stats, columns=self.SUMMARY_COL_LIST)
        return group_stats

    def _calculate_stats(self, sequences,
        f_vers, r_vers, p_vers,
        f_names, r_names, p_names,
        deletions, insertions, substitutions, permute) -> pd.DataFrame:
        """Calculate statistics for sequence set for every possible primer version combination

        Args:
            sequences (dask.Series): Sequences that we want to test the primers against
            Fs (np.array): All forward primers from the group
            Rs (np.array): All reverse primers from the group
            Ps (np.array): All probes from the group
            deletions (int): Number of deletions allowed by the fuzzy search
            insertions (int): Number of insertions allowed by the fuzzy search
            substitutions (int): Number of substitutions allowed by the fuzzy search

        Returns:
            pd.DataFrame: Info about the group of primers with structure specified by COL_LIST const.
        """
        res = []
        header = sequences[0]
        # self.bar.update()
        if permute:
            for f_ver, f_name in zip(f_vers, f_names):
                start, f_match = TOOLS.match_fuzzily(
                    f_ver, sequences[1], deletions, insertions, substitutions)

                for r_ver, r_name in zip(r_vers, r_names):
                    if start is not None:
                        r_start, r_match = TOOLS.match_fuzzily(
                            r_ver, sequences[2], deletions, insertions, substitutions)
                        try:
                            end = len(sequences[1]) - 1 - r_start
                        except TypeError:
                            end = None
                    
                    if start is None or end is None or end<0:
                        amplicon = ""
                        amplicon_length = 0
                        start = None
                        end = None
                        PPC = 0
                    else:
                        amplicon = sequences[1][start:end]
                        amplicon_length = len(amplicon)

                        PPC = TOOLS.calculate_PPC(F_primer=f_ver,
                                                    F_match=f_match,
                                                    R_primer=r_ver,
                                                    R_match=r_match)

                    res.append([f_name, f_ver, r_name, r_ver,
                                        header, amplicon, amplicon_length, start, end, PPC])

        else:
            for primer_set in range(len(f_vers)):
                f_ver = f_vers[primer_set]
                f_name = f_names[primer_set]
                r_ver = r_vers[primer_set]
                r_name = r_names[primer_set]

                start, f_match = TOOLS.match_fuzzily(
                    f_ver, sequences[1], deletions, insertions, substitutions)

                if start is not None:
                    r_start, r_match = TOOLS.match_fuzzily(
                            r_ver, sequences[2], deletions, insertions, substitutions)
                    try:
                        end = len(sequences[1]) - 1 - r_start
                    except TypeError:
                        end = None

                if start is None or end is None or end<0:
                    amplicon = ""
                    amplicon_length = 0
                    start = None
                    end = None
                    PPC = 0
                else:
                    amplicon = sequences[1][start:end]
                    amplicon_length = len(amplicon)

                    PPC = TOOLS.calculate_PPC(F_primer=f_ver,
                                                F_match=f_match,
                                                R_primer=r_ver,
                                                R_match=r_match)

                res.append([f_name, f_ver, r_name, r_ver,
                                    header, amplicon, amplicon_length, start, end, PPC])
        df = pd.DataFrame(res, columns=self.COL_LIST)
        return df

    def _craft_summary(self, group_df: pd.DataFrame, group: str) -> dict:
        """Create summary of the analysis of a group of primers.

        Args:
            group_df (pd.DataFrame): Dataframe with detailed analysis
            group (str): Name of the group of primers

        Returns:
            dict: Info about the group of primers specified 
                  by SUMMARY_COL_LIST const.
        """
        v_stats = {key: [] for key in self.SUMMARY_COL_LIST}
        for fversion in group_df["F Primer Version"].unique():
            for rversion in group_df["R Primer Version"].unique():
                filter_r_version = (group_df["R Primer Version"] == rversion)
                filter_f_version = (group_df["F Primer Version"] == fversion)
                filter_matching = filter_r_version & filter_f_version

                if np.sum(filter_matching)>0:
                    n_seqs = np.sum(filter_matching)
                    seqs_matched = np.sum(filter_matching & (group_df["Amplicon Sense Length"] != 0))

                    mean_ppc = group_df.loc[filter_matching, "PPC"].mean().round(5)
                    
                    v_stats["Primer Group"].append(group)
                    v_stats["F Version"].append(fversion)
                    v_stats["R Version"].append(rversion)
                    v_stats["Mean PPC"].append(mean_ppc)
                    v_stats["Sequences matched(%)"].append(
                        (seqs_matched / n_seqs)*100)
        
        return v_stats