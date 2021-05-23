# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os
from numba import jit, prange
from tqdm import trange
from tqdm import tqdm
import h5py
import gc
import operator
import dask.dataframe as dd
from dask.distributed import Client, progress, LocalCluster
import dask.multiprocessing
import dask.threaded
import sys
import time
import warnings
from pyprimer.modules.Benchmark.benchmark_tools import TOOLS
from pyprimer.utils.essentials import Essentials


class Benchmark(object):

    BENCHMARK_qPCR_COL_LIST = [
                "F Primer Name",
                "F Primer Version",
                "P Probe Name",
                "P Probe Version",
                "R Primer Name",
                "R Primer Version",
                "Sequence Header",
                "Amplicon Sense",
                "Amplicon Sense Length",
                "Amplicon Sense Start",
                "Amplicon Sense End",
                "PPC"]

    SUMMARY_qPCR_COL_LIST = ["Primer Group",
                        "F Version",
                        "P Version",
                        "R Version",
                        "Mean PPC",
                        "Sequences matched(%)"]

    def __init__(self, primer_df, sequence_df, nCores = 4):
        self.primers = primer_df
        # self.sequences = dd.from_pandas(sequence_df, npartitions = nCores)
        self.sequences = sequence_df
        self.nCores = nCores

    def qPCR_performance(self, deletions = 0, insertions = 0, substitutions = 0,
                         savedir = ".", hdf_fname = 'pyprimer_benchmark.h5', csv_fname = "pyprimer_summary.csv",):
        def generate_group_summary(group_df, group, col_list):
            v_stats = dict((key,[]) for key in col_list)
            for fversion in group_df["F Primer Version"].unique():
                for rversion in group_df["R Primer Version"].unique():
                    for pversion in group_df["P Probe Version"].unique():
                        mean_ppc = group_df.loc[(group_df["F Primer Version"] == fversion) & (group_df["R Primer Version"] == rversion) & (group_df["P Probe Version"] == pversion), "PPC"].mean()
                        seqs_matched = len(group_df.loc[(group_df["F Primer Version"] == fversion) & (group_df["R Primer Version"] == rversion) & (group_df["P Probe Version"] == pversion) & (group_df["Amplicon Sense Length"] != 0), "Amplicon Sense Length"])
                        n_seqs = group_df.loc[(group_df["F Primer Version"] == fversion) & (group_df["R Primer Version"] == rversion) & (group_df["P Probe Version"] == pversion), :].shape[0]
                        v_stats["Primer Group"].append(group)
                        v_stats["F Version"].append(fversion)
                        v_stats["P Version"].append(pversion)
                        v_stats["R Version"].append(rversion)
                        v_stats["Mean PPC"].append(mean_ppc)
                        try:
                            percent_matched = (seqs_matched / n_seqs)*100
                        except:
                            percent_matched = 0
                        v_stats["Sequences matched(%)"].append(percent_matched)
            group_stats = pd.DataFrame(v_stats, columns = col_list)
            return group_stats

        def help_analyse(x):
            return analyse(x, Fs, Rs, Ps, self.BENCHMARK_qPCR_COL_LIST,
                           self.deletions, self.insertions, self.substitutions)

        def analyse(sequences_df, Fs, Rs, Ps, col_list, deletions, insertions, substitutions):
            res = []
            for i in range(sequences_df.shape[0]):
                sequences = sequences_df.iloc[i,:].values.tolist()
                for f in Fs:
                    for r in Rs:
                        header = sequences[0]
                        f_name = f[2]
                        r_name = r[2]
                        f_ver = f[5]
                        r_ver = r[5]
                        f_res = TOOLS.match_fuzzily(f_ver, sequences[1], deletions, insertions, substitutions)
                        r_res = TOOLS.match_fuzzily(r_ver, sequences[2], deletions, insertions, substitutions)

                        # if (type(f_res) == type(tuple())) and (type(r_res) == type(tuple())):
                        #     start = f_res[0]
                        #     r_start = r_res[0]
                        #     end = (len(sequences[1]) - 1) - r_start
                        #     f_match = f_ver
                        #     r_match = r_ver
                        #     if start < end:
                        #         amplicon = sequences[1][start:end]
                        #         amplicon_length = len(amplicon)
                        #     else:
                        #         amplicon = ""
                        #         amplicon_length = 0
                        
                        if (f_res == None) or (r_res == None):
                            start = None
                            end = None
                            amplicon = ""
                            amplicon_length = 0
                            f_match = ""
                            r_match = ""
                            p_ver = ""
                            p_match = None
                            PPC = 0
                        
                        else:
                            Forwards = {}
                            if type(f_res) == type(tuple()):
                                Forwards[0] = (f_res[0], f_ver, 0) # (start, match, distance)
                            else:
                                for f_i in range(f_res):
                                    Forwards[f_i] = (f_res[f_i].start, f_res[f_i].matched, f_res[f_i].dist)
                            Reverses = {}
                            if type(r_res) == type(tuple()):
                                Reverses[0] = (r_res[0], r_ver, 0)
                            else:
                                for r_i in range(r_res):
                                    Reverses[r_i] = (r_res[r_i].start, r_res[f_i].matched, r_res[f_i].dist)
                            matches = {}
                            for k_f, v_f in Forwards.items():
                                start = v_f[0]
                                for k_r, v_r in Reverses.items():
                                    end = (len(sequences[1]) - 1) - v_r[0]
                                    if end < start:
                                        matches[f"{k_f}:{k_r}"] = (False, "")
                                    amplicon = sequences[1][start:end]
                                    for p in Ps:
                                        p_ver = p[5]
                                        p_name = p[2]
                                        p_res = TOOLS.match_fuzzily(p_ver, sequences[1], deletions, insertions, substitutions)
                                        if p_res == None:
                                            matches[f"{k_f}:{k_r}"] = (False, p_ver)
                                        else:
                                            matches[f"{k_f}:{k_r}"] = (True, p_ver)
                            target_dist = np.Inf
                            n_match = 0
                            for k, v in matches.items():
                                if v[0] == True:
                                    n_match += 1
                                    klist = k.split(":")
                                    k_f = int(klist[0])
                                    k_r = int(klist[1])
                                    f_good = Forwards[k_f]
                                    r_good = Reverses[k_r]
                                    mean_dist = (f_good[2] + r_good[2] + 1e-6)/2 # 1e-6 for smoothing
                                    if mean_dist < target_dist:
                                        target_dist = mean_dist
                                        start = f_good[0]
                                        f_match = f_good[1]
                                        end = (len(sequences[1]) - 1) - r_good[0]
                                        r_match = r_good[1]
                                        amplicon = sequences[1][start:end]
                                        amplicon_length = len(amplicon)
                                        p_ver = v[1]
                                        p_match = v[1]
                            if n_match == 0:
                                start = None
                                end = None
                                amplicon = ""
                                amplicon_length = 0
                                f_match = ""
                                r_match = ""
                                p_ver = ""
                                p_match = None
                                PPC = 0
                            else:
                                PPC = TOOLS.calculate_PPC(F_primer=f_ver,
                                                        F_match=f_match,
                                                        R_primer=r_ver,
                                                        R_match=r_match)

                        res.append([f_name, f_ver, f_name, p_ver,
                                    r_name, r_ver, header, amplicon,
                                    amplicon_length, start, end, PPC])

            df = pd.DataFrame(res, columns = col_list)
            return df

        self.savedir = savedir
        self.hdf_fname = hdf_fname
        self.csv_fname = csv_fname
        self.deletions = deletions
        self.insertions = insertions
        self.substitutions = substitutions

        unique_groups = self.primers["ID"].unique()
        # bench_df = pd.DataFrame(columns = self.BENCHMARK_qPCR_COL_LIST)
        # bench = dd.from_pandas(bench_df, npartitions = self.nCores)
        summary = pd.DataFrame(columns = self.SUMMARY_qPCR_COL_LIST)
        os.makedirs(self.savedir, exist_ok = True)
        #############################Experimental method##############################
        df_parts = []
        chunk_size = (self.sequences.shape[0] // self.nCores) + 1
        id_min = 0
        id_max = chunk_size
        for chunk in range(self.nCores):
            if id_max < self.sequences.shape[0]:
                part = self.sequences.iloc[id_min:id_max, :]
            else:
                part = self.sequences.iloc[id_min:, :]
            part.reset_index(drop = True, inplace = True)
            id_min += chunk_size
            id_max += chunk_size
            df_parts.append(part)
        cluster = LocalCluster(n_workers = self.nCores, threads_per_worker = 2)
        client = Client(cluster)# "192.168.45.241:8786"
        ##############################################################################
        for group in tqdm(unique_groups):
            Fs = self.primers.loc[(self.primers["ID"] == group) & (self.primers["Type"] == "F"),:].values
            Rs = self.primers.loc[(self.primers["ID"] == group) & (self.primers["Type"] == "R"),:].values
            Ps = self.primers.loc[(self.primers["ID"] == group) & (self.primers["Type"] == "P"),:].values
            print(f"Processing group {group}\n")
            # df_parts = self.sequences.map_partitions(
            #     lambda df: df.apply(help_analyse, axis = 1), meta=('df', None)
            # ).compute(scheduler = "threads")
            #############################Experimental method##############################
            futures = client.map(help_analyse, df_parts, pure=False)
            progress(futures)
            group_df = pd.DataFrame(columns = self.BENCHMARK_qPCR_COL_LIST)
            for element in client.gather(futures):
                group_df.append(element)
            ##############################################################################
            group_df.reset_index(drop = True, inplace = True)
            print("Performance computed, generating group summary\n")
            group_stats = generate_group_summary(group_df, group, self.SUMMARY_qPCR_COL_LIST)
            summary = summary.append(group_stats)
            print("Summary generated, saving group benchmark to HDF file\n")
            # group_df["Amplicon Sense Length"].apply(lambda x: str(x))
            # group_df["Amplicon Sense Start"].apply(lambda x: str(x))
            # group_df["Amplicon Sense End"].apply(lambda x: str(x))
            # group_df["PPC"].apply(lambda x: str(x))
            group_df = group_df.astype(str)
            gdf = dd.from_pandas(group_df, npartitions = self.nCores)
            gdf.to_hdf(
                path_or_buf = os.path.join(self.savedir, self.hdf_fname),
                key = group,
                mode = "a",
                complevel = 0,
                format = "table",
                errors = "ignore",
                min_itemsize = 850,
                scheduler = "processes",
                data_columns = True
            )
            #############################Experimental method##############################
            client.restart()
            ##############################################################################
        summary.to_csv(os.path.join(self.savedir, self.csv_fname), index = False)
        print(f"Benchmark results saved to {os.path.join(self.savedir, self.hdf_fname)}\n")
        print(f"Benchmark summary saved to {os.path.join(self.savedir, self.csv_fname)}\n")


