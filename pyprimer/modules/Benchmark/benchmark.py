# -*- coding: utf-8 -*-
import pandas as pd
import csv
import numpy as np
import pickle
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
import warnings
import tables
from contextlib import suppress
from distributed.comm.core import CommClosedError
import logging
warnings.simplefilter('ignore', tables.NaturalNameWarning)
warnings.simplefilter('ignore', tables.PerformanceWarning)
warnings.simplefilter("ignore", UserWarning)

class Benchmark(object):

    BENCHMARK_qPCR_COL_LIST = [
                "F Primer Name",
                "F Primer Version",
                "P Probe Version",
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

    def __init__(self, primer_df, sequence_df, savedir = "./results", tmpdir = "./tmp", nCores = 4):
        self.primers = primer_df
        self.nCores = nCores
        self.tmpdir = tmpdir
        self.savedir = savedir
        print("Preparing data chunks")
        try:
            flist = os.listdir(self.tmpdir)
            if len(flist) >= 1:
                chunkpaths = [os.path.join(self.tmpdir,f) for f in flist]
                self.chunkpaths = chunkpaths
                del sequence_df
            else:
                raise ValueError()
        except:
            os.makedirs(self.tmpdir, exist_ok = True)
            chunkpaths = []
            chunk_size = (sequence_df.shape[0] // self.nCores) + 1
            id_min = 0
            id_max = chunk_size
            for chunk in range(self.nCores):
                if id_max < sequence_df.shape[0]:
                    part = sequence_df.iloc[id_min:id_max, :]
                elif id_min > sequence_df.shape[0]:
                    break
                else:
                    part = sequence_df.iloc[id_min:, :]
                id_min += chunk_size
                id_max += chunk_size
                fpath = os.path.join(self.tmpdir,f"chunk_{chunk}.csv")
                chunkpaths.append(fpath)
                part.to_csv(fpath, index = False)
            self.chunkpaths = chunkpaths
        

    def qPCR_performance(self, deletions = 0, insertions = 0, substitutions = 0,
                         fname = 'pyprimer_benchmark.feather', csv_fname = "pyprimer_summary.csv",):
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
            group_stats = pd.DataFrame(v_stats)
            return group_stats

        def analyse(sequences_path, Fs, Rs, Ps, col_list, deletions, insertions, substitutions):
            res = []
            with open(sequences_path, "r", newline='') as csvfile:
                seqreader = csv.reader(csvfile, delimiter = ',', quotechar ='"')
                for sequences in seqreader:
                    if sequences[0] == "Header":
                        pass
                    else:
                        for f in Fs:
                            for r in Rs:
                                for p in Ps:
                                    header = sequences[0]
                                    f_name = f[2]
                                    f_ver = f[5]
                                    p_ver = p[5]
                                    r_ver = r[5]
                                    f_res = TOOLS.match_fuzzily(f_ver, sequences[1], deletions, insertions, substitutions)
                                    r_res = TOOLS.match_fuzzily(r_ver, sequences[2], deletions, insertions, substitutions)

                                    if (f_res == None) or (r_res == None):
                                        start = None
                                        end = None
                                        amplicon = ""
                                        amplicon_length = 0
                                        f_match = ""
                                        r_match = ""
                                        p_match = ""
                                        PPC = 0
                                    
                                    else:
                                        Forwards = {}
                                        if type(f_res) == type(tuple()):
                                            Forwards[0] = (f_res[0], f_ver, 0) # (start, match, distance)
                                        else:
                                            for f_i in range(len(f_res)):
                                                Forwards[f_i] = (f_res[f_i].start, f_res[f_i].matched, f_res[f_i].dist)
                                        Reverses = {}
                                        if type(r_res) == type(tuple()):
                                            Reverses[0] = (r_res[0], r_ver, 0)
                                        else:
                                            for r_i in range(len(r_res)):
                                                Reverses[r_i] = (r_res[r_i].start, r_res[r_i].matched, r_res[r_i].dist)
                                        matches = {}
                                        for k_f, v_f in Forwards.items():
                                            start = v_f[0]
                                            for k_r, v_r in Reverses.items():
                                                end = (len(sequences[1]) - 1) - v_r[0]
                                                if end < start:
                                                    matches[f"{k_f}:{k_r}"] = False
                                                amplicon = sequences[1][start:end]
                                                if len(amplicon) > 850:
                                                    matches[f"{k_f}:{k_r}"] = False
                                                else:
                                                    p_res = TOOLS.match_fuzzily(p_ver, amplicon, deletions, insertions, substitutions)
                                                    if p_res == None:
                                                        matches[f"{k_f}:{k_r}"] = False
                                                    else:
                                                        matches[f"{k_f}:{k_r}"] = True
                                        target_dist = np.Inf
                                        n_match = 0
                                        for k, v in matches.items():
                                            if v:
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
                                                    if amplicon_length > 850:
                                                        n_match -= 1
                                                        start = None
                                                        end = None
                                                        amplicon = ""
                                                        amplicon_length = 0
                                                        f_match = ""
                                                        r_match = ""
                                                        PPC = 0
                                        if n_match <= 0:
                                            start = None
                                            end = None
                                            amplicon = ""
                                            amplicon_length = 0
                                            f_match = ""
                                            r_match = ""
                                            PPC = 0
                                        else:
                                            PPC = TOOLS.calculate_PPC(F_primer=f_ver,
                                                                    F_match=f_match,
                                                                    R_primer=r_ver,
                                                                    R_match=r_match)

                                    res.append([f_name, f_ver, p_ver,
                                                r_ver, header, amplicon,
                                                amplicon_length, start, end, PPC])
            res_df = pd.DataFrame(res, columns = col_list)
            del res
            return res_df

        self.fname = fname
        self.csv_fname = csv_fname
        self.deletions = deletions
        self.insertions = insertions
        self.substitutions = substitutions

        unique_groups = self.primers["ID"].unique()
        summary = pd.DataFrame(columns = self.SUMMARY_qPCR_COL_LIST)
        os.makedirs(self.savedir, exist_ok = True)
        print("Running Benchmark")
        cluster = LocalCluster(n_workers = self.nCores, threads_per_worker = 4, silence_logs=logging.ERROR)
        client = Client(cluster, timeout = 120)

        for group in tqdm(unique_groups):
            def help_analyse(x):
                return analyse(x, Fs, Rs, Ps, self.BENCHMARK_qPCR_COL_LIST,
                            self.deletions, self.insertions, self.substitutions)
            Fs = self.primers.loc[(self.primers["ID"] == group) & (self.primers["Type"] == "F"),:].values
            Rs = self.primers.loc[(self.primers["ID"] == group) & (self.primers["Type"] == "R"),:].values
            Ps = self.primers.loc[(self.primers["ID"] == group) & (self.primers["Type"] == "P"),:].values
            print(f"Processing group {group}\n")
            futures = client.map(help_analyse, self.chunkpaths)
            progress(futures)
            result_chunks = client.gather(futures)
            group_df = pd.concat(result_chunks)
            group_df.reset_index(drop = True, inplace = True)
            print("\nPerformance computed, generating group summary\n")
            group_stats = generate_group_summary(group_df, group, self.SUMMARY_qPCR_COL_LIST)
            summary = summary.append(group_stats)
            client.cancel(futures)
            del group_stats
            del result_chunks
            print("Summary generated, saving group benchmark to Feather\n")
            group_df.to_feather(os.path.join(self.tmpdir, f"{group}_"+self.fname), compression = "uncompressed")
            print(f"Benchmark results saved to {os.path.join(self.savedir, self.fname)}\n")
            del group_df

        summary.to_csv(os.path.join(self.savedir, self.csv_fname), index = False)
        print(f"Benchmark results saved to {os.path.join(self.savedir, self.fname)}\n")
        print(f"Benchmark summary saved to {os.path.join(self.savedir, self.csv_fname)}\n")