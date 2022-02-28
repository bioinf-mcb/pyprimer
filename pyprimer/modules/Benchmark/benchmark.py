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
from primer3.bindings import calcEndStability, calcHeterodimer
import logging
warnings.simplefilter('ignore', tables.NaturalNameWarning)
warnings.simplefilter('ignore', tables.PerformanceWarning)
warnings.simplefilter("ignore", UserWarning)
warnings.filterwarnings("ignore")

class Benchmark(object):


    BENCHMARK_qPCR_COL_LIST = [
                "Group ID",
                "Prime ID",
                "F Primer Name",
                "F Primer Version",
                "P Probe Version",
                "R Primer Version",
                "Sequence Header",
                "Amplicon Sense",
                "Amplicon Sense Length",
                "Amplicon Sense Start",
                "Amplicon Sense End",
                "F Heterodimer Structure",
                "P Heterodimer Structure",
                "R Heterodimer Structure",
                "F Heterodimer Tm(Celsius)",
                "P Heterodimer Tm(Celsius)",
                "R Heterodimer Tm (Celsius)",
                "F Heterodimer Gibbs(cal)",
                "P Heterodimer Gibbs(cal)",
                "R Heterodimer Gibbs(cal)",
                "F Prime",
                "P Prime",
                "R Prime",
                "System Prime",
                "Primer Pair Coverage"]

    SUMMARY_qPCR_COL_LIST = [
                        "Primer Group",
                        "F Version",
                        "P Version",
                        "R Version",
                        "Mean PPC",
                        "Sequences matched(%)",
                        "Misprime ratio(%)"]

    def __init__(self, primer_df, sequence_df, n_sequences, savedir = "./results", tmpdir = "./tmp", nCores = 4):
        self.primers = pd.read_csv(primer_df)
        self.nCores = nCores
        self.tmpdir = tmpdir
        self.savedir = savedir
        print("Preparing data chunks")
        self.n_sequences = n_sequences
        try:
            print("Revoking chunks")
            flist = os.listdir(self.tmpdir)
            if len(flist) >= 1:
                chunkpaths = [os.path.join(self.tmpdir,f) for f in flist if "chunk_" in f]
                self.chunkpaths = chunkpaths
            else:
                raise ValueError()
        except:
            print("Revoke failed")
            sequence_df = pd.read_csv(sequence_df)
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
            del sequence_df

    def qPCR_performance(self, deletions = None, insertions = None, substitutions = None, distance = 3,
                         fname = 'pyprimer_benchmark.feather', csv_fname = "pyprimer_summary.csv", anneal_temp_c = 60,
                         mv_conc_mM = 50, dv_conc_mM = 12, dntp_conc_mM = 0.8, dna_conc_nM = 50, max_loop = 30,
                         gibbs_threshold = -9000, template_offset = 10):
        # def generate_group_summary(group, col_list, PPC_s, MisRatio_s, SeqMatch_s):
        #     v_stats = dict((key,[]) for key in col_list)
        #     v_stats["Primer Group"].append(group)
        #     v_stats["Mean PPC"].append(np.mean(np.nan_to_num(PPC_s)))
        #     v_stats["Sequences matched(%)"].append(np.mean(SeqMatch_s)*100)
        #     v_stats["Misprime ratio(%)"].append(np.mean(np.nan_to_num(MisRatio_s)))
        #     group_stats = pd.DataFrame(v_stats)
        #     return group_stats

        def analyse(sequences_path, Fs, Rs, Ps, col_list,
                    deletions, insertions, substitutions,
                    distance, anneal_temp_c, mv_conc_mM,
                    dv_conc_mM, dntp_conc_mM, dna_conc_nM,
                    max_loop, gibbs_threshold, template_offset):
            res = []
            with open(sequences_path, "r", newline='') as csvfile:
                seqreader = csv.reader(csvfile, delimiter = ',', quotechar ='"')
                PPC_ = []
                SeqMatch_ = []
                MisRatio_ = []
                for sequences in seqreader:
                    if sequences[0] == "Header":
                        pass
                    else:
                        PPC_list = []
                        SeqMatch_list = []
                        MisRatio_list = []
                        for f in Fs: 
                            for r in Rs:
                                for p in Ps:
                                    group_id = TOOLS.id_generator()
                                    f_heterodimer = None
                                    p_heterodimer = None
                                    r_heterodimer = None
                                    f_tm = None
                                    p_tm = None
                                    r_tm = None
                                    f_dg = None
                                    p_dg = None
                                    r_dg = None
                                    header = sequences[0]
                                    f_name = f[2]
                                    f_ver = f[5]
                                    p_ver = p[5]
                                    r_ver = r[5]
                                    f_res = TOOLS.match_fuzzily(f_ver, sequences[1], deletions, insertions, substitutions, distance)
                                    r_res = TOOLS.match_fuzzily(r_ver, sequences[2], deletions, insertions, substitutions, distance)

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
                                        if f_res != None:
                                            for f_i in range(len(f_res)):
                                                Forwards[f_i] = (f_res[f_i].start, f_res[f_i].end, f_res[f_i].matched, f_res[f_i].dist)
                                        Reverses = {}
                                        if r_res != None:
                                            for r_i in range(len(r_res)):
                                                Reverses[r_i] = (r_res[r_i].start, r_res[r_i].end, r_res[r_i].matched, r_res[r_i].dist)
                                        Probes = {}
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
                                                    p_res = TOOLS.match_fuzzily(p_ver, amplicon, deletions, insertions, substitutions, distance)
                                                    if p_res == None:
                                                        matches[f"{k_f}:{k_r}"] = False
                                                    else:
                                                        for p_i in range(len(p_res)):
                                                            Probes[p_i] = (p_res[p_i].start, p_res[p_i].end, start)
                                                        matches[f"{k_f}:{k_r}"] = True

                                        n_match = 0
                                        binding_id = 0
                                        for k, v in matches.items():
                                            if v:
                                                n_match += 1
                                                klist = k.split(":")
                                                k_f = int(klist[0])
                                                k_r = int(klist[1])
                                                f_good = Forwards[k_f]
                                                r_good = Reverses[k_r]
                                                for k_p in Probes.keys():
                                                    p_good = Probes[k_p]
                                                    start = f_good[0]
                                                    end = (len(sequences[1]) - 1) - r_good[0]

                                                    f_start = f_good[0] - template_offset
                                                    f_end = f_good[1] + template_offset
                                                    f_match = f_good[2]
                                                    r_start = r_good[0] - template_offset
                                                    r_end = r_good[1] + template_offset
                                                    r_match = r_good[2]
                                                    amplicon_start = p_good[2]
                                                    p_start = (amplicon_start + p_good[0]) - template_offset
                                                    p_end = (amplicon_start + p_good[1]) + template_offset

                                                    f_template = TOOLS.extract_template(sequences[2], f_start, f_end)
                                                    p_template = TOOLS.extract_template(sequences[2], p_start, p_end)
                                                    r_template = TOOLS.extract_template(sequences[1], r_start, r_end)

                                                    f_bind = calcHeterodimer(f_ver, 
                                                                                f_template,
                                                                                mv_conc = mv_conc_mM,
                                                                                dv_conc = dv_conc_mM,
                                                                                dntp_conc = dntp_conc_mM,
                                                                                dna_conc = dna_conc_nM,
                                                                                temp_c = anneal_temp_c,
                                                                                max_loop = max_loop,
                                                                                output_structure = True)
                                                    f_stability = calcEndStability(f_ver,
                                                                                    f_template,
                                                                                    mv_conc = mv_conc_mM,
                                                                                    dv_conc = dv_conc_mM,
                                                                                    dntp_conc = dntp_conc_mM,
                                                                                    dna_conc = dna_conc_nM,
                                                                                    temp_c = anneal_temp_c,
                                                                                    max_loop = max_loop)

                                                    p_bind = calcHeterodimer(p_ver, 
                                                                                p_template,
                                                                                mv_conc = mv_conc_mM,
                                                                                dv_conc = dv_conc_mM,
                                                                                dntp_conc = dntp_conc_mM,
                                                                                dna_conc = dna_conc_nM,
                                                                                temp_c = anneal_temp_c,
                                                                                max_loop = max_loop,
                                                                                output_structure = True)
                                                    p_stability = calcEndStability(p_ver,
                                                                                    p_template,
                                                                                    mv_conc = mv_conc_mM,
                                                                                    dv_conc = dv_conc_mM,
                                                                                    dntp_conc = dntp_conc_mM,
                                                                                    dna_conc = dna_conc_nM,
                                                                                    temp_c = anneal_temp_c,
                                                                                    max_loop = max_loop)
                                                    
                                                    r_bind = calcHeterodimer(r_ver, 
                                                                                r_template,
                                                                                mv_conc = mv_conc_mM,
                                                                                dv_conc = dv_conc_mM,
                                                                                dntp_conc = dntp_conc_mM,
                                                                                dna_conc = dna_conc_nM,
                                                                                temp_c = anneal_temp_c,
                                                                                max_loop = max_loop,
                                                                                output_structure = True)
                                                    r_stability = calcEndStability(r_ver,
                                                                                    r_template,
                                                                                    mv_conc = mv_conc_mM,
                                                                                    dv_conc = dv_conc_mM,
                                                                                    dntp_conc = dntp_conc_mM,
                                                                                    dna_conc = dna_conc_nM,
                                                                                    temp_c = anneal_temp_c,
                                                                                    max_loop = max_loop)

                                                    bool_list = []
                                                    bool_list.append(TOOLS.check_correctness(f_bind,
                                                                                             f_stability,
                                                                                             gibbs_threshold,
                                                                                             anneal_temp_c))
                                                    bool_list.append(TOOLS.check_correctness(p_bind,
                                                                                             p_stability,
                                                                                             gibbs_threshold,
                                                                                             anneal_temp_c))
                                                    bool_list.append(TOOLS.check_correctness(r_bind,
                                                                                             r_stability,
                                                                                             gibbs_threshold,
                                                                                             anneal_temp_c))
                                                    if bool_list[0]:
                                                        f_verdict = True
                                                        f_heterodimer = f_bind.ascii_structure
                                                        f_tm = f_bind.tm
                                                        f_dg = f_bind.dg
                                                    if bool_list[1]:
                                                        p_verdict = True
                                                        p_heterodimer = p_bind.ascii_structure
                                                        p_tm = p_bind.tm
                                                        p_gd = p_bind.dg
                                                    if bool_list[2]:
                                                        r_verdict = True
                                                        r_heterodimer = r_bind.ascii_structure
                                                        r_tm = r_bind.tm
                                                        r_dg = r_bind.dg
                                                    if False in bool_list:
                                                        prime = False
                                                    else:
                                                        prime = True
                                                    amplicon = sequences[1][start:end]
                                                    amplicon_length = len(amplicon)
                                                    if (amplicon_length > 850) or (False in bool_list):
                                                        n_match -= 1
                                                        start = None
                                                        end = None
                                                        amplicon = ""
                                                        amplicon_length = 0
                                                        f_heterodimer = ""
                                                        p_heterodimer = ""
                                                        r_heterodimer = ""
                                                        f_tm = None
                                                        p_tm = None
                                                        r_tm = None
                                                        f_dg = None
                                                        r_dg = None
                                                        p_dg = None
                                                        f_verdict = False
                                                        p_verdict = False
                                                        r_verdict = False
                                                        prime = False
                                                        f_match = ""
                                                        p_match = ""
                                                        r_match = ""
                                                        PPC = 0

                                                    else:
                                                        PPC = TOOLS.calculate_PPC(F_primer=f_ver,
                                                                                  F_match=f_match,
                                                                                  R_primer=r_ver,
                                                                                  R_match=r_match)
                                                    PPC_list.append(PPC)
                                                    MisRatio_list.append(prime)
                                                    res.append([group_id, f"{group_id}_{binding_id}",
                                                                f_name, f_ver, p_ver,
                                                                r_ver, header, amplicon,
                                                                amplicon_length, start, end,
                                                                f_heterodimer, p_heterodimer,
                                                                r_heterodimer, f_tm, p_tm, r_tm,
                                                                f_dg, p_dg, r_dg, f_verdict, p_verdict,
                                                                r_verdict, prime, PPC])
                                                    binding_id += 1
                        PPC_num = np.nan_to_num(PPC_list, copy = True)
                        PPC_.append(np.mean(PPC_num))
                        MisRatio_num = np.nan_to_num(MisRatio_list, copy = True)
                        MisRatio_.append(1 - (np.sum(MisRatio_num)/(len(MisRatio_num))))
                        if True in MisRatio_list:
                            SeqMatch_.append(True)
                        else:
                            SeqMatch_.append(False)

            res_df = pd.DataFrame(res, columns = col_list)
            del res
            return (res_df, PPC_, MisRatio_, SeqMatch_)

        self.fname = fname
        self.csv_fname = csv_fname
        self.deletions = deletions
        self.insertions = insertions
        self.substitutions = substitutions
        self.distance = distance
        self.anneal_temp_c = anneal_temp_c
        self.mv_conc_mM = mv_conc_mM
        self.dv_conc_mM = dv_conc_mM
        self.dntp_conc_mM = dntp_conc_mM
        self.dna_conc_nM = dna_conc_nM
        self.max_loop = max_loop
        self.gibbs_threshold = gibbs_threshold
        self.template_offset = template_offset 

        unique_groups = self.primers["ID"].unique()
        summary_ = pd.DataFrame(columns = self.SUMMARY_qPCR_COL_LIST)
        summary = dd.from_pandas(summary_, npartitions=self.nCores)
        del summary_
        os.makedirs(self.savedir, exist_ok = True)
        print("Running Benchmark")
        cluster = LocalCluster(n_workers = self.nCores, threads_per_worker = 4, silence_logs=logging.ERROR)
        client = Client(cluster, timeout = 120)

        for group in tqdm(unique_groups):
            def help_analyse(x):
                return analyse(x, Fs, Rs, Ps, self.BENCHMARK_qPCR_COL_LIST,
                            self.deletions, self.insertions, self.substitutions,
                            self.distance, self.anneal_temp_c, self.mv_conc_mM,
                            self.dv_conc_mM, self.dntp_conc_mM, self.dna_conc_nM,
                            self.max_loop, self.gibbs_threshold, self.template_offset)
            Fs = self.primers.loc[(self.primers["ID"] == group) & (self.primers["Type"] == "F"),:].values
            Rs = self.primers.loc[(self.primers["ID"] == group) & (self.primers["Type"] == "R"),:].values
            Ps = self.primers.loc[(self.primers["ID"] == group) & (self.primers["Type"] == "P"),:].values
            print(f"Processing group {group}\n")
            futures = client.map(help_analyse, self.chunkpaths)
            progress(futures)
            result_chunks = client.gather(futures)
            results_s = []
            PPC_s = []
            MisRatio_s = []
            SeqMatch_s = []
            for element in result_chunks:
                results_s.append(element[0])
                PPC_s += element[1]
                MisRatio_s += element[2]
                SeqMatch_s += element[3]
            group_df = pd.concat(results_s)
            group_df.reset_index(drop = True, inplace = True)
            print("\nPerformance computed, generating group summary\n")

            # group_stats = generate_group_summary(group, self.SUMMARY_qPCR_COL_LIST,
            #                                      PPC_s, MisRatio_s, SeqMatch_s)
            # summary = summary.append(group_stats)
            client.cancel(futures)
            # del group_stats
            del result_chunks
            print("Summary generated, saving group benchmark to Feather\n")
            group_df.to_feather(os.path.join(self.tmpdir, f"{group}_"+self.fname), compression = "uncompressed")
            print(f"Benchmark results saved to {os.path.join(self.tmpdir, group + '_' + self.fname)}\n")
            del group_df

        # summary.to_csv(os.path.join(self.savedir, self.csv_fname), index = False)
        # print(f"Benchmark summary saved to {os.path.join(self.savedir, self.csv_fname)}\n")
