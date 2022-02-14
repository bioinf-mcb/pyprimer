import os
import re
import pandas as pd
import numpy as np
from time import sleep
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import trange, tqdm
from scipy.stats import entropy
from primer3.bindings import 
from pyprimer.modules.Design.design_tools import TOOLS
from pyprimer.utils.essentials import Essentials

class TOOLS:
    @staticmethod
    def clip_msa(msa_):
            msa = list(SeqIO.parse(msa_, "fasta"))
            gene = msa[-1]
            nuc = "-"
            start = 0
            while nuc == "-":
                start += 1
                nuc = gene.seq.__str__()[start]
            # print("Start index found")
            end = start
            while nuc != "-":
                end +=1
                nuc = gene.seq.__str__()[end]
            # print("End index found")
            # print(f"{start}:{end}")
            for i in range(len(msa)):
                new_seq = Seq(msa[i].seq.__str__()[start:end])
                msa[i].seq = new_seq
            print(len(msa[i].seq.__str__()))
            SeqIO.write(msa, msa_, "fasta")

    @staticmethod
    def construct_seq_matrix(msa):
        msa = list(SeqIO.parse(msa, "fasta"))
        w = len(msa[0].seq.__str__())
        h = len(msa)
        seq_matrix = np.ndarray(shape=(h, w), dtype = object)
        for i in range(h):
            seq_matrix[i,:] = list(msa[i].seq.__str__().upper().replace("-", "N"))
        return seq_matrix, w, h

    @staticmethod
    def calculate_probability(seq_matrix, w, h):
        probability_matrix = np.ndarray(shape=(4, w))
        prob_dict = Essentials.NUC_ENCODING
        transform_dict = Essentials.IUPAC
        for i in trange(seq_matrix.shape[1]):
            cured_list = []
            for j in range(seq_matrix.shape[0]):
                cured_list += transform_dict[seq_matrix[j,i]]
            length = len(cured_list)
            unique, counts = np.unique(cured_list, return_counts = True)
            zipped = dict(zip(unique, counts))
            for k, v in zipped.items():
                probability_matrix[prob_dict[k],i] = v/length
        return probability_matrix
    
    @staticmethod
    def calculate_position_entropy(probability_matrix , w, h):
        entropy_matrix = np.ndarray(shape=(1, w))
        for i in trange(probability_matrix.shape[1]):
            column = probability_matrix[:,i]
            entropy_matrix[0,i] = entropy(column)
        return(entropy_matrix)
    
    @staticmethod
    def entropy_sliding_window(entropy_matrix, roi_start, roi_end, stride = 1, windows = [19,20,21]):
        findings = {
            "start":[],
            "end":[],
            "length":[],
            "entropy":[],
        }
        for i in trange(0, entropy_matrix.shape[1], stride):
            tmp_finds = {k:np.nan for k in windows}
            for w in windows:
                if i+w > entropy_matrix.shape[1] -1:
                    pass
                else:
                    tmp_finds[w] = np.mean(entropy_matrix[:,i:i+w])
            minimal = min(tmp_finds, key = tmp_finds.get)
            findings["start"].append(i)
            findings["end"].append(i+minimal)
            findings["length"].append(minimal)
            findings["entropy"].append(tmp_finds[minimal])
        return pd.DataFrame(findings)

    @staticmethod
    def generate_consensus(probs, start, end, threshold = 0.15):
        prob_dict = Essentials.NUC_ENCODING
            }
        transform_dict = Essentials.IUPAC
        slice_ = probs[:,start:end+1]
        position_dict = {}
        for i in range(slice_.shape[1]):
            nucs = []
            if slice_[0,i] > threshold:
                nucs.append("A")
            elif slice_[1,i] > threshold:
                nucs.append("C")
            elif slice_[2,i] > threshold:
                nucs.append("T")
            elif slice_[3,i] > threshold:
                nucs.append("G")
            position_dict[i] = nucs
        
        positions = list(position_dict.values())
        consensus = []
        for k, v in position_dict.items():
            for k_, v_ in transform_dict.items():
                if set(v) == set(v_):
                    consensus.append(k_)
        return ("").join(consensus)