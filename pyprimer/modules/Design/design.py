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

class Design(object):
    def __init__(self, msa, path_to_genes):
        self.msa = msa
        self.path_to_genes = path_to_genes
        print("Preparing gene dictionary\n")
        gene_list = os.listdir(path_to_genes)
        gene_dict = {}
        for element in gene_list:
            gene_name = element.split(".")[0]
            gene_file = os.path.join(path_to_genes, element)
            gene_dict[gene_name] = gene_file
        self.gene_dict = gene_dict
        print("Gene dictionary prepared\n")

    def align_genes(self, alignment_dir, n_proc = 8):    
        self.alignment_dir = alignment_dir
        self.genes_alignment_dir = os.path.join(self.alignment_dir, "aligned_genes")
        os.makedirs(self.genes_alignment_dir, exist_ok = True)
        flist = os.listdir(self.alignment_dir)
        for f in tqdm(flist):
            for k, v in self.gene_dict.items():
                tmp_path = os.path.join(self.genes_alignment_dir,f"{k}_aligned.fasta")
                cmd = f"mafft --thread {n_proc} --add {v} {f} > {tmp_path}}"
                os.system(cmd)
                sleep(2)
                clip_msa(tmp_path)

    def construct_matrices(self, differential_design = False, separator = "\t" positive_samples = None,
                           negative_samples = None, stride = 1, sliding_windows = [19,20,21],
                           consensus_threshold = 0.15):
        self.differential_design = differential_design
        self.separator = separator
        self.positive_samples = positive_samples
        self.negative_samples = negative_samples
        self.stride = stride
        self.windows = sliding_windows
        self.threshold = consensus_threshold
        dataset = []
        if self.differential_design:
        
        else:



        