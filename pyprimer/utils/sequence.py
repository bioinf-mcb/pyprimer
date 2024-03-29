# -*- coding: utf-8 -*-
import threading
import pandas as pd
import os
import re
from collections import Counter, Iterator
from threading import Thread
import numpy as np
from numba import jit
from enum import Enum
from .essentials import Essentials
import tqdm
from tqdm import trange
from time import sleep

class READ_MODES(Enum):
    CSV = 1
    FASTA = 2
    DIRECTORY = 3

class PCRPrimer(object):
    def __init__(self, mode):
        """
        Initialize the class with specific read mode (formatting instructuons in README.md).
        mode - object of type 'string' that sets the mode of reading primers:
            'csv' - tabularized format of primers and metadata in csv format.
            'directory' - many fasta files with primers in one directory, that must be separated from directory with sequences
        """
        self.mode = mode

    def describe_primers(self,
                         primers_path,
                         Na=50,
                         K=0,
                         Tris=0,
                         Mg=0,
                         dNTPs=0,
                         shift=0,
                         nn_table=None,
                         tmm_table=None,
                         imm_table=None,
                         de_table=None,
                         dnac1=25,
                         dnac2=25,
                         salt_correction=False):
        """
        Method that reads primers and parse them into dataframe with all the metadata
        """
        if self.mode == READ_MODES.CSV:
            primers_df = pd.read_csv(primers_path)
            seriesdict = {"Origin": [],
                          "Target": [],
                          "ID": [],
                          "Name": [],
                          "Sequence": [],
                          "Version": [],
                          "Type": [],
                          "Length": [],
                          "GC(%)": [],
                          "AT(%)": [],
                          "Tm": []}
            for _, row in primers_df.iterrows():
                name_ = row['Name']
                origin_, target_, type_ = row["Group"], '', row["Type"]
                sequence_ = row["Sequence"]
                versions = Essentials.get_all_possible_versions(sequence_)
                length_ = len(row["Sequence"])

                for version_ in versions:
                    seriesdict["Origin"].append(origin_)
                    seriesdict["Target"].append(target_)
                    seriesdict["ID"].append("{}_{}".format(origin_, target_))
                    seriesdict["Name"].append(name_)
                    seriesdict["Sequence"].append(sequence_)
                    seriesdict["Version"].append(version_)

                    gc_ = Essentials.GCcontent(version_)
                    seriesdict["GC(%)"].append(gc_)
                    seriesdict["AT(%)"].append(100 - gc_)

                    tm_ = Essentials.Tm(seq=version_, GC=gc_, 
                                    Na=Na,
                                    K=K,
                                    Tris=Tris,
                                    Mg=Mg,
                                    dNTPs=dNTPs,
                                    shift=shift,
                                    nn_table=nn_table,
                                    tmm_table=tmm_table,
                                    imm_table=imm_table,
                                    de_table=de_table,
                                    dnac1=dnac1,
                                    dnac2=dnac2,
                                    salt_correction=salt_correction)
                    seriesdict["Tm"].append(tm_)
                    seriesdict["Type"].append(type_)
                    seriesdict["Length"].append(length_)
            return pd.DataFrame(seriesdict)

        elif self.mode == READ_MODES.FASTA:
            # TODO for meged fasta file with all primers
            raise NotImplementedError(
                "This variant of ReadPrimers method is yet not implemented")

        elif self.mode == READ_MODES.DIRECTORY:
            return self._describe_dir(primers_path,
                                      Na=Na,
                                      K=K,
                                      Tris=Tris,
                                      Mg=Mg,
                                      dNTPs=dNTPs,
                                      shift=shift,
                                      nn_table=nn_table,
                                      tmm_table=tmm_table,
                                      imm_table=imm_table,
                                      de_table=de_table,
                                      dnac1=dnac1,
                                      dnac2=dnac2,
                                      salt_correction=salt_correction)
        else:
            raise ValueError(
                "Unspecified {} mode, use 'csv', 'directory' or 'fasta' instead".format(self.mode))

    def _describe_dir(self, primers_path, **kwargs):
        primers_df = pd.DataFrame(columns=[
                                  "Origin", "Target", "ID", "Name", "Sequence", "Version", "Type", "Length", "GC(%)", "AT(%)", "Tm"])
        filelist = os.listdir(primers_path)
        groups = {}
        for f in filelist:
            seriesdict = {"Origin": [],
                          "Target": [],
                          "ID": [],
                          "Name": [],
                          "Sequence": [],
                          "Version": [],
                          "Type": [],
                          "Length": [],
                          "GC(%)": [],
                          "AT(%)": [],
                          "Tm": []}

            headers = []
            seqs = []
            with open(os.path.join(primers_path, f), "r") as fh:
                for line in fh:
                    if '>' in line:
                        headers.append(line.strip()[1:])
                    else:
                        seqs.append(line.strip())

            for i in range(len(headers)):
                name_ = headers[i]
                origin_, target_, type_ = headers[i].split("|")
                sequence_ = seqs[i]
                versions = Essentials.get_all_possible_versions(sequence_)
                length_ = len(seqs[i])

                for version_ in versions:
                    seriesdict["Origin"].append(origin_)
                    seriesdict["Target"].append(target_)
                    seriesdict["ID"].append("{}_{}".format(origin_, target_))
                    seriesdict["Name"].append(name_)
                    seriesdict["Sequence"].append(sequence_)
                    seriesdict["Version"].append(version_)

                    gc_ = Essentials.GCcontent(version_)
                    seriesdict["GC(%)"].append(gc_)
                    seriesdict["AT(%)"].append(100 - gc_)

                    tm_ = Essentials.Tm(seq=version_, GC=gc_, **kwargs)
                    seriesdict["Tm"].append(tm_)
                    seriesdict["Type"].append(type_)
                    seriesdict["Length"].append(length_)
            groups[f] = seriesdict
        for key, item in groups.items():
            temp_df = pd.DataFrame(item)
            primers_df = primers_df.append(temp_df)
        return primers_df

class QueuedGenerator(Iterator):
    def __init__(self, generator):
        super().__init__()
        self.gen = generator
        self.queue = []

        thr1 = Thread(target=self._wrapper)
        thr1.start()
        # sleep(0.5)
        # thr2 = Thread(target=self._wrapper)
        # thr2.start()
        self.threads = [thr1]

    def _wrapper(self):
        try:
            self.queue.append(next(self.gen))
        except StopIteration:
            pass

    def __next__(self):
        try:
            thr = self.threads.pop()
            thr.join()
            res = self.queue.pop()

            thr = Thread(target=self._wrapper)
            thr.start()
            self.threads.append(thr)
            return res 
        except IndexError:
            raise StopIteration

class Sequence(object):
    def __init__(self, mode):
        self.mode = mode
    
    def _process_fastalist(self, fastalist):
        seriesdict = {"Header": [], "Sense Sequence": [],
                          "Antisense Sequence": [], "Length": [], "N(%)": []}
            
        idx = 0
        for element in fastalist:
            element_list = element.split("\n")
            header = element_list[0].replace(">", "")
            sequence = ("").join(
                element_list[1:]).replace("-", "N").upper()
            del element_list
            try:
                length_sequence = len(sequence)

                antisense_sequence = Essentials.Antisense(sequence[::-1])
                n_percent = np.multiply(
                    np.divide(Counter(sequence).pop("N", 0), length_sequence), 100)
                seriesdict["Antisense Sequence"].append(antisense_sequence)
                seriesdict["N(%)"].append(n_percent)

                seriesdict["Length"].append(length_sequence)
                seriesdict["Header"].append(header)
                seriesdict["Sense Sequence"].append(sequence)
            except:
                os.makedirs("./logs/",  mode=755, exist_ok=True)
                with open("./logs/Problematic_{}_{}.txt".format(header[-14:], idx), "w") as f:
                    f.write(sequence)
                raise KeyError(
                    'Unexpected character in {} [index {}]'.format(header, idx))
            idx += 1
        return pd.DataFrame(seriesdict)

    def _craft_generator(self, seqs_path, chunk_size):
        def __to_df(fastas):
            seqs_df = pd.DataFrame(
                columns=["Header", "Sense Sequence", "Antisense Sequence", "Length", "N(%)"])
            seqs_df = seqs_df.append(self._process_fastalist(fastas))
            return seqs_df

        fastas = []
        with open(seqs_path, "r") as fh:
            tmp = ''
            for line in fh:
                if line[0]=='>':
                    if tmp != '':
                        fastas.append(tmp)
                        if len(fastas) >= chunk_size:
                            yield __to_df(fastas)
                            fastas = []

                    tmp = line
                else:
                    tmp += line

        fastas.append(tmp)
        yield __to_df(fastas)

    def describe_sequences(self, seqs_path, seqs_no=None, verbose=False, full=True, chunk_size=None):
        """Calculate description of the input sequences

        Args:
            seqs_path (str): Path to the seqences
            seqs_no (int, optional): Take only the first `seqs_no` number of sequences. Defaults to None.
            verbose (bool): Flag that enables the progressbar

        Returns:
            pd.DataFrame: Summary of the sequences
        """
        if self.mode == READ_MODES.CSV:
            # TODO for DataFrame with columns "Header" "Sequence"
            # seqs_df = pd.read_csv(self.seqs_path)
            # self.dataframe = seqs_df
            raise NotImplementedError(
                "This variant of ReadSequences method is yet not implemented")

        elif self.mode == READ_MODES.FASTA:
            if chunk_size is None:
                with open(seqs_path, "r") as fh:
                    fasta = fh.read()

                fastalist = fasta.split("\n>")
                seqs_df = pd.DataFrame(
                    columns=["Header", "Sense Sequence", "Antisense Sequence", "Length", "N(%)"])
                seqs_df = seqs_df.append(self._process_fastalist(fastalist))
                return seqs_df
            else:
                return QueuedGenerator(self._craft_generator(seqs_path, chunk_size))

        elif self.mode == READ_MODES.DIRECTORY:
            raise NotImplementedError(
                "This variant of ReadSequences method is yet not implemented")

        else:
            raise ValueError(
                "Unspecified {} mode, use 'csv', 'directory' or 'fasta' instead".format(self.mode))
