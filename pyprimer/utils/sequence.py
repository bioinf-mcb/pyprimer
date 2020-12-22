# -*- coding: utf-8 -*-
import pandas as pd
import os
import re
from collections import Counter
import numpy as np
from numba import jit
from enum import Enum
from .essentials import Essentials
import tqdm

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
            # TODO for DataFrame with columns "Header" "Sequence"
            # primers_df = pd.read_csv(primers_path)
            # self.dataframe = primers_df
            raise NotImplementedError(
                "This variant of ReadPrimers method is yet not implemented")

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


class Sequence(object):
    def __init__(self, mode):
        self.mode = mode

    def describe_sequences(self, seqs_path, seqs_no=None, verbose=False, full=True):
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
            seqs_df = pd.DataFrame(
                columns=["Header", "Sense Sequence", "Antisense Sequence", "Length", "N(%)"])
            with open(seqs_path, "r") as fh:
                fasta = fh.read()
            fastalist = fasta.split("\n>")
            seriesdict = {"Header": [], "Sense Sequence": [],
                          "Antisense Sequence": [], "Length": [], "N(%)": []}
            
            if seqs_no is not None:
                fastalist = fastalist[:seqs_no]
            if verbose:
                fastalist = tqdm.tqdm(fastalist)
                
            idx = 0
            for element in fastalist:
                element_list = element.split("\n")
                header = element_list[0].replace(">", "")
                sequence = ("").join(
                    element_list[1:]).replace("-", "N").upper()
                del element_list
                try:
                    length_sequence = len(sequence)

                    if full:
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
                    with open("./Problematic_{}_{}.txt".format(header[-14:], idx), "w") as f:
                        f.write(sequence)
                    raise KeyError(
                        'Unexpected character in {} [index {}]'.format(header, idx))
                idx += 1
            seqs_df = seqs_df.append(pd.DataFrame(seriesdict))

            return seqs_df

        elif self.mode == READ_MODES.DIRECTORY:
            raise NotImplementedError(
                "This variant of ReadSequences method is yet not implemented")

        else:
            raise ValueError(
                "Unspecified {} mode, use 'csv', 'directory' or 'fasta' instead".format(self.mode))


if __name__ == "__main__":
    test_primer = PCRPrimer("directory")
    test_primer.DescribePrimers("./data/primers")
    test_sequence = Sequence("fasta")
    test_sequence.DescribeSequences("./data/merged.fasta")
    test_primer.dataframe.to_csv("primers_test_df.csv", index=False)
    test_sequence.dataframe.to_csv("sequence_test_df.csv", index=False)
