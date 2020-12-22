from pyprimer.modules import PPC
from pyprimer.utils.sequence import PCRPrimer, Sequence, READ_MODES
from pyprimer.modules.PPC.ppc_tools import TOOLS
import pandas as pd
import h5py
import pytest

data_dir = "/bi/aim/scratch/afrolova/COVID19/github/pyprimer/data"

@pytest.fixture()
def read_squences():
    test_primer = PCRPrimer(READ_MODES.DIRECTORY)
    primer_df = test_primer.describe_primers(f"{data_dir}/primers")

    test_sequence = Sequence(READ_MODES.FASTA)
    sequences_df = test_sequence.describe_sequences(f"{data_dir}/merged.fasta")

    return primer_df, sequences_df

def test_fuzzy_match():
    literal = TOOLS.match_fuzzily("TCGACATCACC", "TCGACATCACCA")
    assert literal[0] == 0 and literal[1] == "TCGACATCACC"
    
    gappy = TOOLS.match_fuzzily("TGGACATCACC", "TCGACATCACCA")
    assert gappy[0] == 0 and gappy[1] == "TCGACATCACC"

    wrong = TOOLS.match_fuzzily("ACCGTAT", "TCGACATCACCA")
    assert wrong[0] is None and wrong[1] == ''

    edgy = TOOLS.match_fuzzily("TCGAGACC", "ACTTGACATCGACATCACCA")
    assert edgy[0] == 8 and edgy[1] == "TCGACATC"


def test_PPC_calculation(read_squences):
    primer_df, sequences_df = read_squences
    
    test_pcr = PPC(primer_df, sequences_df.sample(100, random_state=42))
    summary = test_pcr.analyse_primers(nCores=8, deletions=1, insertions=0, substitutions=2)
    summary.to_csv("tests/test_tmp/summary_df.csv", index = False)

    assert pd.read_csv("tests/test_tmp/summary_df.csv").equals(pd.read_csv("tests/goal_tmp/summary_df.csv"))


def test_PPC_temp_memory(read_squences):
    primer_df, sequences_df = read_squences
    
    test_pcr = PPC(primer_df, sequences_df.sample(100, random_state=42), memsave=True, tempdir="tests/test_tmp")
    _ = test_pcr.analyse_primers(nCores=8, deletions=1, insertions=0, substitutions=2)

    for key in h5py.File("tests/goal_tmp/PCRBenchmark.h5", "r").keys():
        goal_df = pd.read_hdf("tests/goal_tmp/PCRBenchmark.h5", key=key)
        tested_df = pd.read_hdf("tests/test_tmp/PCRBenchmark.h5", key=key)
    assert tested_df.equals(goal_df)


if __name__=="__main__":
    test_PPC_temp_memory()