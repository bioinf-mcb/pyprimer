import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

from pyprimer.modules.Benchmark import benchmark
import os
import pandas as pd
import psutil
import time

if __name__ == "__main__":
    resdir = "/storage/BINF/PawelLab/mkowalski/mkowalski/primers/data/PROPER_PYPRIMER/MCB_design/results/"
    tmpdir = "/storage/BINF/PawelLab/mkowalski/mkowalski/primers/data/PROPER_PYPRIMER/MCB_design/tmp/"
    primers_df = "/storage/BINF/PawelLab/mkowalski/mkowalski/primers/data/PROPER_PYPRIMER/MCB_design/primers.csv"
    sequences_df = "/storage/BINF/PawelLab/mkowalski/mkowalski/primers/data/PROPER_PYPRIMER/test_sequence_df.csv"
    runner = benchmark.Benchmark(primers_df, sequences_df, n_sequences = 1526, savedir = resdir, tmpdir = tmpdir, nCores = 32)
    runner.qPCR_performance(distance = 5, anneal_temp_c = 60, mv_conc_mM = 50, dv_conc_mM = 12)
    print("done")
