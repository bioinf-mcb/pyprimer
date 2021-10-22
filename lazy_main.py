from pyprimer.modules.Benchmark import benchmark
from pyprimer.utils.sequence import Sequence, READ_MODES, PCRPrimer
from pyprimer.utils import collect_results
import os
import pandas as pd
import time
from tqdm import tqdm
import threading

if __name__ == "__main__":
    resdir = "./results_test"
    tmpdir = "./tmp"

    # Describing primers from csv
    primers_df = PCRPrimer(READ_MODES.CSV).describe_primers("./artic_primers.csv")

    # primers_df = "./tests/primers_df.csv"
    sequence_path = "/data/SARS-CoV-2/sample/sequences.fasta"
    # total_seqs = 3.6e6
    total_seqs = 2767

    chunks = 1000

    sequences_gen = Sequence(READ_MODES.FASTA).describe_sequences(sequence_path, chunk_size = chunks)

    # Regular computations
    # runner = benchmark.Benchmark(primers_df, sequences_gen, savedir = f'{resdir}', tmpdir = tmpdir, nCores = 20)
    # runner.qPCR_performance(deletions = 1, insertions = 1, substitutions = 1)

    # Lazy loading
    with tqdm(enumerate(sequences_gen), total=total_seqs//chunks) as t:
        for idx, sequences_df in t:
            runner = benchmark.Benchmark(primers_df, sequences_df, savedir = f'{resdir}/{idx}', tmpdir = tmpdir, nCores = 8, verbose=False)
            runner.qPCR_performance(deletions = 1, insertions = 1, substitutions = 1, pbar=t)

            with open(f'{resdir}/{idx}/info.csv', 'w') as f:
                f.write(f"samples_num,{len(sequences_df)}")

    collect_results(resdir)
    print("done")
