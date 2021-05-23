from pyprimer.modules.Benchmark import benchmark
import os
import pandas as pd
import psutil

if __name__ == "__main__":
    resdir = "/storage/BINF/PawelLab/mkowalski/mkowalski/GitHub/"
    tmpdir = os.path.join(resdir, "tmp")
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

    print("Reading primers\n")
    primers_df = pd.read_csv(f"{resdir}primers_df.csv")

    print("Reading sequences\n")
    sequences_df = pd.read_csv(f"{resdir}sequence_df.csv")#sequence_df.csv
    # sequences_df = sequences_df.iloc[:1000,:]

    #print("Sampling sequences\n")
    #sequences_df = sequences_df.iloc[:10001,:]
    #sequences_df.to_csv(f"{resdir}sample_set.csv", index = False)
    print("Running benchmark")
    runner = benchmark.Benchmark(primers_df, sequences_df, savedir = resdir, tmpdir = tmpdir, nCores = psutil.cpu_count())
    runner.qPCR_performance()
    print("done")
