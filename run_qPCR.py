from pyprimer.modules.Benchmark import benchmark
import os
import pandas as pd

if __name__ == "__main__":
    resdir = "/klaster/scratch/mkowalski/primers/data/NEW_CLADES/results/Validated/"
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

    print("Reading primers\n")
    primers_df = pd.read_csv(f"{resdir}primers_df.csv")

    print("Reading sequences\n")
    sequences_df = pd.read_csv(f"{resdir}sample_set.csv")

    #print("Sampling sequences\n")
    #sequences_df = sequences_df.iloc[:10001,:]
    #sequences_df.to_csv(f"{resdir}sample_set.csv", index = False) 
    print("Running benchmark")
    
    runner = benchmark.Benchmark(primers_df, sequences_df, nCores = 16)
    runner.qPCR_performance(savedir = resdir)
    print("done")
