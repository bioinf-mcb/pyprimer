from pyprimer.modules.Benchmark import Benchmark
import os
import pandas as pandas

if __name__ == "__main__":
    resdir = "/klaster/scratch/mkowalski/primers/data/NEW_CLADES/results/Validated/"
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

    print("Reading primers\n")
    primers_df = pd.read_csv(f"{resdir}primers_df.csv")

    print("Reading sequences\n")
    sequences_df = pd.read_csv(f"{resdir}sequence_df.csv")

    print("Sampling sequences\n")
    sequences_df = sequences_df.iloc[:10001,:]

    print("Running benchmark")
    
    runner = Benchmarl(primer_df, sequence_df, nCores = 16)
    runner.qPCR_performance(savedir = resdir)
    print("done")
