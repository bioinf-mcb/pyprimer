from pyprimer.utils.sequence import PCRPrimer, Sequence, READ_MODES
from pyprimer.utils.sequence import Sequence
from pyprimer.modules.PPC import PPC
import os
import dask.dataframe as dd
import pandas as pd

if __name__ == "__main__":
   resdir = "/klaster/scratch/mkowalski/primers/data/NEW_CLADES/results/Validated/"
   os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

   print("Describing primers\n")
   # test_primer = PCRPrimer(READ_MODES.DIRECTORY)
   # primers_df = test_primer.describe_primers("/klaster/scratch/mkowalski/primers/data/primers/")
   primers_df = pd.read_csv(f"{resdir}primers_df.csv")
   
   print("Describing sequences\n")
   # test_sequence = Sequence(READ_MODES.FASTA)
   # sequences_df = test_sequence.describe_sequences("/klaster/scratch/mkowalski/primers/data/NEW_CLADES/FROZEN_DATASET_NOGAPS_5pN.fasta")
   sequences_df = pd.read_csv(f"{resdir}sequence_df.csv")

   print("Calculating performance\n")
   test_pcr = PPC(
      primers_df,
      sequences_df,
      memsave = True,
      tempdir = resdir,
      fname = "FROZEN_RESULTS_DEBUGGED.h5"
      )

   summary = test_pcr.analyse_primers(
      nCores=32,
      deletions = 0,
      insertions = 0,
      substitutions = 0
      )

   print("Dumping to CSV\n")
   summary.to_csv(
      f"{resdir}FROZEN_RESULTS_SUMMARY_DEBUGGED.csv",
      index = False
      )
   # sequence_df.to_csv(
   #    f"{resdir}sequences.csv",
   #    index = False
   # )

   # primer_df.to_csv(
   #    f"{resdir}primers.csv",
   #    index = False
   # )
