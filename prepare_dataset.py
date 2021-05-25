from pyprimer.utils import sequence
import psutil
from tqdm import trange
import pandas as pd
import os
nproc = psutil.cpu_count()
# processor = sequence.Sequence(sequence.READ_MODES.FASTA)
# sequence_df = processor.describe_sequences(seqs_path = "../FROZEN_DATASET_NOGAPS_5pN.fasta", verbose = True)
# sequence_df.to_csv("sequence_df.csv", index = False)
sequence_df = pd.read_csv("../sequence_df.csv")
os.makedirs("../tmp", exist_ok = True)
chunk_size = (sequence_df.shape[0] // nproc) + 1
id_min = 0
id_max = chunk_size
for chunk in trange(nproc):
    if id_max < sequence_df.shape[0]:
        part = sequence_df.iloc[id_min:id_max, :]
    elif id_min > sequence_df.shape[0]:
        break
    else:
        part = sequence_df.iloc[id_min:, :]
    id_min += chunk_size
    id_max += chunk_size
    fpath = os.path.join("../tmp",f"chunk_{chunk}.csv")
    part.to_csv(fpath, index = False)
print("done")