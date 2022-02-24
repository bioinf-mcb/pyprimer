from utils import sequence
import gc
gc.enable()

mode = sequence.READ_MODES.DIRECTORY
primers = sequence.PCRPrimer(mode = mode)
primers_df = primers.describe_primers("/storage/BINF/PawelLab/mkowalski/mkowalski/primers/data/PROPER_PYPRIMER/MCB_design/primers/")
primers_df.to_csv("/storage/BINF/PawelLab/mkowalski/mkowalski/primers/data/PROPER_PYPRIMER/MCB_design/primers.csv", index = False)

# del primers_df
# gc.collect()

# mode = sequence.READ_MODES.FASTA
# data = sequence.Sequence(mode = mode)
# sequence_df = data.describe_sequences(seqs_path = "/storage/BINF/PawelLab/mkowalski/mkowalski/primers/data/PROPER_PYPRIMER/FROZEN_DATASET_5pN.fasta", verbose = True)
# sequence_df.to_csv("/storage/BINF/PawelLab/mkowalski/mkowalski/primers/data/PROPER_PYPRIMER/sequence_df.csv", index = False)
