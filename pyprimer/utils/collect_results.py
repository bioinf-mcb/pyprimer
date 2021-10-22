#%%
import glob
import pandas as pd
import os
import shutil
import tqdm

# %%
def _collect_feathers(group_names, resdir, suffix):
    with tqdm.tqdm() as t:
        for name in group_names:
            dframes = []
            idx = 0
            while True:
                files = glob.glob(f"{resdir}/{idx}/{name}*.feather")
                if len(files) == 0:
                    break
                dframes.append(pd.read_feather(files[0]))
                idx += 1
            t.update()

            collected = pd.concat(dframes)
            collected = collected.reset_index(drop=True)
            collected.to_feather(f"{resdir}/{name}{suffix}", compression = "uncompressed")

def _collect_summaries(resdir):
    primer_df = pd.read_csv(f"{resdir}/0/pyprimer_summary.csv")
    num = pd.read_csv(f"{resdir}/0/info.csv", header=None)[1].values[0]
    primer_df["Mean PPC"] *= num
    primer_df["Sequences matched(%)"] *= num

    total_num = num

    idx = 1
    with tqdm.tqdm() as t:
        while True:
            try:
                directory = f"{resdir}/{idx}"
                summary = pd.read_csv(f"{directory}/pyprimer_summary.csv")
                num = pd.read_csv(f"{directory}/info.csv", header=None)[1].values[0]
                primer_df["Mean PPC"] += summary["Mean PPC"] * num
                primer_df["Sequences matched(%)"] += summary["Sequences matched(%)"] * num
                total_num += num
                idx += 1
                t.update()
            except FileNotFoundError:
                break

    primer_df["Mean PPC"] = round(primer_df["Mean PPC"]/total_num, 3)
    primer_df["Sequences matched(%)"] = round(primer_df["Sequences matched(%)"]/total_num, 3)
    primer_df.to_csv(f"{resdir}/pyprimer_summary.csv", index = False)

def collect_results(resdir, suffix="_pyprimer_benchmark.feather"):
    """ Collects partial results and merges them into a single set of files.
    It will be saved into the `resdir` directory, and the rest of the files will be deleted.

    Args:
        resdir(`str`): location of the files
        suffix(`str`): suffix of benchmark filenames (default: '_pyprimer_benchmark.feather')
    """

    group_names = []

    for fname in glob.glob(f"{resdir}/0/*.feather"):
        fname = os.path.split(fname)[-1]
        group_name = fname.replace(suffix, "")
        group_names.append(group_name)

    print("Collecting feathers")
    print(group_names)
    _collect_feathers(group_names, resdir, suffix)

    print("Collecting summaries")
    _collect_summaries(resdir)
    idx = 0
    while True:
        try:
            shutil.rmtree(f"{resdir}/{idx}")
            idx += 1
        except FileNotFoundError:
            break