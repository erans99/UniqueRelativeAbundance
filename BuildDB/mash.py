import glob
import os
import time
import pandas


def make_sketches(mash_exe, fs, in_dir, mash_out):
    cnt = 0
    for f in fs:
        com = "%s sketch -s 1e4 %s" % (mash_exe, f)
        print("\nRunning (%d of %d):" % (cnt, len(fs)))
        print(com)
        print("at", time.ctime())
        os.system(com)
        cnt += 1
    print("Merging msh files")
    com = "%s paste %s %s/*.msh" % (mash_exe, mash_out, in_dir)
    os.system(com)
    print("Removing tmp file")
    com = "rm -f %s/*.msh" % in_dir
    os.system(com)


def make_dists_csv(df_dists, csv_out):
    if os.path.isfile(csv_out):
        print("Dists matrix exists")
        return
    print("Making dists matrix")
    fs = list(set(df_dists.f1.values))
    res = {}
    for i, f in enumerate(fs):
        tmp = df_dists[df_dists.f1 == f][['f2', 'dist']]
        res[f] = tmp.set_index('f2')['dist']
        if i % 10 == 0:
            print("At %d of %d" % (i, len(fs)))
    res = pandas.DataFrame(res)
    res = res[res.index]
    res.to_csv(csv_out)
    print("Done")
    return


def run():
    singles_path = "./genomes" #path for all representatives. Note: needs write permissions!
    mash_exe = "/usr/bin/mash-Linux64-v2.0/mash"
    f_singles = os.path.join(singles_path, "*.fa")
    out_dir = "./mash" # where all mash output files are created
    mash_out = os.path.join(out_dir, "full_ref.msh")
    res_out = os.path.join(out_dir, "dists.txt")
    csv_out = os.path.join(out_dir, "dists_matrix.csv")

    fs = glob.glob(f_singles)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    if os.path.isfile(mash_out):
        print("%s exists. Not calculating mash" % mash_out)
    else:
        make_sketches(mash_exe, fs, singles_path, mash_out)
    if os.path.isfile(res_out):
        print("%s exists. Not calculating distances" % res_out)
    else:
        print("Creating dists file")
        com = "%s dist %s %s > %s" % (mash_exe, mash_out, mash_out, res_out)
        os.system(com)
        print("dists file complete")

    df_dists = pandas.read_csv(res_out, delim_whitespace=True, header=None)
    df_dists.columns = ['f1', 'f2', 'dist', 'p_val', 'matches']
    df_dists = df_dists[df_dists.f1 != df_dists.f2]

    df_dists['f1'] = df_dists['f1'].apply(os.path.basename, 1)
    df_dists['f2'] = df_dists['f2'].apply(os.path.basename, 1)
    make_dists_csv(df_dists, csv_out)


if __name__ == '__main__':
    run()
