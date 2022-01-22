import os
import numpy as np

fw = open("precisions.dat", "w")

# listdir = os.listdir(".")

# print(listdir)

# for folder in listdir:
for folder in range(1,7):
    folder = str(folder)
    if os.path.isdir(folder):
        precisions_path = f"{folder}/precisions"
        b_fpath = f"{precisions_path}/avg_p.dat"
        fo = open(b_fpath, "r")
        foo = fo.readlines()
        for line in foo:
            fw.write(f"{folder} {line}")
fw.close()