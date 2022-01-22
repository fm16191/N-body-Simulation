import os
import numpy as np

fw = open("benchs.dat", "w")

# listdir = os.listdir(".")

# print(listdir)

# for folder in listdir:
for folder in range(7):
    folder = str(folder)
    # if folder.startswith("."): continue
    if os.path.isdir(folder):
        bench_path = f"{folder}/benchs"
        precisions_path = f"{folder}/precisions"
        for b_file in os.listdir(bench_path):
            if not b_file.endswith('.dat'): continue
            compiler = b_file[2:-4]
            b_fpath = f"{bench_path}/{b_file}"
            print(b_fpath)
            fo = open(b_fpath, "r")
            foo = fo.readlines()
            fo.close()
            data = []
            for line in foo:
                data.append(np.float32(line.split(" ")[-1][:-1]))
            min = np.min(data)
            max = np.max(data)
            mean = np.mean(data)
            print(min, max, mean)
            fw.write(f"{folder} {compiler} {min} {max} {mean}\n")

fw.close()
