import numpy as np
from sys import argv
from os.path import isdir
from os import listdir

def compare(f0_filename, f1_filename, output_filename=None):
    f0 = open(f0_filename)
    f00 = f0.readlines()
    f0.close

    f1 = open(f1_filename)
    f11 = f1.readlines()
    f1.close

    print(f"Comparing {f0_filename} - {f1_filename} files\n")

    iterations = 0
    diff_max = 0.0
    diff_o2 = 0
    d_min = list()
    d_max = list()
    avg_precision_order_1 = list()
    avg_precision_order_2 = list()
    for i, line in enumerate(f00):
        if line[0] == " " or line[0] == "-" or line[0]== "\033" or line[0]== "\n": continue
        if line[0] == "=":
            it = int(line[3:-4])
            print(f"Processing iteration {it}", end="\r")
            d_min.append(99999999999999)
            d_max.append(0.0)
            iterations = 0
            if(it != 0):
                if d_min[-1] == 99999999999999:
                    d_min[-1] = 0.0
                avg_precision_order_1.append(diff_max/(3*iterations))
                avg_precision_order_2.append(diff_o2/(3*iterations))
            continue
        line0 = line.split("|")[2][:-1]
        line1 = f11[i].split("|")[2][:-1]
        if line0 != line1:
            var0 = 0.0
            var1 = 0.0
            diff = 0
            o2 = 0
            for i in range(3):
                var0 = np.float64(line0.split("_")[i])
                var1 = np.float64(line1.split("_")[i])
                diff = abs(var1-var0)
                if (diff < d_min[-1] and diff != 0.0): d_min[-1] = diff
                if (diff > d_max[-1]): d_max[-1] = diff
                diff_max += diff
                o2 += np.power((var1-var0),2)
            diff_o2 += np.sqrt(o2)
        iterations += 1

    avg_precision_order_1.append(diff_max/(3*iterations))
    avg_precision_order_2.append(diff_o2/(3*iterations))
    if d_min[-1] == 99999999999999:
        d_min[-1] = 0.0

    print(f"Less precision difference :\t\t{np.min(d_min)}")
    print(f"Max precision difference :\t\t{np.max(d_max)}")
    print(f"Average precision loss at order 1 :\t{np.mean(avg_precision_order_1)}")
    print(f"Average precision loss at order 2 :\t{np.mean(avg_precision_order_2)}")


    if output_filename:
        # Write in seperate file
        fo = open(file_path.replace("/p_","/avg_p_"), "w")
        for i in range(len(d_min)):
            fo.write(f"{i} {d_min[i]} {d_max[i]} {avg_precision_order_1[i]} {avg_precision_order_2[i]}\n")
        fo.close
        # Write in one file
        fo = open(output_filename, "a")
        fo.write(f"{file_path.replace('.dat','').split('p_')[-1]} {np.min(d_min)} {np.max(d_max)} {np.mean(avg_precision_order_1)} {np.mean(avg_precision_order_2)}\n")
        fo.close()

# print(f"diff_max : {diff_max} - iterations : {iterations}")
# print("ratio :", diff_max/iterations)

# given precision : 32b long mantissa
# np.float32 : 8b long mantissa
# float = np.float64 : 17b long mantissa
# float = np.longdouble : 19b long mantissa

baseline_filename = "../0/out.dat"
f1_filename = "out.dat"

if len(argv) == 1:
    f0_filename = baseline_filename

if len(argv) == 2:
    f0_filename = baseline_filename
    f1_filename = argv[1]

if len(argv) == 3:
    f0_filename = argv[1]
    f1_filename = argv[2]

# if len(argv) == 4:
#     f0_filename = argv[1]
#     f1_filename = argv[2]
#     output_filename = argv[3]

if isdir(f1_filename):
    folder = f1_filename
    if folder[-1] == "/":
        folder = folder[:-1]
    output_filename = f"{folder}/avg_p.dat"
    for file in listdir(folder):
        if not file.startswith("p"): continue
        # output_filename = f"{folder}/avg_{file}"
        file_path = f"{folder}/{file}"
        compare(baseline_filename, file_path, output_filename)

else:
    compare(f0_filename, f1_filename)



