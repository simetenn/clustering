from clustering import CHalos
import argparse
import os
import re
import time
import numpy as np

b = 0.28
minNrHalos = 10
linkingLength = 2

def pairFiles(foldername):
    pattern = re.compile(r"(.*)H(\d_position.xls)$")

    file_pairs = []

    for root, dirs, files in os.walk(foldername):

        for filename in files:
            result = pattern.search(filename)
            if result is not None:
                filename2 = re.sub(pattern, r"\1V\2", filename)
                if filename2 in files:

                    file_pairs.append((root, os.path.join(root, filename), os.path.join(root, filename2)))

    return file_pairs



def compare(filename1, filename2, outdir):
    tmp = outdir.split("/")



    halos1 = CHalos(filename1, b, minNrHalos)
    halos1.linkingLength = linkingLength
    halos1.FOF()

    halos1.printInformation()

    halos2 = CHalos(filename2, b, minNrHalos)
    halos2.linkingLength = linkingLength
    halos2.FOF()

    halos2.printInformation()

def analysis(foldername):
    start_time = time.time()

    file_pairs = pairFiles(foldername)
    compare(file_pairs[0][1], file_pairs[0][2], file_pairs[0][0])


    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert all xls files to csv")
    parser.add_argument("foldername")

    args = parser.parse_args()

    analysis(args.foldername)
