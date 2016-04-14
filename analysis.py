from clustering import CHalos
import os
import re
import time
import numpy as np
import shutil
import glob
import pylab as plt
from uncertainpy import prettyPlot
import scipy.ndimage.filters
import scipy.interpolate


b = 0.28
minNrHalos = 10
linkingLength = 2

nrBins = 5
sigma = 2

data_folder = "data"
output_dir = "results"



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



def allFOF(foldername, output_dir="results"):
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    os.makedirs(os.path.join(output_dir, "figures"))

    pattern = re.compile(r"^(\d{4}.*)(H|V)(\d_position.csv)$")

    start_time = time.time()


    sizes = {}
    nrParticlesInHalos = {}
    nrHalos = {}
    nrParticles = {}
    percentageInHalos = {}

    for root, dirs, files in os.walk(foldername):
        for filename in files:
            if filename.endswith(".csv"):
                result = pattern.search(filename)

                halos = CHalos(os.path.join(root, filename), b, minNrHalos)
                halos.linkingLength = linkingLength
                halos.FOF()

                # For testing purposes only
                # halos.nrParticlesInHalos = 0
                # halos.nrParticles = 0
                # halos.nrHalos = 0
                # halos.percentageInHalos = 0


                tmp = root.split("/")
                key = tmp[1] + "_" + tmp[2] + "_" + result.group(1) + "_" + result.group(2)

                # f = open(os.path.join(output_dir, tmp[1] + "_" + tmp[2] + "_" + result.group(1) + "_" + result.group(2) + ".dat"), "a")
                # f.write("{},{},{},{}\n".format(halos.nrParticles,
                #                                halos.nrHalos,
                #                                halos.nrParticlesInHalos,
                #                                halos.percentageInHalos))
                # f.close()

                if key in sizes:
                    sizes[key].append(halos.getSortedSizes())
                    nrParticlesInHalos[key].append(halos.nrParticlesInHalos)
                    nrHalos[key].append(halos.nrHalos)
                    nrParticles[key].append(halos.nrParticles)
                    percentageInHalos[key].append(halos.percentageInHalos)
                else:
                    sizes[key] = [halos.getSortedSizes()]
                    nrParticlesInHalos[key] = [halos.nrParticlesInHalos]
                    nrHalos[key] = [halos.nrHalos]
                    nrParticles[key] = [halos.nrParticles]
                    percentageInHalos[key] = [halos.percentageInHalos]



    for key in sizes:
        f = open(os.path.join(output_dir, key + ".txt"), "w")
        f.write("                    Mean  Var\n")
        f.write("nrParticles:        {} {}\n".format(np.mean(nrParticles[key]), np.var(nrParticles[key])))
        f.write("nrHalos:            {} {}\n".format(np.mean(nrHalos[key]), np.var(nrHalos[key])))
        f.write("nrParticlesInHalos: {} {}\n".format(np.mean(nrParticlesInHalos[key]), np.var(nrParticlesInHalos[key])))
        f.write("percentageInHalos:  {} {}\n".format(np.mean(percentageInHalos[key]), np.var(percentageInHalos[key])))
        f.close()




    # PLotting
    for key in sizes:

        mean = np.mean(sizes[key], 0)
        var = np.var(sizes[key], 0)
        # prettyPlot(np.mean(sizes[key], 0), xrange(1, len(sizes[key][0]) + 1), xlabel="nrParticles", ylabel="nr of cluster with number of particles > nrParticles", title="Cumulative cluster size", color=0)
        # prettyPlot(np.mean(sizes[key], 0), xlabel="nrParticles", ylabel="nr of cluster with number of particles > nrParticles", title="Cumulative cluster size", color=0)
        # prettyPlot(np.var(sizes[key], 0), xlabel="nrParticles", ylabel="nr of cluster with number of particles > nrParticles", title="Cumulative cluster size", color=2, new_figure=False)
        # plt.yscale('log')
        # prettyPlot(xrange(1, len(sizes[key][0]) + 1),mean, ylabel="nrParticles", xlabel="nr of cluster with number of particles > nrParticles", title="Cumulative cluster size", color=0)


        #


        color1 = 0
        color2 = 8

        ax, tableau20 = prettyPlot(xrange(1, len(mean) + 1), mean, "Cumulative cluster size", "nr of cluster with number of particles > nrParticles", "nrParticles, mean", color1)
        ax2 = ax.twinx()
        ax2.tick_params(axis="y", which="both", right="on", left="off", labelright="on",
                        color=tableau20[color2], labelcolor=tableau20[color2], labelsize=14)
        ax2.set_ylabel("nrParticles, variance", color=tableau20[color2], fontsize=16)
        ax.spines["right"].set_edgecolor(tableau20[color2])

        ax2.set_xlim([1, len(mean) + 1])
        ax2.set_ylim([min(var), max(var)])

        ax2.plot(xrange(1, len(var) + 1), var,
                 color=tableau20[color2], linewidth=2, antialiased=True)

        ax.tick_params(axis="y", color=tableau20[color1], labelcolor=tableau20[color1])
        ax.set_ylabel("nr of cluster with number of particles > nrParticles, mean", color=tableau20[color1], fontsize=16)
        ax.spines["left"].set_edgecolor(tableau20[color1])
        plt.savefig(os.path.join(output_dir, "figures", key + ".png"))





    # Calculate fractionalDifference between H and V
    pattern = re.compile(r"(.*)(H)$")

    for keyH in sizes:
        result = pattern.search(keyH)
        if result is not None:
            keyV = re.sub(pattern, r"\1V", keyH)
            if keyV in sizes:

                meanH = np.mean(sizes[keyH], 0)
                meanV = np.mean(sizes[keyV], 0)

                x, diff = fractionalDifference(meanH, meanV)

                prettyPlot(x, diff, ylabel="Fractional difference H/V", xlabel="nr of particles", title="Fractional difference", color=0)
                plt.savefig(os.path.join(output_dir, "figures", result.group(1) + "fractional.png"))



                prettyPlot(xrange(1, len(meanH) + 1), meanH, ylabel="nrParticles", xlabel="nr of cluster with number of particles > nrParticles", title="Cumulative cluster size", color=0)
                prettyPlot(xrange(1, len(meanV) + 1), meanV, ylabel="nrParticles", xlabel="nr of cluster with number of particles > nrParticles", title="Cumulative cluster size", color=2, new_figure=False)
                plt.ylim([min(min(meanH), min(meanV)), max(max(meanH), max(meanV))])
                plt.xlim([1, max(len(meanH), len(meanV))])
                plt.legend(["H", "V"])
                plt.savefig(os.path.join(output_dir, "figures", result.group(1) + "compare.png"))



    print "--- %s seconds ---" % (time.time() - start_time)


# TODO look closer at this
def fractionalDifference(meanH, meanV):
    histHdata = (plt.histogram(meanH, nrBins)[0])
    histVdata = (plt.histogram(meanV, nrBins)[0])

    x = (plt.histogram(meanH, nrBins)[1][1:])

    massH = np.zeros(nrBins)
    massV = np.zeros(nrBins)
    massH[-1] = histHdata[-1]
    massV[-1] = histVdata[-1]
    for i in range(nrBins-2, -1, -1):
        massH[i] = histHdata[i] + massH[i+1]
        massV[i] = histVdata[i] + massV[i+1]

    tmp = np.abs(massV - massH)/massH.astype(float)
    diff = (scipy.ndimage.filters.gaussian_filter(tmp, sigma))
    return x, diff


def compare(output_dir="results"):
    f = open(os.path.join(output_dir, "comparison.txt"))
    for filename in glob.glob(output_dir):
        data = np.loadtxt(filename, delimiter=",")

    f.close()



if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="convert all xls files to csv")
    # parser.add_argument("-f", "--full")


    allFOF(data_folder)
