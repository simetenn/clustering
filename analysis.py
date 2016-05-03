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
import sys
import pickle
import argparse

b = 0.28
minNrHalos = 10
linkingLength = 1.5

nrBins = 5
sigma = 2

data_folder = "data"
output_dir = "results"


sys.setrecursionlimit(50000)

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


def save_obj(obj, folder, name):
    with open(os.path.join(folder, name + '.pkl'), 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name, folder):
    with open(os.path.join(folder, name + '.pkl'), 'rb') as f:
        return pickle.load(f)


def allFOF(foldername, analysed_results_dir="obj"):

    if not os.path.isdir(analysed_results_dir):
        # shutil.rmtree(analysed_results_dir)
        os.makedirs(analysed_results_dir)

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


    save_obj(sizes, "sizes", analysed_results_dir)
    save_obj(nrParticlesInHalos, "nrParticlesInHalos", analysed_results_dir)
    save_obj(nrHalos, "nrHalosnrHalos")
    save_obj(percentageInHalos, "percentageInHalos", analysed_results_dir)
    save_obj(nrParticles, "nrParticles", analysed_results_dir)

    print "--- %s seconds ---" % (time.time() - start_time)


def calculateMean(datasett):
    t_max = []
    interpolations = []
    legend = []
    i = 0
    c = 0
    # plt.figure()
    for data in datasett:
        t = xrange(1, len(data) + 1)

        if len(t) > len(t_max):
            t_max = t

        interpolations.append(scipy.interpolate.InterpolatedUnivariateSpline(t, data, k=3))

        # prettyPlot(t, data, new_figure=False, color=c)
        legend.append("{}".format(i))
        i += 1
        c += 2


    # plt.show()
    tmp_mean = []
    for interpolation in interpolations:
        tmp_mean.append(interpolation(t_max))

    # Check if interpolation is equal original values
    for i in range(len(tmp_mean)):
        print range(1, len(datasett[i]) + 1)[:10]
        print datasett[i][:10]
        plt.loglog(datasett[i], xrange(1, len(datasett[i]) + 1), "*")
        plt.loglog(tmp_mean[i], t_max)

        plt.show()

    mean = np.mean(tmp_mean, 0)
    std = np.std(tmp_mean, 0)

    return mean, std



def results(output_dir="results", analysed_results_dir="obj"):
    output_dir = os.path.join(output_dir, analysed_results_dir)
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    os.makedirs(os.path.join(output_dir, "figures"))


    sizes = load_obj("sizes", analysed_results_dir)
    nrParticlesInHalos = load_obj("nrParticlesInHalos", analysed_results_dir)
    nrHalos = load_obj("nrHalosnrHalos", analysed_results_dir)
    percentageInHalos = load_obj("percentageInHalos", analysed_results_dir)
    nrParticles = load_obj("nrParticles", analysed_results_dir)


    for key in sizes:
        f = open(os.path.join(output_dir, key + ".txt"), "w")
        f.write("                    Mean  std\n")
        f.write("nrParticles:        {} {}\n".format(np.mean(nrParticles[key]), np.std(nrParticles[key])))
        f.write("nrHalos:            {} {}\n".format(np.mean(nrHalos[key]), np.std(nrHalos[key])))
        f.write("nrParticlesInHalos: {} {}\n".format(np.mean(nrParticlesInHalos[key]), np.std(nrParticlesInHalos[key])))
        f.write("percentageInHalos:  {} {}\n".format(np.mean(percentageInHalos[key]), np.std(percentageInHalos[key])))
        f.close()




    # PLotting
    for key in sizes:
        mean, std = calculateMean(sizes[key])

        color1 = 0
        color2 = 8

        # ax, tableau20 = prettyPlot(xrange(1, len(mean) + 1), mean, "Cumulative cluster size", "nr of cluster with number of cells > nrParticles", "nr cells, mean", color1)


        # prettyPlot(xrange(1, len(mean) + 1), std, "Cumulative cluster size", "nr of cluster with number of particles > nrParticles", "nrParticles, mean", color2, new_figure=False)

        # plt.yscale('log')
        # ax2 = ax.twinx()
        # ax2.tick_params(axis="y", which="both", right="on", left="off", labelright="on",
        #                 color=tableau20[color2], labelcolor=tableau20[color2], labelsize=14)
        # ax2.set_ylabel("nrParticles, variance", color=tableau20[color2], fontsize=16)
        # ax.spines["right"].set_edgecolor(tableau20[color2])
        #
        # ax2.set_xlim([1, len(mean) + 1])
        # ax2.set_ylim([min(std), max(std)])
        #
        # ax2.plot(xrange(1, len(std) + 1), std,
        #          color=tableau20[color2], linewidth=2, antialiased=True)
        #
        # ax.tick_params(axis="y", color=tableau20[color1], labelcolor=tableau20[color1])
        # ax.set_ylabel("nr of cluster with number of particles > nrParticles, std", color=tableau20[color1], fontsize=16)
        # ax.spines["left"].set_edgecolor(tableau20[color1])
        #
        # plt.yscale('log')
        # plt.xscale('log')
        plt.savefig(os.path.join(output_dir, "figures", key + ".png"))










    sys.exit()



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





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform a clustering analysis")
    parser.add_argument("analysed_results_dir", help="Analysed results folder")
    parser.add_argument("-f", "--fof", action="store_true", help="Perform the friend of friends analysis")
    parser.add_argument("-r", "--results", action="store_true", help="Analyse the results")

    #crash at data/eksitatorisk/visuell/1500 V1 WFA+PSD95 V6_position.csv

    args = parser.parse_args()

    if args.fof:
        allFOF(data_folder, analysed_results_dir=args.analysed_results_dir)

    if args.results:
        results(analysed_results_dir=args.analysed_results_dir)
