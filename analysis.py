from clustering import CHalos
import os
import re
import time
import numpy as np
import shutil
import glob
import pylab as plt
from uncertainpy import prettyPlot, prettyBar
import scipy.ndimage.filters
import sys
import pickle
import argparse

b = 0.28
minNrHalos = 10
linkingLength = 1.5

nr_bins = 10
sigma = 1

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


def save_obj(obj, name, folder):
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
    save_obj(nrHalos, "nrHalosnrHalos", analysed_results_dir)
    save_obj(percentageInHalos, "percentageInHalos", analysed_results_dir)
    save_obj(nrParticles, "nrParticles", analysed_results_dir)

    print "--- %s seconds ---" % (time.time() - start_time)



def calculateMean(datasett):
    data_max = 0
    data_min = 100000
    for data in datasett:
        if max(data) > data_max:
            data_max = max(data)

        if min(data) < data_min:
            data_min = min(data)

    bins = np.linspace(data_min, data_max, nr_bins)

    tmp_mean = []
    for data in datasett:
        tmp_mean.append(np.histogram(data, bins=bins)[0])

    mean = np.mean(tmp_mean, 0)
    var = np.var(tmp_mean, 0)

    tmp = (bins[1:] - bins[:-1])/2.
    size = bins[1:] - tmp

    width = bins[1] - bins[0]

    return size, mean, var, bins, width





def plotAllSizes(sizes):
    # Plotting all data
    for key in sizes:
        c = 0
        i = 0

        t_max = 0
        t_min = 10000

        data_max = 0
        data_min = 10000

        legend = []
        for data in sizes[key]:
            t = range(1, len(data) + 1)
            prettyPlot(data, t, "nr cells", "nr of cluster with number of cells > nrParticles", "Cumulative cluster size", color=c, new_figure=False)

            if max(data) > data_max:
                data_max = max(data)
            if min(data) < data_min:
                data_min = min(data)
            if max(t) > t_max:
                t_max = max(t)
            if min(t) < t_min:
                t_min = min(t)

            legend.append("Datasett {}".format(i))
            i += 1
            c += 2

        plt.yscale('log')
        plt.xscale('log')
        plt.ylim([t_min, t_max])
        plt.xlim([data_min, data_max])
        plt.legend(legend)
        plt.savefig(os.path.join(output_dir, "figures", key + ".png"))
        plt.clf()


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


    # Plotting all data
    for key in sizes:
        c = 0
        i = 0

        t_max = 0
        t_min = 10000

        data_max = 0
        data_min = 10000

        legend = []
        for data in sizes[key]:
            t = range(1, len(data) + 1)
            prettyPlot(data, t, "Cumulative cluster size", "Cluster size", "Nr of clusters", color=c, new_figure=False)

            if max(data) > data_max:
                data_max = max(data)
            if min(data) < data_min:
                data_min = min(data)
            if max(t) > t_max:
                t_max = max(t)
            if min(t) < t_min:
                t_min = min(t)

            legend.append("Datasett {}".format(i))
            i += 1
            c += 2

        plt.yscale('log')
        plt.xscale('log')
        plt.ylim([t_min, t_max])
        plt.xlim([data_min, data_max])
        plt.legend(legend)

        plt.savefig(os.path.join(output_dir, "figures", key + ".png"))
        plt.clf()


    # sys.exit()
    for key in sizes:
        size, mean, var, bins, width = calculateMean(sizes[key])


        cumsum = np.cumsum(mean[::-1])[::-1]
        sumstd = np.sqrt(np.cumsum(var[::-1])[::-1])
        # prettyPlot(size, std, "Cumulative cluster size, mean", "nr cells, mean", "nr of cluster with number of cells > nrParticles", color=2)
        # prettyPlot(size, cumsum, "Cumulative cluster size, mean", "nr cells, mean", "nr of cluster with number of cells > nrParticles", color=0, new_figure=False)

        # plt.legend(["Mean", "Standard deviation"])

        colors = np.zeros(len(cumsum), dtype=int)
        ax, color = prettyBar(cumsum, index=bins[:-1], colors=colors, error=sumstd, width=width, linewidth=2)
        ax.set_xticks(bins-width/2)
        ax.set_xticklabels(np.round(bins, 0).astype(int), fontsize=14, rotation=0)

        plt.yscale('log')
        plt.ylabel("Nr of clusters")
        plt.xlabel("Cluster size", fontsize=16)
        plt.title("Cumulative cluster size, mean")
        # plt.xscale('log')
        # plt.show()
        plt.savefig(os.path.join(output_dir, "figures", key + "mean.png"))









    # Calculate fractionalDifference between H and V
    pattern = re.compile(r"(.*)(H)$")

    for keyH in sizes:
        result = pattern.search(keyH)
        if result is not None:
            keyV = re.sub(pattern, r"\1V", keyH)
            if keyV in sizes:


                data_max = 0
                data_min = 100000
                for data in sizes[keyH]:
                    if max(data) > data_max:
                        data_max = max(data)

                    if min(data) < data_min:
                        data_min = min(data)

                for data in sizes[keyV]:
                    if max(data) > data_max:
                        data_max = max(data)

                    if min(data) < data_min:
                        data_min = min(data)

                bins = np.linspace(data_min, data_max, nr_bins)


                tmp_meanH = []
                for data in sizes[keyH]:
                    tmp_meanH.append(np.histogram(data, bins=bins)[0])

                meanH = np.mean(tmp_meanH, 0)
                varH = np.var(tmp_meanH, 0)


                tmp_meanV = []
                for data in sizes[keyV]:
                    tmp_meanV.append(np.histogram(data, bins=bins)[0])

                meanV = np.mean(tmp_meanV, 0)
                varV = np.var(tmp_meanV, 0)


                tmp = (bins[1:] - bins[:-1])/2.
                size = bins[1:] - tmp


                cumsumH = np.cumsum(meanH[::-1])[::-1]
                cumsumV = np.cumsum(meanV[::-1])[::-1]
                sumvarV = np.sqrt(np.cumsum(varH[::-1])[::-1])
                sumvarH = np.sqrt(np.cumsum(varV[::-1])[::-1])
                sumstdV = np.sqrt(sumvarH)
                sumstdH = np.sqrt(sumvarV)

                # print sumH # nr of clusters
                # print sumV # nr of clusters
                # print t # Size of clusters

                diff = (cumsumV - cumsumH)/cumsumV


                # nonzeroV = np.where(cumsumV != 0)
                # nonzeroH = np.where(cumsumH != 0)
                # prettyPlot(size[nonzeroV], cumsumV[nonzeroV], "Cumulative cluster size, mean", "Cluster size", "Nr of cluster", color=0)
                # prettyPlot(size[nonzeroH], cumsumH[nonzeroH], "Cumulative cluster size, mean", "Cluster size", "Nr of cluster", color=2, new_figure=False)


                # # plt.plot(t, diff)
                # plt.legend(["V", "H"])
                # plt.yscale('log')
                # plt.xscale('log')
                # plt.savefig(os.path.join(output_dir, "figures", keyH[:-1] + "compare.png"))



                width = bins[1] - bins[0]
                colors = np.zeros(len(cumsumV), dtype=int)
                ax, color = prettyBar(cumsumV, index=bins[:-1], colors=colors, error=sumstdV, width=width, linewidth=2)
                ax, color = prettyBar(cumsumH, index=bins[:-1], colors=colors+4, error=sumstdH, width=width, linewidth=2, new_figure=False, alpha=0.6, error_kw=dict(ecolor=color[4], lw=2, capsize=10, capthick=2))
                ax.set_xticks(bins-width/2)
                ax.set_xticklabels(np.round(bins, 0).astype(int), fontsize=14, rotation=0)

                plt.yscale('log')
                plt.legend(["V", "H"])
                plt.ylabel("Nr of clusters")
                plt.xlabel("Cluster size", fontsize=16)
                plt.title("Cumulative cluster size, mean")

                # plt.xscale('log')
                # plt.show()
                plt.savefig(os.path.join(output_dir, "figures", keyH[:-1] + "compare.png"))


                width = bins[1] - bins[0]
                colors = np.zeros(len(cumsumV), dtype=int)
                ax, color = prettyBar(diff, index=bins[:-1], colors=colors, error=sumstdV, width=width, linewidth=2)
                ax.set_xticks(bins-width/2)
                ax.set_xticklabels(np.round(bins, 0).astype(int), fontsize=14, rotation=0)

                plt.yscale('log')
                plt.legend(["V", "H"])
                plt.ylabel("Nr of clusters")
                plt.xlabel("Cluster size", fontsize=16)
                plt.title("Cumulative cluster size, mean")

                prettyPlot(size, diff, "Fractional difference, (V-H)/V", "CLuster size", "Fractional difference nr of cluster", color=0)


                # plt.plot(t, diff)
                plt.savefig(os.path.join(output_dir, "figures", keyH[:-1] + "difference.png"))


                # sys.exit(1)

                #     massH = np.zeros(nrBins)
                #     massV = np.zeros(nrBins)
                #     massH[-1] = histHdata[-1]
                #     massV[-1] = histVdata[-1]
                #
                #     for i in range(nrBins-2, -1, -1):
                #         massH[i] = histHdata[i] + massH[i+1]
                #         massV[i] = histVdata[i] + massV[i+1]
                #
                #     tmp = np.abs(massV - massH)/massH.astype(float)


                #
                # meanH = np.mean(sizes[keyH], 0)
                # meanV = np.mean(sizes[keyV], 0)
                #
                # x, diff = fractionalDifference(meanH, meanV)
                #
                # prettyPlot(x, diff, ylabel="Fractional difference H/V", xlabel="nr of particles", title="Fractional difference", color=0)
                # plt.savefig(os.path.join(output_dir, "figures", result.group(1) + "fractional.png"))
                #
                #
                #
                # prettyPlot(xrange(1, len(meanH) + 1), meanH, ylabel="nrParticles", xlabel="nr of cluster with number of particles > nrParticles", title="Cumulative cluster size", color=0)
                # prettyPlot(xrange(1, len(meanV) + 1), meanV, ylabel="nrParticles", xlabel="nr of cluster with number of particles > nrParticles", title="Cumulative cluster size", color=2, new_figure=False)
                # plt.ylim([min(min(meanH), min(meanV)), max(max(meanH), max(meanV))])
                # plt.xlim([1, max(len(meanH), len(meanV))])
                # plt.legend(["H", "V"])
                # # plt.savefig(os.path.join(output_dir, "figures", result.group(1) + "compare.png"))
                # plt.show()
                #



#
# # TODO look closer at this
# def fractionalDifference(meanH, meanV):
#     histHdata = (plt.histogram(meanH, nrBins)[0])
#     histVdata = (plt.histogram(meanV, nrBins)[0])
#
#     x = (plt.histogram(meanH, nrBins)[1][1:])
#
#     massH = np.zeros(nrBins)
#     massV = np.zeros(nrBins)
#     massH[-1] = histHdata[-1]
#     massV[-1] = histVdata[-1]
#
#     for i in range(nrBins-2, -1, -1):
#         massH[i] = histHdata[i] + massH[i+1]
#         massV[i] = histVdata[i] + massV[i+1]
#
#     tmp = np.abs(massV - massH)/massH.astype(float)
#     diff = (scipy.ndimage.filters.gaussian_filter(tmp, sigma))
#
#     return x, diff
#




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
