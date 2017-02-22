from __future__ import division

from clustering import CHalos
from prettyplot import prettyPlot, prettyBar, get_colormap, spines_color, set_style, create_figure, get_current_colormap


import uncertainties.unumpy as unp

import os
import re
import time
import numpy as np
import shutil
import pylab as plt
import sys
import pickle
import argparse
import scipy
# import seaborne as sns

b = 0.28
minnrHalos = 10
linkingLength = 2

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
    with open(os.path.join(folder, name + ".pkl"), "wb") as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name, folder):
    with open(os.path.join(folder, name + ".pkl"), "rb") as f:
        return pickle.load(f)



def calculateLinkingLength(foldername):
    average_distance = []
    linking_lengths = []

    for root, dirs, files in os.walk(foldername):
        for filename in files:
            if filename.endswith(".csv"):
                halos = CHalos(os.path.join(root, filename), b, minnrHalos)

                average_distance.append(halos.average_distance)
                linking_lengths.append(halos.calculateLinkingLength())


    f = open(os.path.join("distances.txt"), "w")
    f.write("Average linking length: {}\n".format(sum(linking_lengths)/float(len(linking_lengths))))
    f.write("Average distance:       {}".format(sum(average_distance)/float(len(average_distance))))
    f.close()



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
    save_obj(nrHalos, "nrHalos", analysed_results_dir)
    save_obj(percentageInHalos, "percentageInHalos", analysed_results_dir)
    save_obj(nrParticles, "nrParticles", analysed_results_dir)

    print("--- {} seconds ---".format((time.time() - start_time)))









class Analysis:
    def __init__(self, output_dir="results", analysed_results_dir="obj"):
        self.analysed_results_dir = analysed_results_dir

        self.output_dir = os.path.join(output_dir, self.analysed_results_dir)
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.isdir(os.path.join(self.output_dir, "figures")):
            os.makedirs(os.path.join(self.output_dir, "figures"))

        self.load()


    def load(self):
        self.sizes = load_obj("sizes", self.analysed_results_dir)
        self.nrParticlesInHalos = load_obj("nrParticlesInHalos", self.analysed_results_dir)
        self.nrHalos = load_obj("nrHalos", self.analysed_results_dir)
        self.percentageInHalos = load_obj("percentageInHalos", self.analysed_results_dir)
        self.nrParticles = load_obj("nrParticles", self.analysed_results_dir)



    def _calculateMean(self, datasett):
        data_max = 0
        data_min = 100000
        for data in datasett:
            if max(data) > data_max:
                data_max = max(data)

            if min(data) < data_min:
                data_min = min(data)

        bins = np.linspace(data_min, data_max, nr_bins)

        mean = []
        for data in datasett:
            mean.append(np.histogram(data, bins=bins)[0])


        tmp = (bins[1:] - bins[:-1])/2.
        size = bins[1:] - tmp

        width = bins[1] - bins[0]

        return size, np.array(mean), bins, width




    def save_to_file(self):
        for key in self.sizes:
            f = open(os.path.join(self.output_dir, key + ".txt"), "w")
            f.write("                    Mean  stderror\n")
            f.write("nrParticles:        {} {}\n".format(np.mean(self.nrParticles[key]), scipy.stats.sem(self.nrParticles[key])))
            f.write("nrHalos:            {} {}\n".format(np.mean(self.nrHalos[key]), scipy.stats.sem(self.nrHalos[key])))
            f.write("self.nrParticlesInHalos: {} {}\n".format(np.mean(self.nrParticlesInHalos[key]), scipy.stats.sem(self.nrParticlesInHalos[key])))
            f.write("percentageInHalos:  {} {}\n".format(np.mean(self.percentageInHalos[key]), scipy.stats.sem(self.percentageInHalos[key])))
            f.close()



    def createInterpolation(self):
        self.interpolations = {}
        self.cumulative = {}
        self.unique_sizes = {}
        self.counts = {}

        max_size = []
        for key in self.sizes:
            for experiment in self.sizes[key]:
                max_size.append(max(experiment))
        self.max_size = max(max_size)



        for key in self.sizes:
            self.interpolations[key] = []
            self.cumulative[key] = []
            self.unique_sizes[key] = []
            self.counts[key] = []

            for experiment in self.sizes[key]:

                sizes, counts = np.unique(experiment, return_counts=True)
                #
                # print("=========================================")
                # extended_sizes = np.arange(sizes[-1] + 1, self.max_size + 6, 5)
                # extended_counts = np.zeros(len(extended_sizes))
                #
                # extended_sizes = np.concatenate((sizes, extended_sizes))
                # extended_counts = np.concatenate((counts, extended_counts))
                # # sizes = np.concatenate((sizes, extended_sizes))
                # # counts = np.concatenate((counts, extended_counts))

                cumulative = np.cumsum(counts[::-1])[::-1]
                #
                # print(extended_sizes)
                # print(extended_counts)
                # print sizes
                # print counts
                # # print cumulative
                # # print self.max_size
                # print cumulative

                interpolation = scipy.interpolate.InterpolatedUnivariateSpline(sizes, cumulative, k=1, ext=1)

                self.interpolations[key].append(interpolation)
                self.cumulative[key].append(cumulative)
                self.unique_sizes[key].append(sizes)
                self.counts[key].append(counts)


    def plotFractionalDifference(self):
        pattern = re.compile(r"(.*)(H)$")

        for keyH in self.sizes:
            result = pattern.search(keyH)
            if result is not None:
                keyV = re.sub(pattern, r"\1V", keyH)
                if keyV in self.sizes:


                    # Find min max
                    max_size = []
                    min_size = []
                    for size in self.unique_sizes[keyH]:
                        max_size.append(size.max())
                        min_size.append(size.min())

                    for size in self.unique_sizes[keyV]:
                        max_size.append(size.max())
                        min_size.append(size.min())


                    max_cumulative = []
                    min_cumulative = []
                    for cumulative in self.cumulative[keyV]:
                        max_cumulative.append(cumulative.max())
                        min_cumulative.append(cumulative.min())


                    for cumulative in self.cumulative[keyH]:
                        max_cumulative.append(cumulative.max())
                        min_cumulative.append(cumulative.min())


                    max_size = max(max_size)
                    min_size = min(min_size)

                    max_cumulative = max(max_cumulative)
                    min_cumulative = min(min_cumulative)

                    size = np.linspace(min_size, max_size, max_size - min_size + 1)


                    # Interpolate
                    cumulativeH = []
                    for inter in self.interpolations[keyH]:
                        cumulativeH.append(inter(size))

                    cumulativeH = np.array(cumulativeH)


                    cumulativeV = []
                    for inter in self.interpolations[keyV]:
                        cumulativeV.append(inter(size))

                    cumulativeV = np.array(cumulativeV)


                    # Calculate mean and error
                    meanH = np.mean(cumulativeH, axis=0)
                    stderrH = scipy.stats.sem(cumulativeH, 0)

                    meanV = np.mean(cumulativeV, axis=0)
                    stderrV = scipy.stats.sem(cumulativeV, 0)

                    # fractional_difference = abs(meanH - meanV)/meanV

                    meanHerr = unp.uarray(meanH, stderrH)
                    meanVerr = unp.uarray(meanV, stderrV)
                    # fractional_difference = abs(meanHerr - meanVerr)L/meanVerr

                    fractional_difference_value = np.zeros(len(size))
                    fractional_difference_err = np.zeros(len(size))

                    for i, (H, V) in enumerate(zip(meanHerr, meanVerr)):
                        if unp.nominal_values(H) == 0 or unp.nominal_values(V) == 0:
                            fractional_difference_value[i] = np.nan
                            fractional_difference_err[i] = np.nan
                        else:
                            tmp_fractional_difference = abs(H - V)/V
                            fractional_difference_value[i] = unp.nominal_values(tmp_fractional_difference)
                            fractional_difference_err[i] = unp.std_devs(tmp_fractional_difference)



                    mask = ~np.isnan(fractional_difference_value)

                    prettyPlot(size, fractional_difference_value,
                               "Fractional difference, $\\frac{{|H - V|}}{{V}}$\n {}".format(keyH[:-3]).replace("_", " "),
                               "Cluster size", "Nr of clusters",
                               style="seaborn-white")

                    colors = get_current_colormap()

                    plt.fill_between(size,
                                     fractional_difference_value - fractional_difference_err,
                                     fractional_difference_value + fractional_difference_err,
                                     alpha=0.5, color=colors[0])
                    plt.legend(["fractional difference", "standard error of the fractional difference"])



                    plt.savefig(os.path.join(self.output_dir, "figures", "fractional-difference_" + keyH[:-3] + ".png"))
                    plt.close()


    def plotCumulativeInterpolation(self):

        for key in self.sizes:

            legend = []
            nr_colors = len(self.unique_sizes[key])

            # Find min max
            max_size = []
            min_size = []
            for size in self.unique_sizes[key]:
                max_size.append(size.max())
                min_size.append(size.min())

            max_size = max(max_size)
            min_size = min(min_size)

            #
            # max_cumulative = []
            # min_cumulative = []


            # array of "integers" from min up to inlcuding max
            size = np.linspace(min_size, max_size, max_size - min_size + 1)

            for i, inter in enumerate(self.interpolations[key]):

                cumulative = inter(size)
                prettyPlot(size, cumulative,
                           "Cumulative cluster size, interpolated",
                           "Cluster size", "Nr of clusters",
                           color=i,
                           nr_colors=nr_colors,
                           new_figure=False,
                           style="seaborn-white")

                legend.append("Datasett {}".format(i + 1))

            #     max_cumulative.append(cumulative.max())
            #     min_cumulative.append(cumulative.min())
            #
            #
            # max_cumulative = max(max_cumulative)
            # min_cumulative = min(min_cumulative)
            #
            #
            # plt.xlim([min_size, max_size])
            # plt.ylim([min_cumulative, max_cumulative])

            plt.yscale("log")
            plt.xscale("log")

            plt.legend(legend)

            plt.savefig(os.path.join(self.output_dir, "figures", "cumulative_" + key + "_interpolated.png"))
            plt.close()


    def plotGrid(self):

        nr_plots = len(self.sizes.keys())/2
        grid_size = np.ceil(np.sqrt(nr_plots))
        grid_x_size = int(grid_size)
        grid_y_size = int(np.ceil(nr_plots/float(grid_x_size)))

        nr_colors = 2
        set_style(style="seaborn-white", nr_colors=nr_colors)
        fig, axes = plt.subplots(nrows=grid_y_size, ncols=grid_x_size)

        # Add a larger subplot to use to set a common xlabel and ylabel
        j = 0

        # Calculate fractionalDifference between H and V
        pattern = re.compile(r"(.*)(H)$")


        for keyH in self.sizes:
            result = pattern.search(keyH)
            if result is not None:
                keyV = re.sub(pattern, r"\1V", keyH)
                if keyV in self.sizes:


                    # Find min max
                    max_size = []
                    min_size = []
                    for size in self.unique_sizes[keyH]:
                        max_size.append(size.max())
                        min_size.append(size.min())

                    for size in self.unique_sizes[keyV]:
                        max_size.append(size.max())
                        min_size.append(size.min())


                    max_cumulative = []
                    min_cumulative = []
                    for cumulative in self.cumulative[keyV]:
                        max_cumulative.append(cumulative.max())
                        min_cumulative.append(cumulative.min())


                    max_size = max(max_size)
                    min_size = min(min_size)

                    max_cumulative = max(max_cumulative)
                    min_cumulative = min(min_cumulative)


                    nx = j % grid_x_size
                    ny = int(np.floor(j/float(grid_x_size)))

                    if grid_y_size == 1:
                        ax = axes[nx]
                    else:
                        ax = axes[ny][nx]


                    for size, cumulative in zip(self.unique_sizes[keyV], self.cumulative[keyV]):
                        prettyPlot(size, cumulative, color=0, nr_colors=nr_colors, ax=ax)

                    for size, cumulative in zip(self.unique_sizes[keyH], self.cumulative[keyH]):
                        prettyPlot(size, cumulative, color=1, nr_colors=nr_colors, ax=ax)
                    #
                    tmp_title = keyH[:-3].split("_")
                    title = tmp_title[0] + " " + tmp_title[1] + "\n" + tmp_title[2]
                    ax.set_title(title, fontsize=10)

                    ax.set_yscale("log")
                    ax.set_xscale("log")

                    ax.set_xlim([min_size, max_size])
                    ax.set_ylim([min_cumulative, max_cumulative])


                    for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_rotation(-25)

                    for tick in ax.xaxis.get_minor_ticks():
                        tick.label.set_rotation(-25)

                    # ax.grid(False)
                    # ax.set_facecolor("white")


                    j += 1



        for i in xrange(j, grid_x_size*grid_y_size):
            nx = i % grid_x_size
            ny = int(np.floor(i/float(grid_x_size)))

            if grid_y_size == 1:
                ax = axes[nx]
            else:
                ax = axes[ny][nx]

            ax.axis("off")



        colormap = get_current_colormap()

        lineV = ax.plot(min_size-1, min_cumulative-1, color=colormap[0])[0]
        lineH = ax.plot(min_size-1, min_cumulative-1, color=colormap[1])[0]

        plt.figlegend([lineV, lineH], ["V", "H"], (0.75, 0.2))

        plt.suptitle("Cumulative cluster size", fontsize=18)
        plt.tight_layout()
        plt.subplots_adjust(top=0.9, bottom=0.12, left=0.1, right=0.9)
        fig.text(0.5, 0.02, "Cumulative cluster size", ha="center", size=16)
        fig.text(0.02, 0.5, "Cluster size", va="center", rotation="vertical", size=16)
        plt.savefig(os.path.join(self.output_dir, "figures", "cumulative_grid.png"))
        plt.close()



    def plotCompare(self):

        pattern = re.compile(r"(.*)(H)$")

        for keyH in self.sizes:
            result = pattern.search(keyH)
            if result is not None:
                keyV = re.sub(pattern, r"\1V", keyH)
                if keyV in self.sizes:


                    # Find min max
                    max_size = []
                    min_size = []
                    for size in self.unique_sizes[keyH]:
                        max_size.append(size.max())
                        min_size.append(size.min())

                    for size in self.unique_sizes[keyV]:
                        max_size.append(size.max())
                        min_size.append(size.min())


                    max_cumulative = []
                    min_cumulative = []
                    for cumulative in self.cumulative[keyV]:
                        max_cumulative.append(cumulative.max())
                        min_cumulative.append(cumulative.min())


                    max_size = max(max_size)
                    min_size = min(min_size)

                    max_cumulative = max(max_cumulative)
                    min_cumulative = min(min_cumulative)


                    nr_colors = 2
                    create_figure(style="seaborn-white", nr_colors=nr_colors)
                    prettyPlot([min_size-1], [min_cumulative-1], "Cumulative cluster size", "Cluster size", "Nr of clusters", color=0, nr_colors=nr_colors, new_figure=False)
                    prettyPlot([min_size-1], [min_cumulative-1], "Cumulative cluster size", "Cluster size", "Nr of clusters", color=1, nr_colors=nr_colors, new_figure=False)
                    plt.legend(["V", "H"])



                    for size, cumulative in zip(self.unique_sizes[keyV], self.cumulative[keyV]):
                        prettyPlot(size, cumulative, "Cumulative cluster size", "Cluster size", "Nr of clusters", color=0, nr_colors=nr_colors, new_figure=False)

                    for size, cumulative in zip(self.unique_sizes[keyH], self.cumulative[keyH]):
                        prettyPlot(size, cumulative, "Cumulative cluster size", "Cluster size", "Nr of clusters", color=1, nr_colors=nr_colors, new_figure=False)


                    plt.xlim([min_size, max_size])
                    plt.ylim([min_cumulative, max_cumulative])

                    plt.yscale("log")
                    plt.xscale("log")

                    plt.savefig(os.path.join(self.output_dir, "figures", keyH[:-1] + "combined.png"))
                    plt.close()

    def plotMean(self):
        for key in self.sizes:

            # Find min max
            max_size = []
            min_size = []
            for size in self.unique_sizes[key]:
                max_size.append(size.max())
                min_size.append(size.min())

            max_size = max(max_size)
            min_size = min(min_size)



            # array of "integers" from min up to inlcuding max
            size = np.linspace(min_size, max_size, max_size - min_size + 1)

            cumulative = []

            for inter in self.interpolations[key]:
                cumulative.append(inter(size))

            cumulative = np.array(cumulative)

            mean = np.mean(cumulative, axis=0)
            p_05 = np.percentile(cumulative, 5, 0)
            p_95 = np.percentile(cumulative, 95, 0)
            stderr = scipy.stats.sem(cumulative, 0)

            prettyPlot(size, mean,
                       "Mean cluster size, {}".format(key).replace("_", " "),
                       "Cluster size", "Nr of clusters",
                       style="seaborn-white")

            colors = get_current_colormap()
            #
            # plt.fill_between(size, p_05, p_95,
            #                  alpha=0.5, color=colors[0])\
            # plt.legend(["mean", "90\% confidence interval"])

            plt.fill_between(size, mean - stderr, mean + stderr,
                             alpha=0.5, color=colors[0])
            plt.legend(["mean", "standard error of the mean"])


            # plt.xlim([min_size, max_size])
            plt.ylim([min(mean), max(mean)])
            #
            plt.yscale("log")
            plt.xscale("log")


            plt.savefig(os.path.join(self.output_dir, "figures", "mean_" + key + ".png"))
            plt.close()



    def plotCumulative(self):
        for key in self.unique_sizes:

            legend = []
            nr_colors = len(self.unique_sizes[key])


            # Find min max
            max_size = []
            min_size = []
            max_cumulative = []
            min_cumulative = []

            i = 1
            create_figure(style="seaborn-white", nr_colors=nr_colors)
            for size, cumulative in zip(self.unique_sizes[key], self.cumulative[key]):

                prettyPlot(size, cumulative,
                           "Cumulative cluster size",
                           "Cluster size",
                           "Nr of clusters",
                           nr_colors=nr_colors,
                           new_figure=False,
                           style="seaborn-white")


                legend.append("Datasett {}".format(i))

                max_size.append(size.max())
                min_size.append(size.min())

                max_cumulative.append(cumulative.max())
                min_cumulative.append(cumulative.min())

                i+=1

            max_size = max(max_size)
            min_size = min(min_size)
            max_cumulative = max(max_cumulative)
            min_cumulative = min(min_cumulative)


            plt.xlim([min_size, max_size])
            plt.ylim([min_cumulative, max_cumulative])
            plt.yscale("log")
            plt.xscale("log")
            plt.legend(legend)

            plt.savefig(os.path.join(self.output_dir, "figures", "cumulative_" + key + ".png"))
            plt.close()


    def results_new(self):
        self.createInterpolation()

        self.save_to_file()
        self.plotCumulative()
        self.plotCumulativeInterpolation()
        self.plotMean()
        self.plotGrid()
        self.plotCompare()
        self.plotFractionalDifference()







if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform a clustering analysis")
    parser.add_argument("analysed_results_dir", help="Analysed results folder")
    parser.add_argument("-f", "--fof", action="store_true", help="Perform the friend of friends analysis")
    parser.add_argument("-r", "--results", action="store_true", help="Analyse the results")
    parser.add_argument("-l", "--linking_lengths", action="store_true", help="Calculate the linking length")

    args = parser.parse_args()

    if args.fof:
        allFOF(data_folder, analysed_results_dir=args.analysed_results_dir)

    if args.results:
        analysis = Analysis(analysed_results_dir=args.analysed_results_dir)
        analysis.results()


    if args.linking_lengths:
        calculateLinkingLength(data_folder)
