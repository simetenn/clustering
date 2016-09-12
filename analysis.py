from clustering import CHalos
from uncertainpy import prettyPlot, prettyBar, get_colormap, spines_edge_color, set_style

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

b = 0.28
minNrHalos = 10
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
    with open(os.path.join(folder, name + '.pkl'), 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name, folder):
    with open(os.path.join(folder, name + '.pkl'), 'rb') as f:
        return pickle.load(f)



def calculateLinkingLength(foldername):
    average_distance = []
    linking_lengths = []

    for root, dirs, files in os.walk(foldername):
        for filename in files:
            if filename.endswith(".csv"):
                halos = CHalos(os.path.join(root, filename), b, minNrHalos)

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

    mean = []
    for data in datasett:
        mean.append(np.histogram(data, bins=bins)[0])


    tmp = (bins[1:] - bins[:-1])/2.
    size = bins[1:] - tmp

    width = bins[1] - bins[0]

    return size, np.array(mean), bins, width




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
        f.write("                    Mean  stderror\n")
        f.write("nrParticles:        {} {}\n".format(np.mean(nrParticles[key]), scipy.stats.sem(nrParticles[key])))
        f.write("nrHalos:            {} {}\n".format(np.mean(nrHalos[key]), scipy.stats.sem(nrHalos[key])))
        f.write("nrParticlesInHalos: {} {}\n".format(np.mean(nrParticlesInHalos[key]), scipy.stats.sem(nrParticlesInHalos[key])))
        f.write("percentageInHalos:  {} {}\n".format(np.mean(percentageInHalos[key]), scipy.stats.sem(percentageInHalos[key])))
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
        nr_hues = len(sizes[key])
        for data in sizes[key]:
            t = range(1, len(data) + 1)

            prettyPlot(data, t, "Cumulative cluster size", "Cluster size", "Nr of clusters", color=c, nr_hues=nr_hues, new_figure=False)

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
            c += 1

        plt.yscale('log')
        plt.xscale('log')
        plt.ylim([t_min, t_max])
        plt.xlim([data_min, data_max])
        plt.legend(legend)

        plt.savefig(os.path.join(output_dir, "figures", key + ".png"))
        # plt.clf()
        plt.close()


    #Calculate and plot mean values
    for key in sizes:
        size, mean, bins, width = calculateMean(sizes[key])

        cumsumAll = np.cumsum(mean[:, ::-1], 1)[:, ::-1]
        cumsum = np.mean(cumsumAll, 0)
        sumstd = scipy.stats.sem(cumsumAll, 0)
        # prettyPlot(size, std, "Cumulative cluster size, mean", "nr cells, mean", "nr of cluster with number of cells > nrParticles", color=2)
        # prettyPlot(size, cumsum, "Cumulative cluster size, mean", "nr cells, mean", "nr of cluster with number of cells > nrParticles", color=0, new_figure=False)

        # plt.legend(["Mean", "Standard deviation"])
        max_sizes = []
        for size in sizes[key]:
            max_sizes.append(max(size))
        max_size = max(max_sizes)

        bins = (bins/max_size)*100
        width = (width/max_size)*100

        ax = prettyBar(cumsum, index=bins[:-1], color=0, nr_hues=2, error=sumstd, width=width, linewidth=2)
        ax.set_xticks(bins-width/2)
        ax.set_xticklabels(np.round(bins, 0).astype(int), fontsize=14, rotation=0)

        plt.yscale('log')
        plt.ylabel("Nr of clusters")
        plt.xlabel("Cluster size [\% of max cluster size {}]".format(max_size), fontsize=16)
        plt.title("Cumulative cluster size, mean")

        plt.savefig(os.path.join(output_dir, "figures", key + "mean.png"))
        plt.close()



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





                # plot cumulative cluster size H and V from same animal in the same plot
                t_max = 0
                t_min = 10000
                nr_hues = 2
                prettyPlot([data_min-1], [t_min-1], "Cumulative cluster size", "Cluster size", "Nr of clusters", color=0, nr_hues=nr_hues, new_figure=False)
                prettyPlot([data_min-1], [t_min-1], "Cumulative cluster size", "Cluster size", "Nr of clusters", color=1, nr_hues=nr_hues, new_figure=False)
                plt.legend(["V", "H"])


                for data in sizes[keyV]:
                    t = range(1, len(data) + 1)

                    prettyPlot(data, t, "Cumulative cluster size", "Cluster size", "Nr of clusters", color=0, nr_hues=nr_hues, new_figure=False)

                    if max(t) > t_max:
                        t_max = max(t)
                    if min(t) < t_min:
                        t_min = min(t)

                    i += 1

                for data in sizes[keyH]:
                    t = range(1, len(data) + 1)

                    prettyPlot(data, t, "Cumulative cluster size", "Cluster size", "Nr of clusters", color=1, nr_hues=nr_hues, new_figure=False)

                    if max(t) > t_max:
                        t_max = max(t)
                    if min(t) < t_min:
                        t_min = min(t)

                    i += 1

                plt.yscale('log')
                plt.xscale('log')
                plt.ylim([t_min, t_max])
                plt.xlim([data_min, data_max])


                plt.savefig(os.path.join(output_dir, "figures", keyH[:-1] + "combined.png"))
                plt.close()







                meanH = []
                for data in sizes[keyH]:
                    meanH.append(np.histogram(data, bins=bins)[0])

                meanH = np.array(meanH)


                meanV = []
                for data in sizes[keyV]:
                    meanV.append(np.histogram(data, bins=bins)[0])

                meanV = np.array(meanV)

                width = bins[1] - bins[0]
                cumsumAllH = np.cumsum(meanH[:, ::-1], 1)[:, ::-1]
                cumsumAllV = np.cumsum(meanV[:, ::-1], 1)[:, ::-1]


                cumsumH = np.mean(cumsumAllH, 0)
                cumsumV = np.mean(cumsumAllV, 0)
                sumstdV = np.std(cumsumAllV, 0)
                sumstdH = np.std(cumsumAllH, 0)


                diff = 1 - cumsumH/cumsumV.astype(float)

                stddiff = diff*np.sqrt(sumstdV**2/cumsumV**2 + sumstdH**2/cumsumH**2)

                #plot two mean values against each other
                bins = (bins/data_max)*100

                width = bins[1] - bins[0]
                ax = prettyBar(cumsumV, index=bins[:-1], color=0, nr_hues=2, error=sumstdV, width=width, linewidth=2)
                ax = prettyBar(cumsumH, index=bins[:-1], color=4, nr_hues=6, error=sumstdH, width=width, linewidth=2,
                               new_figure=False, alpha=0.6,
                               error_kw=dict(ecolor=get_colormap()[4], lw=2, capsize=10, capthick=2))
                ax.set_xticks(bins-width/2)
                ax.set_xticklabels(np.round(bins, 0).astype(int), fontsize=14, rotation=0)

                plt.yscale('log')
                plt.legend(["V", "H"])
                plt.ylabel("Nr of clusters")
                plt.xlabel("Cluster size [\% of max cluster size {}]".format(data_max), fontsize=16)
                plt.title("Cumulative cluster size, mean")

                plt.savefig(os.path.join(output_dir, "figures", keyH[:-1] + "compare.png"))
                plt.close()




                # Plot fractional difference
                width = bins[1] - bins[0]
                ax = prettyBar(diff, index=bins[:-1], color=0, nr_hues=2, error=stddiff, width=width, linewidth=2)
                ax.set_xticks(bins-width/2)
                ax.set_xticklabels(np.round(bins, 0).astype(int), fontsize=14, rotation=0)

                plt.ylabel("Fractional difference nr of cluster")
                plt.xlabel("Cluster size [\% of max cluster size {}]".format(data_max), fontsize=16)
                plt.title("Fractional difference, (V-H)/V")

                # prettyPlot(size, diff, "Fractional difference, (V-H)/V", "CLuster size", "Fractional difference nr of cluster", color=0)

                plt.savefig(os.path.join(output_dir, "figures", keyH[:-1] + "difference.png"))
                plt.close()






def plot_grid(output_dir="results", analysed_results_dir="obj"):
    output_dir = os.path.join(output_dir, analysed_results_dir)
    if not os.path.isdir(os.path.join(output_dir, "figures")):
        os.makedirs(os.path.join(output_dir, "figures"))


    sizes = load_obj("sizes", analysed_results_dir)



    nr_plots = len(sizes.keys())/2
    grid_size = np.ceil(np.sqrt(nr_plots))
    grid_x_size = int(grid_size)
    grid_y_size = int(np.ceil(nr_plots/float(grid_x_size)))

    fig, axes = plt.subplots(nrows=grid_y_size, ncols=grid_x_size)

    # Add a larger subplot to use to set a common xlabel and ylabel
    set_style("white")
    ax = fig.add_subplot(111, zorder=-10)
    spines_edge_color(ax, edges={"top": "None", "bottom": "None",
                                 "right": "None", "left": "None"})
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax.set_xlabel("Cumulative cluster size")
    ax.set_ylabel("Cluster size")
    j = 0



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


                nx = j % grid_x_size
                ny = int(np.floor(j/float(grid_x_size)))

                if grid_y_size == 1:
                    ax = axes[nx]
                else:
                    ax = axes[ny][nx]


                t_max = 0
                t_min = 10000
                nr_hues = 2
                prettyPlot([data_min-1], [t_min-1], color=0, nr_hues=nr_hues, new_figure=False, ax=ax)
                prettyPlot([data_min-1], [t_min-1], color=1, nr_hues=nr_hues, new_figure=False, ax=ax)
                plt.legend(["V", "H"])


                for data in sizes[keyV]:
                    t = range(1, len(data) + 1)

                    prettyPlot(data, t, color=0, nr_hues=nr_hues, ax=ax)

                    if max(t) > t_max:
                        t_max = max(t)
                    if min(t) < t_min:
                        t_min = min(t)

                for data in sizes[keyH]:
                    t = range(1, len(data) + 1)

                    prettyPlot(data, t, color=1, nr_hues=nr_hues, ax=ax)

                    if max(t) > t_max:
                        t_max = max(t)
                    if min(t) < t_min:
                        t_min = min(t)

                tmp_title = keyH[:-3].split("_")
                title = tmp_title[0] + " " +  tmp_title[1] + "\n" + tmp_title[2]
                ax.set_title(title, fontsize=10)
                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_ylim([t_min, t_max])
                ax.set_xlim([data_min, data_max])

                for tick in ax.get_xticklabels():
                    tick.set_rotation(-30)

                ax.tick_params(labelsize=10)

                j += 1

    for i in xrange(j, grid_x_size*grid_y_size):
        nx = i % grid_x_size
        ny = int(np.floor(i/float(grid_x_size)))

        if grid_y_size == 1:
            ax = axes[nx]
        else:
            ax = axes[ny][nx]

        ax.axis("off")


    plt.suptitle("Cumulative cluster size", fontsize=18)
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.savefig(os.path.join(output_dir, "figures", "cumulative_grid.png"))
    plt.close()




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
        results(analysed_results_dir=args.analysed_results_dir)
        plot_grid(analysed_results_dir=args.analysed_results_dir)

    if args.linking_lengths:
        calculateLinkingLength(data_folder)
