import matplotlib
matplotlib.use('GTKAgg')

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import argparse
import time
import sys

from mpl_toolkits.mplot3d import Axes3D

sys.setrecursionlimit(50000)

def system(cmd, verbose=1):
    """Run system command cmd using subprocess module."""
    if verbose >= 1:
        print cmd

    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        if output:
            if verbose >= 2:
                print output


    except subprocess.CalledProcessError as e:
        if verbose >= 1:
            print """Command: \n%s \nfailed""" % cmd
            print 'Return code:', e.returncode
            print e.output


class CHalo:
    def __init__(self, ind=None, ID=-1):
        if ind is None:
            self.indices = []
            self.size = 0
        else:
            self.indices = ind
            self.size = len(ind)
            self.ID = ID

    def __add__(self, other):
        if type(other) == type(self):
            self.indices.extend(other.indices)
        else:
            self.indices.append(other)
        self.size += other.size

    def __radd__(self, other):
        if type(other) == type(self):
            self.indices.extend(other.indices)
        else:
            self.indices.append(other)
            self.size += other.size

    def __str__(self):
        return "".join(str(i) + " " for i in self.indices)



class CHalos:
    def __init__(self, name, b=0.28, minNrHalos=5, f=0.7):
        self.name = name

        self.halos = []
        self.nrHalos = len(self.halos)
        self.minNrHalos = minNrHalos
        self.b = b
        self.f = f
        self.nrLinking = 5000
        self.average_distance = 0

        self.loadPositions()

        self.nrParticlesInHalos = 0
        self.percentageInHalos = 0

        self.nrLinking = self.nrParticles

    def initialize(self):
        def limits(i):
            return (np.min(self.positions[:, i]), np.max(self.positions[:, i]))

        #Simulation details
        self.xLimits = limits(0)
        self.yLimits = limits(1)
        self.zLimits = limits(2)


        self.diffX = self.xLimits[1] - self.xLimits[0]
        self.diffY = self.yLimits[1] - self.yLimits[0]
        self.diffZ = self.zLimits[1] - self.zLimits[0]

        self.volume = self.diffX*self.diffY*self.diffZ
        self.nrParticles = self.positions.shape[0]

        self.average_distance = pow(self.volume/self.nrParticles, 1./3)
        self.linkingLength = self.b*pow(self.volume/self.nrParticles, 1./3)



    def loadPositions(self):

        print "Loading: " + self.name
        if self.name.split(".")[-1] == "xls":
            system("libreoffice --headless --convert-to csv " + self.name.replace(" ", "\\ "), 0)

        self.positions = np.loadtxt(self.name[:-3] + "csv", delimiter=",", skiprows=2, dtype="string")
        self.positions[:, 3] = "-1"
        self.positions = self.positions[:, :4].astype(float, copy=False)

        self.initialize()



    def distance2(self, particle1, particle2):
        return np.sum((self.positions[particle2, :3] - self.positions[particle1, :3])**2)



    def distance(self, particle1, particle2):
        return np.sqrt(np.sum((self.positions[particle2, :3] - self.positions[particle1, :3])**2))



    def findNeighbors(self, particle, allFriends):
        allFriends.append(particle)
        friends = []

        #Write this using numpy
        for nextParticle in self.nextParticles:
            if self.distance(particle, nextParticle) < self.linkingLength:
                    friends.append(nextParticle)

        for friend in friends:
            self.nextParticles.remove(friend)  # Might slow down


        for friend in friends:
            self.findNeighbors(friend, allFriends)


    def getSortedSizes(self):
        sizes = []
        for halo in self.halos:
            sizes.append(halo.size)

        # sizes = [10, 213, 22, 123, 10, 432]
        # sizes = np.array(sizes)*np.random.normal(1, 1)

        return np.sort(sizes)[::-1]


    def FOF(self):
        """
        Finding halos with the FOF-method
        """
        print "Calculating FOF"

        self.nextParticles = range(0, self.positions.shape[0])
        while (len(self.nextParticles) > 0):
            allFriends = []

            nextParticle = self.nextParticles[0]
            del self.nextParticles[0]
            self.findNeighbors(nextParticle, allFriends)

            if len(allFriends) >= self.minNrHalos:
                self.halos.append(CHalo(allFriends, len(self.halos)))

        self.update()



    def update(self):
        self.nrHalos = len(self.halos)

        for halo in xrange(self.nrHalos):
            self.positions[self.halos[halo].indices, -1] = halo

        self.calculateNrParticlesInHalos()
        self.percentageInHalos = float(self.nrParticlesInHalos)/self.nrParticles




    def calculateLinkingLength(self):
        """
        The linking length is chosen such that a fraction f of all particles
        is linked together with atleast one other particle

        For large groups > NrLinking we only calculate this for NrLinking particles
        By defult NrLinking = 10000.
        """
        tmpNrParticles = self.nrParticles
        delta = 1

        if self.nrParticles > self.nrLinking:
            tmpNrParticles = self.nrLinking
            delta = self.nrParticles/self.nrLinking

        LinkingLengths = []
        # Use all other particles or only the nrlinking subset
        next_particles = range(0, self.positions.shape[0])
        for i in xrange(delta, tmpNrParticles, delta):
            prevTmpLinkingLength = self.distance(0, i)
            for next_particle in next_particles:
                if next_particle != i:
                    tmpLinkingLength = self.distance(next_particle, i)
                    if tmpLinkingLength < prevTmpLinkingLength:
                        prevTmpLinkingLength = tmpLinkingLength
                        LinkingLengths.append(prevTmpLinkingLength)

        LinkingLengths.sort()
        return LinkingLengths[int(len(LinkingLengths)*self.f)]


    def calculateNrParticlesInHalos(self):
        s = 0
        for halo in self.halos:
            s += halo.size
        self.nrParticlesInHalos = s


    def save(self, name, delimiter=" "):
        f = open(name, "w")
        for halo in self.halos:
            for particle in halo.indices:
                line = ""
                for position in self.positions[particle]:
                    line += delimiter + str(position)
                    line += "\n"
                    f.write(line)
                    f.close()



    def printInformation(self):
        print "--------------------------------------------------------"
        print self.name
        print "--------------------------------------------------------"
        print "Nr of Particles:			 ", self.nrParticles
        print "Nr of Halos:			 ", self.nrHalos
        print "Nr of Particles in Halos:	 ", self.nrParticlesInHalos
        print "Percentage of Particles in Halos: ", self.percentageInHalos
        print "--------------------------------------------------------"



    def plotHalos(self, name):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        tmpHalo = CHalo()
        for halo in self.halos:
            tmpHalo + halo



        ax.scatter(self.positions[tmpHalo.indices, 0],
                   self.positions[tmpHalo.indices, 1],
                   self.positions[tmpHalo.indices, 2],
                   s=8,
                   c=self.positions[tmpHalo.indices, -1],
                   linewidths=0.1)

        ax.set_xlabel('X-position [um]')
        ax.set_ylabel('Y-position [um]')
        ax.set_zlabel('Z-position [um]')
        ax.set_title("Halos")
        ax.set_xlim(self.positions[:, 0].min(), self.positions[:, 0].max())
        ax.set_ylim(self.positions[:, 1].min(), self.positions[:, 1].max())
        ax.set_zlim(self.positions[:, 2].min(), self.positions[:, 2].max())

        plt.savefig(name)
        # plt.show()



    def plotParticles(self, name):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        color = self.positions[:, -1].copy()
        color[self.positions[:, -1] >= 0] = 0

        ax.scatter(self.positions[:, 0],
                   self.positions[:, 1],
                   self.positions[:, 2],
                   s=8,
                   c=color,
                   linewidths=0.1)

        ax.set_xlabel('X-position [um]')
        ax.set_ylabel('Y-position [um]')
        ax.set_zlabel('Z-position [um]')
        ax.set_title("Particles")

        ax.set_xlim(self.positions[:, 0].min(), self.positions[:, 0].max())
        ax.set_ylim(self.positions[:, 1].min(), self.positions[:, 1].max())
        ax.set_zlim(self.positions[:, 2].min(), self.positions[:, 2].max())

        plt.savefig(name)
        # plt.show()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="FOF analysis on a single file")
    parser.add_argument("filename")

    args = parser.parse_args()

    b = 0.28
    minNrHalos = 10
    linkingLength = 2

    halos = CHalos(args.filename, b, minNrHalos)
    start_time = time.time()
    halos.linkingLength = linkingLength
    halos.FOF()
    print("--- %s seconds ---" % (time.time() - start_time))


    savename = args.filename
    savename = savename.replace(" ", "-")
    savename = savename.replace("/", "-")

    #halos.save("outHalos.dat")
    halos.printInformation()
    halos.plotParticles(savename + "-particles.png")
    halos.plotHalos(savename + "-halos.png")
