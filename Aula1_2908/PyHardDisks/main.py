from collections import namedtuple
from random import random
from math import sqrt
#import matplotlib.pyplot as plt
import time

L, R = 50, 1
nplacements = 1000
tolerance = 100
disks = []
Disk = namedtuple('Disk', 'x y')

def main():
    start_time = time.time()
    run()
    print("--- %s seconds ---" % (time.time() - start_time))

def run():
    #file = open("PyNumerOfDisksPlaced.txt", "w")
    #for i in range(nplacements):
    #    file.write("{}\n".format(placeDisks()))
    #file.close()

    data = [placeDisks() for i in range(nplacements)]
    #data = np.loadtxt("placingcircles.dat")
    #hist(x, bins=None, range=None, density=None, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, color=None, label=None, stacked=False, normed=None, hold=None, data=None, **kwargs)
    plt.hist(data, facecolor='green', alpha=0.75)
    plt.grid(True)
    plt.savefig("/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS COMPUTACIONAIS C/metcompc/Aula1_2908/Results/PyNumberOfDisksPlaced.svg")

def placeDisks():
    disks = []
    attempts = 0
    while attempts < tolerance:
        disk = Disk(x = randomPosition(), y = randomPosition())
        if not overlapping(disk):
            disks.append(disk)
        else:
            attempts += 1
    return len(disks)


#improvement: create class disk and define equality
def overlapping(disk):
    for other in disks:
        xDiff, yDiff = abs(other.x - disk.x), abs(other.y - disk.y)
        if sqrt(pow(xDiff, 2) + pow(yDiff, 2)) <= 2*R:
            return True
    return False

def randomPosition():
    return (L-2)*random() + 1

main()