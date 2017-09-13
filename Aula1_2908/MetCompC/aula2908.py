from collections import namedtuple
from random import random,seed
import numpy as np
import matplotlib.pyplot as plt

L, R, nruns = 50, 1, 1000
tolerance = 1000
disks = []
Disk = namedtuple('Disk', 'x y')

#improvement: create class disk and define equality
def overlapping(disk):
    for other in disks:
        x_sup, y_sup = abs(other.x - disk.x), abs(other.y - disk.y)
        if x_sup < 2*R or y_sup < 2*R:
            return True
    return False

def randomPosition():
    return (L-2)*random() + 1
#seed()
def place_disks(i):
    attempts = 0
    while attempts < tolerance:
        disk = Disk(x = randomPosition(), y = randomPosition())
        if not overlapping(disk):
            disks.append(disk)
        else:
            attempts += 1
    return len(disks)

def run():
    #data = [place_disks(i) for i in range(nruns)]
    data = np.loadtxt("placingcircles.dat")
    #print(*data)
    plt.hist(data, 50, normed=1, facecolor='green', alpha=0.75)
    plt.grid(True)
    plt.savefig("placingcircles.svg")
    #plot.show()

def aula1Difusao():
    ts = [2,10,50,100,200,1000,2000]
    for i in range(5):
        k = (i+1)*0.1
        file = "aula1difusaoK{:.3f}.dat".format(k)
        print(file)
        x,y = np.loadtxt(file,unpack=True)
        for i in range(7):
            plt.plot(x[(1001*i):(1001*(i+1))],y[(1001*i):(1001*(i+1))], label="t={}".format(ts[i]))
        plt.legend()
        plt.title("K={:.1f}".format(k))
        plt.savefig(file.replace(".dat", ".svg"));
        plt.gcf().clear()
aula1Difusao()

#x = np.arange(0, 5, 0.1);
#y = np.sin(x)
#plt.plot(x, y)
