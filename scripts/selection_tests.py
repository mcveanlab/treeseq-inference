import msprime
import collections
import math
from scipy import stats
import matplotlib.pyplot as pyplot
from SFS import incrementalSFS

def H(n):
    """Returns an approximate value of n-th harmonic number.

       http://en.wikipedia.org/wiki/Harmonic_number
    """
    # Euler-Mascheroni constant
    gamma = 0.57721566490153286060651209008240243104215933593992
    return gamma + math.log(n) + 0.5/n - 1./(12*n**2) + 1./(120*n**4)

def Theta_W(sfs, n):
    """
    Calculate the (mutationless) Watterson estimator
    """
    return sum(sfs.values()) / H(n)

def Theta_pi(sfs, n):
    """
    Calculate the (mutationless) nucleotide diversity
    2 sum(i(n−i)*ξi) / n / (n-1)
    """
    return 2 * sum([i*(n-i)*sfs[i] for i in sfs.keys()]) / n / (n-1)

def Theta_H(sfs, n):
    """
    Calculate the theta estimate used in Fay & Wu's H
    2 * sum(i^2 ξi) / n / (n-1)
    """
    return 2 * sum([(i**2) * sfs[i] for i in sfs.keys()]) / n / (n-1)
    
#globals    
T_W  = []
T_pi = []
T_H  = []
D    = []
FayWuH=[]

def save_estimates(start, end, sfs, n):
    watt = Theta_W(sfs, n)
    ndiv = Theta_pi(sfs, n)
    tH = Theta_H(sfs, n)
    T_W.append((watt, end-start))
    T_pi.append((ndiv, end-start))
    T_H.append((tH, end-start))
    D.append(((end+start)/2, ndiv - watt))
    FayWuH.append(((end+start)/2, ndiv - tH))


recombination_rate=1e-6
mutation_rate=1e-6
ts = msprime.simulate(10000, 
        Ne=5000, length=10000, recombination_rate=recombination_rate)        
#or try generating an identical neutral ftprime simulation using 
# python3 ./examples/examples.py -T 10000 -N 5000 -r 1e-6 -L 10000 -a 0.0001 -b 0.0001 -k 5000 -t neutral_ts
ts = msprime.load("../ftprime/neutral_ts")

print("Length: ", ts.get_sequence_length(), ", n samples: ", ts.num_samples)
print(stats.describe([tree.get_total_branch_length() for tree in ts.trees()]))

incrementalSFS(ts, save_estimates)



print("Watterson: ", sum([w[0]*w[1] for w in T_W])/ts.get_sequence_length())
print("Pairwise: ", sum([p[0]*p[1] for p in T_pi])/ts.get_sequence_length())
print("H: ", sum([p[0]*p[1] for p in T_H])/ts.get_sequence_length())

pyplot.plot(*zip(*D))
pyplot.plot(*zip(*FayWuH))
pyplot.legend(["Tajima's D", "Fay & Wu's H"])
pyplot.xlim(0,ts.get_sequence_length())
pyplot.marker="o"
pyplot.show()