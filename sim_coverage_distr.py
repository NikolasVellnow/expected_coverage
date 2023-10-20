
import itertools as it
from timeit import default_timer as timer

import matplotlib.pyplot as plt
import numpy as np
import numpy.random as rnd


class Read:
    """ Class to represent sequencing reads """

    def __init__(self, length=150, maps_to_grc=False):
        self.length = length
        self.maps_to_grc = maps_to_grc

    def __str__(self):
        return f"{self.length},{self.maps_to_grc}"


class Genome_m1:
    """ Class to represent reference genome according to model 1 (null model)"""

    def __init__(self, size=int(1.5*10**9)):
        self.size = size
        self.site_depths = np.zeros(size, dtype=np.int16)
        self.site_grc_mappings = np.full(size, False)

    def __str__(self) -> str:
        site_list = list(zip(self.site_grc_mappings, self.site_depths))
        site_list = map(str, site_list)
        return "\n".join(site_list)

    def get_depth(self, position):
        """ Returns site object at specified position"""
        return self.site_depths[position]

    def set_depth(self, position, new_site):
        """ Updates site object at specified position"""
        self.site_depths[position] = new_site

    def map_reads(self, num_reads, read_len=150):
        """ Simulates a read mapping process """
        tmp_read = Read(length=read_len)
        for read in range(0, num_reads):
            rnd_pos = rnd.randint(0, self.size - (tmp_read.length-1))
            read_sites = np.arange(rnd_pos, rnd_pos+tmp_read.length)
            self.site_depths[read_sites] += 1

    def plot_coverage_hist(self, bins, read_length=150):
        """ Plots a histogram of simulated read coverage """
        # We exclude the beginning and end of genome to exclude edge effects
        plt.hist(
            self.site_depths[read_length: (self.size-read_length):1], bins=bins, density=False)
        plt.show()


class Genome_m2:
    """ Class to represent reference genome according to model 2"""

    def __init__(self, size=int(1.5*10**9), p=0.3):
        self.size = size
        self.p = p
        self.site_depths = np.zeros(size, dtype=np.int16)
        self.site_grc_mappings = np.full(size, False)
        self.G_G = int(p * size)
        # GRC sequences are at beginning of the chromosomes...
        self.grc_sites = np.arange(0, self.G_G)

        self.site_grc_mappings[self.grc_sites] = True

    def __str__(self) -> str:
        site_list = list(zip(self.site_grc_mappings, self.site_depths))
        site_list = map(str, site_list)
        return "\n".join(site_list)

    def get_depth(self, position):
        """ Returns site object at specified position"""
        return self.site_depths[position]

    def set_depth(self, position, new_site):
        """ Updates site object at specified position"""
        self.site_depths[position] = new_site

    def get_grc_mapping(self, position):
        """ Returns boolean whether grc maps to site or not """
        return self.site_grc_mappings[position]

    def map_reads(self, num_reads, read_len=150):
        """ Simulates a read mapping process """
        tmp_read = Read(length=read_len)
        for read in range(0, num_reads):
            rnd_float = rnd.random()  # generate random float between 0.0 and 1.0
            # Read has probability of p/(1+p) to come from GRC
            if rnd_float < (self.p/(1+self.p)):
                tmp_read = Read(length=read_len, maps_to_grc=True)
                # TODO self.grc.sites ...
                rnd_pos = rnd.randint(0, self.G_G - (tmp_read.length-1))
            else:
                tmp_read = Read(length=read_len)
                rnd_pos = rnd.randint(0, self.size - (tmp_read.length-1))
            read_sites = np.arange(rnd_pos, rnd_pos+tmp_read.length)
            self.site_depths[read_sites] += 1

    # TODO Write method to save the data

    def plot_coverage_hist(self, bins, read_length=150):
        """ Plots a histogram of simulated read coverage """
        # We exclude the beginning and end of genome to exclude edge effects
        plt.hist(
            self.site_depths[read_length: (self.size-read_length):1], bins=bins, density=False)
        plt.show()


# variables for testing
L = 150
n_reads = 400000
n_bins = np.arange(20, 140, 1)
G_R = 1000000


t0 = timer()
g2 = Genome_m2(size=G_R, p=0.3)
t1 = timer()
print(f'Initialize genome in {t1-t0} s')


t2 = timer()
g2.map_reads(num_reads=n_reads, read_len=L)
t3 = timer()
print(f'Map reads in {t3-t2} s')

# print(g2)


g2.plot_coverage_hist(bins=n_bins, read_length=L)

"""
t0 = timer()
g1 = Genome_m1(10)
t1 = timer()
print(f'Initialize genome in {t1-t0} s')


t2 = timer()
g1.map_reads(num_reads=20, read_len=3)
t3 = timer()
print(f'Map reads in {t3-t2} s')



# 436797999
# 40000000


g1.plot_coverage_hist(range(0, 100))

"""
