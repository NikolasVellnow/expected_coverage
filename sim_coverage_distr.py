
from dataclasses import dataclass
from timeit import default_timer as timer
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import numpy.random as rnd
from scipy.stats import variation


@dataclass
class Read:
    """ Class to represent sequencing reads """
    length: int = 150
    maps_to_grc: bool = False


class Genome:
    """ Class to represent reference genome generally"""

    def __init__(self, size=100000):
        self.size = size
        # numpy-array to save whether GRC maps at site (True) or not (False)
        self.site_grc_mappings = np.full(size, None)

    def __str__(self) -> str:
        site_list = list(zip(range(0, self.size, 1), self.site_grc_mappings))
        site_list = map(str, site_list)
        return "\n".join(site_list)

    def get_size(self) -> int:
        """Returns genome size"""
        return self.size

    def get_site_grc_mappings(self):
        """Returns np.array of booleans for wether the GRC maps to a site or not"""
        return self.site_grc_mappings


class GenomeM1(Genome):
    """Class to represent genome according to model 1"""

    def __init__(self, size=10000):
        super().__init__(size)
        self.site_grc_mappings = np.full(size, False)


class GenomeM2(Genome):
    """Class to represent genome according to model 2"""

    def __init__(self, size=100000, rel_size_p=0.3):
        super().__init__(size)
        self.rel_size_p = rel_size_p
        self.site_grc_mappings = np.full(self.size, False)
        self.G_G = int(self.rel_size_p * self.size)
        # GRC sequences are at beginning of the chromosomes...
        self.grc_sites = np.arange(0, self.G_G)
        self.site_grc_mappings[self.grc_sites] = True


class Genome_m2:
    """ Class to represent reference genome according to model 2"""

    def get_grc_mapping(self, position):
        """ Returns boolean whether grc maps to site or not """
        return self.site_grc_mappings[position]

    def map_reads(self, num_reads, read_len=150):
        """ Simulates a read mapping process """
        tmp_read = Read(length=read_len)
        for _ in range(0, num_reads):
            rnd_float = rnd.random()  # generate random float between 0.0 and 1.0
            # Read has probability of p/(1+p) to come from GRC
            if rnd_float < (self.rel_size_p/(1+self.rel_size_p)):
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

    def cov_coeff_var(self, read_length=150):
        "Calculates the coefficient of variation od the coverage values"
        coeff_var = variation(
            self.site_depths[read_length: (self.size-read_length):1], axis=0, ddof=1)
        return coeff_var


class GenomeMapping:
    "Class to represent read mapping process generally"

    def __init__(self, genome: Genome):
        self.genome = genome
        # numpy-array to save how often each site has been sequenced
        self.site_depths = np.zeros(self.genome.get_size(), dtype=np.int16)

    def __str__(self) -> str:
        complete_site_list = list(zip(range(0, self.genome.get_size(), 1),
                                      self.genome.get_site_grc_mappings(), self.site_depths))
        if self.genome.get_size() <= 20:
            site_list = map(str, complete_site_list)
        else:
            first_part = complete_site_list[0:10]
            second_part = complete_site_list[self.genome.get_size(
            )-10: self.genome.get_size()]
            first_part.extend(second_part)
            site_list = map(str, first_part)

        return "\n".join(site_list)

    def map_reads(self, num_reads, read_len=150):
        """ Simulates a read mapping process """
        tmp_read = Read(length=read_len)
        for _ in range(0, num_reads):
            rnd_pos = rnd.randint(
                0, self.genome.get_size() - (tmp_read.length-1))
            read_sites = np.arange(rnd_pos, rnd_pos+tmp_read.length)
            self.site_depths[read_sites] += 1

    def plot_coverage_hist(self, bins, read_length=150):
        """ Plots a histogram of simulated read coverage """
        # We exclude the beginning and end of genome to exclude edge effects
        plt.hist(
            self.site_depths[read_length: (self.genome.get_size()-read_length):1], bins=bins, density=False)
        plt.show()

    def cov_coeff_var(self, read_length=150):
        "Calculates the coefficient of variation of the coverage values"
        coeff_var = variation(
            self.site_depths[read_length: (self.genome.get_size()-read_length):1], axis=0, ddof=1)
        return coeff_var


class GenomeMappingM1(GenomeMapping):
    """ Class to represent read mapping process according to model 1 (null model)"""

    # def __init__(self, genome: GenomeM1):
    #     self.genome = genome
    #     # numpy-array to save how often each site has been sequenced
    #     self.site_depths = np.zeros(self.genome.get_size(), dtype=np.int16)

    # def __str__(self) -> str:
    #     site_list = list(zip(range(0, self.genome.get_size(), 1),
    #                      self.genome.get_site_grc_mappings(), self.site_depths))
    #     site_list = map(str, site_list)
    #     return "\n".join(site_list)

    def map_reads(self, num_reads, read_len=150):
        """ Simulates a read mapping process """
        tmp_read = Read(length=read_len)
        for _ in range(0, num_reads):
            rnd_pos = rnd.randint(
                0, self.genome.get_size() - (tmp_read.length-1))
            read_sites = np.arange(rnd_pos, rnd_pos+tmp_read.length)
            self.site_depths[read_sites] += 1

    def plot_coverage_hist(self, bins, read_length=150):
        """ Plots a histogram of simulated read coverage """
        # We exclude the beginning and end of genome to exclude edge effects
        plt.hist(
            self.site_depths[read_length: (self.genome.get_size()-read_length):1], bins=bins, density=False)
        plt.show()

    def cov_coeff_var(self, read_length=150):
        "Calculates the coefficient of variation of the coverage values"
        coeff_var = variation(
            self.site_depths[read_length: (self.genome.get_size()-read_length):1], axis=0, ddof=1)
        return coeff_var


class GenomeMappingM2:
    """ Class to represent reference genome according to model 2"""

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
        for _ in range(0, num_reads):
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

    def cov_coeff_var(self, read_length=150):
        "Calculates the coefficient of variation od the coverage values"
        coeff_var = variation(
            self.site_depths[read_length: (self.size-read_length):1], axis=0, ddof=1)
        return coeff_var


def main():
    "main function"
    # variables for testing
    L = 150
    n_reads = 400000
    # n_bins = 5  # np.arange(20, 140, 1)
    G_R = 1000000
    rel_size_p = 0.3

    read_1 = Read()
    print(read_1)

    t0 = timer()
    g1 = GenomeM1(size=G_R)
    t1 = timer()
    print(f'Initialize genome in {t1-t0} s')

    gm1 = GenomeMappingM1(genome=g1)
    t2 = timer()
    gm1.map_reads(num_reads=n_reads, read_len=L)
    t3 = timer()
    print(f'Map reads in {t3-t2} s')

    print("GenomeM1:\n", g1)
    print("GenomeMapping:\n", gm1)

   # print(f'Coefficient of variation: {g2.cov_coeff_var(read_length=L)}')

    # g2.plot_coverage_hist(bins=n_bins, read_length=L)

    # TODO create useful test cases, e.g. several instances that are independtenly manipulated etc.


if __name__ == "__main__":
    main()
