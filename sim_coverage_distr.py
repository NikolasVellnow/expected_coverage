
from dataclasses import dataclass
from timeit import default_timer as timer

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
        self.size = int(size)
        # numpy-array to save whether GRC maps at site (True) or not (False)
        self.site_grc_mappings = np.full(self.size, None)

    def __str__(self) -> str:
        complete_site_list = list(
            zip(range(0, self.size, 1), self.site_grc_mappings))
        if self.size <= 20:
            site_list = map(str, complete_site_list)
        else:
            center = int((self.size/2))
            first_part = complete_site_list[0:5]
            middle_part = complete_site_list[self.size -
                                             center: (self.size - center + 5)]
            last_part = complete_site_list[self.size-5: self.size]
            first_part.extend(middle_part)
            first_part.extend(last_part)
            site_list = map(str, first_part)
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
        self.site_grc_mappings = np.full(self.size, False)


class GenomeM2(Genome):
    """Class to represent genome according to model 2"""

    def __init__(self, size=100000, rel_size_p=0.3):
        super().__init__(size)
        self.rel_size_p = rel_size_p
        self.site_grc_mappings = np.full(self.size, False)
        self.size_grc_region = int(self.rel_size_p * self.size)
        # GRC sequences are at beginning of the chromosomes...
        self.grc_sites = np.arange(0, self.size_grc_region)
        self.site_grc_mappings[self.grc_sites] = True

    def get_rel_size_p(self) -> float:
        return self.rel_size_p

    def get_size_grc_region(self) -> int:
        return self.size_grc_region


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
            center = int((self.genome.get_size()/2))
            first_part = complete_site_list[0:5]
            middle_part = complete_site_list[self.genome.get_size() -
                                             center: (self.genome.get_size() - center + 5)]
            last_part = complete_site_list[self.genome.get_size(
            )-5: self.genome.get_size()]
            first_part.extend(middle_part)
            first_part.extend(last_part)
            site_list = map(str, first_part)
        return "\n".join(site_list)

    def map_reads(self, num_reads, read_len=150):
        """ Simulates a read mapping process generally """
        num_reads = int(num_reads)
        # tmp_read = Read(length=read_len) # Do i need the Read class?
        for _ in range(0, num_reads):
            rnd_pos = rnd.randint(
                0, self.genome.get_size() - (read_len-1))
            read_sites = np.arange(rnd_pos, rnd_pos+read_len)
            self.site_depths[read_sites] += 1

    def plot_coverage_hist(self, bins, read_length=150):
        """ Plots a histogram of simulated read coverage """
        # We exclude the beginning and end of genome to exclude edge effects
        read_length = int(read_length)
        plt.hist(
            self.site_depths[read_length: (self.genome.get_size()-read_length):1], bins=bins, density=False)
        plt.xlabel("Simulated sequencing coverage")
        plt.ylabel("Counts")
        plt.show()

    def cov_coeff_var(self, read_length=150):
        "Calculates the coefficient of variation of the coverage values"
        # We exclude the beginning and end of genome to exclude edge effects
        coeff_var = variation(
            self.site_depths[read_length: (self.genome.get_size()-read_length):1], axis=0, ddof=1)
        return coeff_var


class GenomeMappingM1(GenomeMapping):
    """ Class to represent read mapping process according to model 1 (null model)"""

    def plot_coverage_hist(self, bins, line_mean: float, read_length=150):
        """ Plots a histogram of simulated read coverage """
        # We exclude the beginning and end of genome to exclude edge effects
        read_length = int(read_length)
        plt.hist(
            self.site_depths[read_length: (self.genome.get_size()-read_length):1], bins=bins, density=False)
        plt.axvline(line_mean, color='k', linestyle='dashed', linewidth=2)
        plt.xlabel("Simulated sequencing coverage")
        plt.ylabel("Counts")
        plt.show()


class GenomeMappingM2(GenomeMapping):
    """ Class to represent reference genome according to model 2"""

    def __init__(self, genome: GenomeM2):
        super().__init__(genome)
        self.genome = genome

    def map_reads(self, num_reads, read_len=150):
        """ Simulates a read mapping process """
        num_reads = int(num_reads)
        tmp_read = Read(length=read_len)
        for _ in range(0, num_reads):
            rnd_float = rnd.random()  # generate random float between 0.0 and 1.0
            # Read has probability of p/(1+p) to come from GRC
            p_from_grc = self.genome.get_rel_size_p()/(1+self.genome.get_rel_size_p())
            if rnd_float < (p_from_grc):
                # tmp_read = Read(length=read_len, maps_to_grc=True) # Do I need the Read class?
                rnd_pos = rnd.randint(
                    0, self.genome.get_size_grc_region() - (read_len-1))
            else:
                # tmp_read = Read(length=read_len)
                rnd_pos = rnd.randint(
                    0, self.genome.get_size() - (read_len-1))
            read_sites = np.arange(rnd_pos, rnd_pos+tmp_read.length)
            self.site_depths[read_sites] += 1

    # TODO Write method to save the data

    def plot_coverage_hist(self, bins, line_mean: tuple[float, float], read_length=150):
        """ Plots a histogram of simulated read coverage """
        # We exclude the beginning and end of genome to exclude edge effects
        read_length = int(read_length)
        plt.hist(
            self.site_depths[read_length: (self.genome.get_size()-read_length):1], bins=bins, density=False)
        plt.axvline(line_mean[0], color='k', linestyle='dashed', linewidth=2)
        plt.axvline(line_mean[1], color='k', linestyle='dashed', linewidth=2)
        plt.xlabel("Simulated sequencing coverage")
        plt.ylabel("Counts")
        plt.show()

    # def cov_coeff_var(self, read_length=150):
    #     "Calculates the coefficient of variation of the coverage values"
    #     coeff_var = variation(
    #         self.site_depths[read_length: (self.genome.get_size()-read_length):1], axis=0, ddof=1)
    #     return coeff_var


def main():
    "main function"
    # variables for testing
    L = 150
    n_reads = 4.0 * 10**5
    n_bins = np.arange(20, 140, 1)
    G_R = 1.0 * 10**6
    rel_size_p = 0.3

    t0 = timer()
    g1 = GenomeM1(size=G_R)
    t1 = timer()
    print(f'Initialize GenomeM1 in {t1-t0} s')

    gm1 = GenomeMappingM1(genome=g1)

    t2 = timer()
    gm1.map_reads(num_reads=n_reads, read_len=L)
    t3 = timer()
    print(f'GenomeMappingM1 in: {t3-t2} s')

    print("GenomeM1:\n", g1)
    print("GenomeMapping1:\n", gm1)

    print(
        f'Coefficient of variation of Coverage of GenomeM1: {gm1.cov_coeff_var(read_length=L)}')

    t4 = timer()
    g2 = GenomeM2(size=G_R, rel_size_p=rel_size_p)
    t5 = timer()
    print(f'Initialize GenomeM2 in {t5-t4} s')

    gm2 = GenomeMappingM2(genome=g2)

    t6 = timer()
    gm2.map_reads(num_reads=n_reads, read_len=L)
    t7 = timer()
    print(f'GenomeMappingM2 in: {t7-t6} s')

    print("GenomeM2:\n", g2)
    print("GenomeMapping2:\n", gm2)

    print(
        f'Coefficient of variation of Coverage of GenomeM2: {gm2.cov_coeff_var(read_length=L)}')

    # gm1.plot_coverage_hist(bins=n_bins, line_mean=60, read_length=L)
    gm2.plot_coverage_hist(
        bins=n_bins, line_mean=(46.15, 92.31), read_length=L)

    # TODO create useful test cases, e.g. several instances that are independtenly manipulated etc.


if __name__ == "__main__":
    main()
