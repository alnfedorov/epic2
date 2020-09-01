import logging
from math import pi, log, exp, ceil


class Background_island_probscore_statistics:
    #  External genomeLength and gapSize are in units of bps
    #  Internal genomeLength and gapSize are in units of windows
    #  only look at enrichment!
    def __init__(self, total_tags, windowSize, gapSize, window_pvalue,
                 genomeLength, bin_size):
        self.tag_density = total_tags * 1.0 / genomeLength
        self.window_size = windowSize
        # In bps.
        assert (gapSize % windowSize == 0)
        # gap is in bps
        self.gap_size = gapSize / windowSize
        # maximum number of windows allowed in a gap
        self.genome_length = int(ceil(float(genomeLength) / windowSize))
        self.average = self.tag_density * windowSize
        self.bin_size = bin_size

        # Precalculate the poisson, cumulative poisson values up to max (500, 2*self.average).
        self.max_index = max(500, int(2 * self.average))

        #print self.average, self.max_index;
        #self.fact=[];
        self.poisson_value = []
        self.window_score = []
        self.window_scaled_score = []
        self.precalculate_poisson(self.max_index)

        # gap_contribution needs min_tags_in_window
        # So the position of this line is critical.
        self.min_tags_in_window = 0
        sf = 1

        while (sf > window_pvalue):
            #print self.min_tags_in_window, sf;
            # print("poisson", self.min_tags_in_window, self.poisson_value[self.min_tags_in_window])
            sf -= self.poisson_value[self.min_tags_in_window]
            self.min_tags_in_window += 1
        #An alternative approach that uses the scipy package,
        #poisson.sf (n, lambda) = \sum_{i= n+1}^{\infty} p(i, lambda)
        #self.min_tags_in_window = int(self.average);
        #while (scipy.stats.poisson.sf(self.min_tags_in_window-1) > window_pvalue):
        #	self.min_tags_in_window += 1;

        #print "Window read count threshold: ", self.min_tags_in_window;

        self.gap_contribution = self.gap_factor()
        self.boundary_contribution = self.boundary()

        # print("self.gap_contribution", self.gap_contribution)
        # print("self.boundary_contribution", self.boundary_contribution)

        self.cumulative = []
        # new method, first fill the lowest score.
        prob = self.boundary_contribution * self.poisson_value[
            self.min_tags_in_window]
        score = -log(self.poisson_value[self.min_tags_in_window])
        #scaled_score = int(score/self.bin_size);
        scaled_score = int(round(score / self.bin_size))
        self.island_expectation = [0] * (scaled_score + 1)

        # print("island_expectation length", (scaled_score+1))
        self.island_expectation[scaled_score] = prob * self.genome_length
        # print("self.island_expectation[scaled_score]", self.island_expectation[scaled_score])
        # print("scaled_score", scaled_score)

        self.island_expectation[
            0] = self.boundary_contribution * self.genome_length / self.gap_contribution

        # self.root = self.find_asymptotics_exponent()
        #print "Exponent for Asymptotics: ", self.root;

    def precalculate_poisson(self, index):
        assert len(self.window_score) == len(self.window_scaled_score) == len(self.poisson_value)
        if len(self.poisson_value) > index:
            return

        for index in range(len(self.poisson_value), index + 1):
            prob = self.poisson(index, self.average)
            self.poisson_value.append(prob)
            if index < self.average:  # only want to look at enrichment
                self.window_score.append(0)
                self.window_scaled_score.append(0)
            else:
                if prob > 0:
                    self.window_score.append(-log(prob))
                    scaled_score = int(round(-log(prob) / self.bin_size))
                    self.window_scaled_score.append(scaled_score)
                else:  #prob is too small and deemed 0 by the system
                    self.window_score.append(1000)
                    scaled_score = int(round(1000 / self.bin_size))
                    self.window_scaled_score.append(scaled_score)

        self.max_index = len(self.poisson_value)

    def factorial(self, m):
        value = 1.0
        if m != 0:
            while m != 1:
                value = value * m
                m = m - 1
        return value

    # Return the log of a factorial, using Srinivasa Ramanujan's approximation
    def factln(self, m):
        if m < 20:
            value = 1.0
            if m != 0:
                while m != 1:
                    value = value * m
                    m = m - 1
            return log(value)
        else:
            return m * log(m) - m + log(m * (1 + 4 * m *
                                             (1 + 2 * m))) / 6.0 + log(pi) / 2

    def poisson(self, i, average):
        if i < 20:
            return exp(-average) * average ** i / self.factorial(i)
        else:
            exponent = -average + i * log(average) - self.factln(i)
            return exp(exponent)

    """
		gap is in the unit of windows. In each window in the gap, the
		window could have 0, 1, min_tags_in_windows-1 tags.
		say gap = 1, min_tags_in_window= 2, gap_factor = 1 +
		poission(0,a) + poisson(1, a), where 1 represents no gap,
		poisson(0,a) represents a window with 0 tag,
		poisson(1,a) represents a window with 1 tag,
		The gap contribution from each window is not independent
	"""

    def single_gap_factor(self):
        my_gap_factor = 0
        for i in range(self.min_tags_in_window):
            my_gap_factor += self.poisson_value[i]
        return my_gap_factor

    # gap contribution is bigger than 1
    def gap_factor(self):
        if self.gap_size == 0:
            return 1
        else:
            i = 1
            gap_contribution = 1
            # contribution from no gap
            my_gap_factor = self.single_gap_factor()
            # print("self.gap_size " * 100, self.gap_size)
            for i in range(1, int(self.gap_size + 1)):
                gap_contribution += pow(my_gap_factor, i)
            return gap_contribution

    def boundary(self):
        """
		The condition for boundary is a continuous region of
		unqualified windows longer than gap
		"""
        temp = self.single_gap_factor()
        temp = pow(temp, self.gap_size + 1)
        return temp * temp
        # start & end

    #forward method that memorize the calculated results.
    def background_island_expectation(self, scaled_score):
        current_max_scaled_score = len(self.island_expectation) - 1
        if scaled_score > current_max_scaled_score:
            #index is the scaled_score
            for index in range(current_max_scaled_score + 1, scaled_score + 1):
                temp = 0.0
                # i is the number of tags in the added window
                i = self.min_tags_in_window
                while int(round(index - self.window_score[i] / self.bin_size)) >= 0:
                    if self.poisson_value[i] == 0:
                        break
                    temp += self.poisson_value[i] * \
                            self.island_expectation[int(round(index - self.window_score[i] / self.bin_size))]
                    i += 1
                    if i >= self.max_index:
                        self.precalculate_poisson(i + 100)

                temp *= self.gap_contribution
                self.island_expectation.append(temp)
        return self.island_expectation[scaled_score]

    def generate_cumulative_dist(self, outfile=""):
        """
		Generate cumulative distribution: a list of tuples (bins, hist).
		"""
        self.cumulative = [0] * len(self.island_expectation)
        partial_sum = 0.0
        for index in range(1, len(self.island_expectation) + 1):
            complimentary = len(self.island_expectation) - index
            partial_sum += self.island_expectation[complimentary]
            # The end is outside of the index
            self.cumulative[complimentary] = partial_sum

        if outfile != "":
            fixpoint = int(len(self.island_expectation) / 2)
            outf = open(outfile, "w")
            outline = "# Score" + "\t" + "Expect # islands" + "\t" + "Cumulative # Islands" + "\t" + "Asymptotics" + "\n"
            outf.write(outline)
            for index in range(len(self.island_expectation)):
                outline = str(index * self.bin_size) + "\t" + str(
                    self.island_expectation[index]) + "\t" + str(
                    self.cumulative[index]) + "\n"
                outf.write(outline)
            outf.close()

    def find_island_threshold(self, e_value_threshold):
        """
		average is the average number of tags in a window:
		opt.tag_density * opt.window_size
		This one allows single-window islands.
		Returns the island threshold
		"""
        threshold = .0000001 * e_value_threshold
        current_scaled_score = len(self.island_expectation) - 1
        current_expectation = self.island_expectation[-1]
        assert (current_expectation ==
                self.island_expectation[current_scaled_score])
        interval = int(1 / self.bin_size)
        if len(self.island_expectation) > interval:
            partial_cumu = sum(self.island_expectation[-interval:-1])
        else:
            partial_cumu = sum(self.island_expectation)
        while (partial_cumu > threshold or partial_cumu < 1e-100):
            current_scaled_score += interval
            current_expectation = self.background_island_expectation(
                current_scaled_score)
            if len(self.island_expectation) > interval:
                partial_cumu = sum(self.island_expectation[-interval:-1])
            else:
                partial_cumu = sum(self.island_expectation)
            #for index in  xrange(len(self.island_expectation)):
            #print  index*self.bin_size, self.island_expectation[index];

        self.generate_cumulative_dist()
        for index in range(len(self.cumulative)):
            if self.cumulative[index] <= e_value_threshold:
                score_threshold = index * self.bin_size
                break
        return score_threshold


# Species:  hg38
# Window_size:  200
# Gap size:  600
# E value is: 1000.0
# Total read count: 20227554.0
# Genome Length:  3088286401
# Effective genome Length:  2625043440
# Window average: 1.54112146807
# opt.evalue: 1000.0
# Total read count: 20227554.0
# window_size: 200
# opt.gap: 600
# Window pvalue: 0.2
# genome_length: 2625043440
# ('self.average=', 1.541121468069877)
# ('self.poisson_value[0]=', 0.2141408146291121)
# ('window_pvalue', 0.2)
# ('poisson', 0, 0.2141408146291121)
# ('poisson', 1, 0.33001700661489664)
# ('poisson', 2, 0.25429814686118796)
# ('poisson', 3, 0.13063477780605437)
# ('self.gap_contribution', 3.5943004243263443)
# ('self.boundary_contribution', 0.5552199656493835)
# ('island_expectation length', 2990)
# ('self.island_expectation[scaled_score]', 366781.3884759084)
# ('scaled_score', 2989)
# Minimum num of tags in a qualified window:  4
# Generate the enriched probscore summary graph and filter the summary graph to get rid of ineligible windows
# Determine the score threshold from random background
# opt.evalue 1000.0
# The score threshold is:  19.112
# Make and write islands
# Total number of islands:  74071

# Find candidate islands exhibiting clustering ...
# /mnt/work/endrebak/software/anaconda/envs/py27/bin/python /home/endrebak/code/epic_paper/SICER/src/find_islands_in_pr.py -s hg38 -b /mnt/scratch/projects/epic_bencmarks/data/sicer_results/0h_exp1/no_bigwig/0h_exp1_chip-W200.graph -w 200 -g 600 -t 0.85 -e 1000 -f /mnt/scra
# tch/projects/epic_bencmarks/data/sicer_results/0h_exp1/no_bigwig/0h_exp1_chip-W200-G600.scoreisland
# Species:  hg38
# Window_size:  200
# Gap size:  600
# E value is: 1000.0
# Total read count: 2792406.0
# Genome Length:  3088286401
# Effective genome Length:  2625043440
# Window average: 0.212751222128
# opt.evalue: 1000.0
# Total read count: 2792406.0
# window_size: 200
# opt.gap: 600
# Window pvalue: 0.2
# genome_length: 2625043440
# ('self.average=', 0.21275122212834696)
# ('self.poisson_value[0]=', 0.8083572135909076)
# ('window_pvalue', 0.2)
# ('poisson', 0, 0.8083572135909076)
# ('self.gap_contribution', 2.990012655388544)
# ('self.boundary_contribution', 0.18231673960386235)
# ('island_expectation length', 1761)
# ('self.island_expectation[scaled_score]', 411536.58828177786)
# ('scaled_score', 1760)
# Minimum num of tags in a qualified window:  1
# Generate the enriched probscore summary graph and filter the summary graph to get rid of ineligible windows
# ('read_count', 1.0, 'start', 1949600)
# ('read_count', 1.0, 'start', 1949800)
# ('read_count', 1.0, 'start', 1950000)
# ('read_count', 3.0, 'start', 1950200)
# ('read_count', 1.0, 'start', 1950600)
# ('read_count', 1.0, 'start', 1950800)
# ('read_count', 1.0, 'start', 1951000)
# ('read_count', 1.0, 'start', 1951600)
# ('read_count', 1.0, 'start', 1951800)
# ('read_count', 1.0, 'start', 1952600)
# ('read_count', 1.0, 'start', 1953400)
# ('read_count', 1.0, 'start', 1953600)
# ('read_count', 1.0, 'start', 1954000)
# Determine the score threshold from random background
# opt.evalue 1000.0
# The score threshold is:  27.363
# Make and write islands
# Total number of islands:  16431


def compute_score_threshold(chip_counts, window_size, effective_genome_length,
                            gap_size, e_value):
    tag_density = chip_counts / effective_genome_length
    background = Background_island_probscore_statistics(
        chip_counts, window_size, gap_size, 0.2, effective_genome_length,
        0.001)
    score_threshold = background.find_island_threshold(e_value)

    logging.info("Score threshold: {}\n".format(score_threshold))

    min_tags_in_window = background.min_tags_in_window

    average_window_readcount = chip_counts * (
            window_size / float(effective_genome_length))

    return score_threshold, min_tags_in_window, average_window_readcount
