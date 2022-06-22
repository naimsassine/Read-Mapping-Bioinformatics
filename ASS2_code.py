from operator import itemgetter
import numpy
import time
import matplotlib.pyplot as plt
import pdb

class Genome_Hashing(object):
  # This class is built to hash the reference genome, 
  # and enable a faster access
  def __init__(self, gen, len_kmer, interval):
    # method that initialises the hash table and constructs it by going through the genome 'gen'
    self.gen = gen
    self.len_kmer = len_kmer
    self.interval = interval
    self.index = {}
    for i in range(0, len(gen) - len_kmer + 1, interval):
        kmr = gen[i: i + len_kmer]
        if kmr in self.index:
            self.index[kmr].append(i) 
        else:
            self.index[kmr] = [i]

  def query(self, p):
      # function that will query the demander k-mer and output the different positions
    return self.index.get(p[:self.len_kmer], [])


class SaE_read_mapping(object):

  def __init__(self, genome_path, reads_path, kmer_size, seed_int):
    # we initialise the different attributes that we will be using later
    self.gen_filepath = genome_path
    self.reads_filepath = reads_path

    self.seed_length = kmer_size
    self.seed_interval = seed_int

    self.gen_name = ""
    self.gen_sequence = ""

    self.read = dict()

  def parse_genome(self):
    # parse the genome file
    self.gen_name = "dm6"
    file = open(self.gen_filepath)
    file.readline()
    self.gen_sequence = file.read().replace("\n", "").upper()
    file.close()

  def parse_reads(self): 
    # parse the reads, create the dicitonary with the wanted info
    with open(self.reads_filepath) as f:
      line = f.readline().rstrip()
      while line:
        if line.startswith("@") or line.startswith("+"):
          prev = line[0]
          if line.startswith("@"):
            name = line[1:].split("/")[0]
        else:
          if prev.startswith("@"):
            strand = line
          elif prev.startswith("+"):
            self.read[name] = dict()
            self.read[name]["strand"] = strand
            self.read[name]["seeds"] = list()
            self.read[name]["positions_in_ref"] = list()

        line = f.readline().rstrip()

  def detect_seed_in_genome(self, read_string, index):
    # match seeds to k-mers in the genome
    seed_info = list()
    for i in range (0, len(read_string), self.seed_interval):
      end = min(len(read_string), i + self.seed_length)
      cur_seed = read_string[i: end]
      positions_list = index.query(cur_seed)
      if len(positions_list) > 0:
        seed_info.append((positions_list, i, end))
    
    return seed_info

  def fill_info_seeds(self, index_forward):
    # once the seeds matched, we fill in the infromation for each seed of each read, that we will later need
    for key, value in self.read.items():
      list_of_seeds = list()
      list_of_seeds.append(self.detect_seed_in_genome(value["strand"], index_forward))

      max_len, max_list, max_ind = max((len(x), x, list_of_seeds.index(x)) for x in list_of_seeds)
      
      for seed_entry in max_list:
        for position in seed_entry[0]:
          seed_info = dict()
          seed_info["pos_ref"] = int(position)
          seed_info["start_read"] = seed_entry[1]
          seed_info["end_read"] = seed_entry[2]
            
          self.read[key]["seeds"].append(seed_info)
  
  def seed(self):
    # Seed function
    index_forward = Genome_Hashing(self.gen_sequence, self.seed_length, self.seed_interval)

    self.fill_info_seeds(index_forward)

  def extend_seed(self, read_key, candidate, extension):
    # extend the seed and the k-mer of the genome to generate a score
    value = self.read[read_key]
    strand = value["strand"] 
    pos_in_ref = candidate["pos_ref"]
    read_start = candidate["start_read"]
    read_end = candidate["end_read"]

    added_left_pos_read = 0
    added_right_pos_read = 0

    pos_genome_left = 0
    pos_genome_right = 0

    added_left = "A"
    added_right = "A"

    if read_start - extension < 0:
      if pos_in_ref - extension < 0:
        added_left_pos_read = 0
        pos_genome_left = 0
      else:
        added_left_pos_read = 0
        pos_genome_left = pos_in_ref - read_start
    else:
      if pos_in_ref - extension < 0:
        added_left_pos_read = read_start - pos_in_ref
        pos_genome_left = 0
      else:
        added_left_pos_read = read_start - extension
        pos_genome_left = pos_in_ref - extension
   

    if read_end + extension > len(strand):
      if pos_in_ref + (read_end - read_start) + extension >= len(self.gen_sequence):
        added_right_pos_read = len(strand) - 1
        pos_genome_right = len(self.gen_sequence) -1
      else:
        added_right_pos_read = len(strand) - 1
        pos_genome_right = pos_in_ref + (len(strand) - read_start)
    else:
      if pos_in_ref + (read_end - read_start) + extension >= len(self.gen_sequence):
        added_right_pos_read = len(self.gen_sequence) - pos_in_ref
        pos_genome_right = len(self.gen_sequence) -1
      else:
        added_right_pos_read = read_start + extension
        pos_genome_right = pos_in_ref + extension



    score = 0
    if strand[added_left_pos_read] != self.gen_sequence[pos_genome_left]:
      score = score - 1

    if strand[added_right_pos_read] != self.gen_sequence[pos_genome_right]:
      score = score - 1 

    return score

  def find_best_candidates(self, read_key):
    # outputs the best seeds of each read
    # [{'end_read': 11, 'pos_ref': 1135686, 'start_read': 1}, ...]
    seeds = self.read[read_key]["seeds"]
    candidates = {}
    number_of_iterations = len(self.read[read_key]["strand"]) - self.seed_length
    j = 0
    while j < number_of_iterations:
      j = j+1
      for i in range(len(seeds)):
        if i in candidates:
          candidates[i] += self.extend_seed(read_key, seeds[i], j)
        else:
          candidates[i] = self.extend_seed(read_key, seeds[i], j)
    
    # we sort the candidates by the ones who have the least errors
    sorted_candidates = sorted(candidates.items(), key=itemgetter(1))[::-1]
    # we take the first 5 percent of candidates
    num = int(len(sorted_candidates)*0.05)
    chosen_candidates = sorted_candidates[:num]

    # then we remove the candidates that have too many errors
    # we set as margin : 3 errors
    final_candidates = []
    for cand in chosen_candidates:
      if cand[1] > -3:
        final_candidates.append(cand)

    return final_candidates

  def extend(self):
    # return a dictionary with : read : the seeds and their particularities
    result = {}
    for key, value in self.read.items():
      seeds = value["seeds"]
      n_bl = len(seeds)
      candidates = []
      if n_bl == 0:
        candidates = []
      else:
        candidates = self.find_best_candidates(key)
      for match in candidates:
        if key not in result : 
          result[key] = [[seeds[match[0]], match[1]]]
        else:
          result[key] = result[key] + [[seeds[match[0]], match[1]]]
    
    return result

# small dataset
reads_filepath = "/content/Bioinfo_ASS2/10k_reads.fastq"
genome_filepath = "/content/Bioinfo_ASS2/dm6_chr2L.fa"
reduced = SaE_read_mapping(genome_filepath, reads_filepath, 10, 1)
reduced.parse_genome()
reduced.parse_reads()
reduced.seed()
results_reduced = reduced.extend()


# % of reads mapped = numberofreadsmapped/numberofuniquereads
perc_reads_mapped = (len(results_reduced)/len(reduced.read)) * 100

# percentage of multireads
multi_reads = 0
for key in results_reduced:
  contains = 0
  for seed in results_reduced[key]:
    if seed[1] == 0:
      contains += 1
  if contains > 1 :
    multi_reads += 1
  

# number of full matches 
full_matches = 0
for key in results_reduced:
  contains = False
  for seed in results_reduced[key]:
    if seed[1] == 0:
      contains = True
  if contains == True:
    full_matches += 1


# entire dataset
reads_filepath_full = "/content/Bioinfo_ASS2/all_reads.fastq"
genome_filepath_full = "/content/Bioinfo_ASS2/dm6.fa"
full = SaE_read_mapping(genome_filepath_full, reads_filepath_full, 7, 30)
full.parse_genome()
full.parse_reads()
full.seed()
results_full = full.extend()