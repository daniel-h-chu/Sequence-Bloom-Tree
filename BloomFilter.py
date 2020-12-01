class BloomFilter(object):
    def __init__(self, bloom_filter_length, hash_functions):
        self.bloom_filter_length = bloom_filter_length
        self.hash_functions = hash_functions
        self.array = [0] * bloom_filter_length

    def insert(self, kmer):
        for hash_function in self.hash_functions:
            self.array[hash_function(kmer)] = 1

    def query(self, kmer):
        for hash_function in self.hash_functions:
            if self.array[hash_function(kmer)] == 0:
                return False
        return True
