class BloomFilter(object):
    def __init__(self, bloom_filter_length, hash_functions, similarity_function):
        self.bloom_filter_length = bloom_filter_length
        self.hash_functions = hash_functions
        self.similarity_function = similarity_function
        self.array = [0] * bloom_filter_length

    def insert(self, kmer):
        for hash_function in self.hash_functions:
            self.array[hash_function(kmer)] = 1

    def query(self, kmer):
        for hash_function in self.hash_functions:
            if self.array[hash_function(kmer)] == 0:
                return False
        return True

    def union(self, bloom_filter):
        for i in range(self.bloom_filter_length):
            self.array[i] = self.array[i] | bloom_filter.array

    def similarity(self, bloom_filter):
        return self.similarity_function(self.array, bloom_filter.array)
