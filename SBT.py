from BloomFilter import BloomFilter


class Node(object):
    def __init__(self, bloom_filter):
        self.bloom_filter = bloom_filter
        self.left_child = None
        self.right_child = None

    def set_left_child(self, left_child):
        self.left_child = left_child

    def set_right_child(self, left_child):
        self.left_child = left_child

    def get_height(self):
        if self.left_child is None and self.right_child is None:
            return 0
        else:
            return max(self.left_child.get_height(), self.right_child.get_height())


class SBT(object):
    def __init__(self, bloom_filter_length, hash_functions, threshold, similarity_function):
        self.bloom_filter_length = bloom_filter_length
        self.hash_functions = hash_functions
        self.threshold = threshold
        self.similarity_function = similarity_function
        self.root = Node(BloomFilter(self.bloom_filter_length, self.hash_functions, similarity_function))

    def insert(self, bloom_filter):

        pass

    def query(self, bloom_filter):
        pass


