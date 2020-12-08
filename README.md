README

Daniel Chu - dchu12@jhu.edu

Introduction:

This module contains an base implementation of a Sequence Bloom Tree (SBT) along with several different 
variants of the SBT and variants of insertion and querying algorithms. 


How to Use:

Run main.py to simulate the entire benchmarking process (Reading sequencing data -> Inserting sequencing 
data into SBT -> Querying sequences from the SBT -> Saving the SBT). The amount of time spent in each step, 
the precision (TP / (TP + FP)) of the queries is also reported along with the uncompressed size of the SBT
after saving. Paramters in dictionary p can be adjusted to change  the benchmarking process or change the 
SBT implementation.

Run pipelined_main.py to automate the running of several benchmarking simulations (i.e. running main.py with 
different params)
    
    
Parameter Descriptions: 

    bloom_filter_length: int (Must be positive)
        Size of the bloom filters used in the SBT

    k: int (Must be positive)
        Size of k-mer

    bits_to_check: int (Must be positive) or None
        Number of bits to compare between filters when computing similarity if we are using clustering. If
        the parameter is set to None, then all the bits will be compared

    num_sequences: int (Must be positive)           
        Number of sequences to insert into the SBT

    threshold: float (Must be between 0 and 1 inclusive)              
        Fraction of queried kmers that must be present in a bloom filter in order to containue
        querying or to return the filter as a hit

    sequence_len: int (Must be positive)               
        How many bps of each sequence we want to insert into the SBT

    query_size: int (Must be positive)   
        How many bps of each sequence we want to query from the SBT

    num_queries: int (Must be positive)   
        How many queries we want to perform 

    node_class: str (Must be one of ["Base", "SSBT", "HowDet"])
        Type of SBT to use. "Base" generated a base SBT, "SSBT" generated a Split-SBT,
        and "HowDet" generated a HowDet-SBT.
        
    insert_method: str (Must be one of ["Greedy", "Cluster1", "Cluster2"])
        Insertion method to use. "Greedy" inserts nodes 1 by 1 by traversing the tree down the most similar
        child. "Cluster1" inserts all nodes at the same time by computing the pairwise similarity between the nodes
        and creating a parent node between the two most similar nodes and repeat until we have 1 node left. "Cluster2"
        runs similarly to "Cluster1" but all nodes are paired together before the parents are considered for pairing
        again.
        
    query_method: str (Must be one of ["Normal", "Fast"])
        Query method to use. "Normal" hashes the kmers at every filter we query and we check whether or not the
        index that the kmer hashes to tells us that the kmer is present. "Fast" hashes the kmers only once and 
        instead keeps track of a a list of indices that the kmers hash to.
        
    similarity_function: func (Must be one of [hamming, and_hamming, cosine, jaccard, manhattan, euclidian, dice])
        Similarity function to use when inserting nodes. Nodes being more similar result in similarity_function 
        returning a more positive. and_hamming is recommended for SSBT and HowDet. cosine is recommended for Base
        
    hash_functions: list(func) (Must be [hash])
        List of hash functions to use inside the bloom filters. Since it is difficult to construct a hash function
        that is independent and as fast as python's hash(), it is recommended to use only python's hash functino
        
    hash_fraction: float (Must be between 0 and 1 inclusive)
        Proportion of kmers that are hashed into the bloom filter. If hash_fraction is less than
        one, then some kmers are not inserted into the bloom filter. Otherwise, all kmers are inserted. This
        parameter can be used to simualte fractional hash functions (e.g. 1 hash function and a hash fraction of
        1/2 gives you 1/2 of a hash functino)
        
    print_sbt: bool
        If true, then we print the SBT after all the benchmarking metrics are reported

    print_type: str (Must be one of ["Bits", "Names"])
        If print_sbt is true, then we print either the bits of the filters themselves (print_type="Bits") or 
        we print the experiment name corresponding to each filter (print_type="Names")
        
    sequence_prefix: str
        The prefix of your genome files. For example 
        If your genome files are named "file/genome0", "file/genome1", ... then sequence_prefix = "file/genome"
        
    pandas_location: str
        Where your experiment outputs should be saved as an excel file

    sbt_location: str
        Where the SBT should be saved

    benchmark_name: str
        What to name this experiment
        
File Descriptions:

    SBT/SBT.py
        This file contains the SBT class implementation. Variants of SBT are implemented based on what kind of 
        Node the SBT uses (i.e. if the SBT uses BaseNode, then we get a Base SBT and if the SBT uses SSBTNode, 
        then we get a Split-SBT). The SBT calls the Node's insertion and querying methods which are implemented
        based on algorithms described in several papers. 

    SBT/BaseNode.py
        This file contains the BaseNode class implementation. The node developed based on the SBT described in
        Solomon & Kingsford (2015)

    SBT/SSBTNode.py
        This file contains the SSBTNode class implementation. The node developed based on the Split-SBT described in
        Solomon & Kingsford (2018)       

    SBT/HowDetNode.py
        This file contains the HowDetNode class implementation. The node developed based on the HowDet-SBT described in
        Harris & Medvedev (2019)     

    main.py
        This file contains calls to util.py that execute general process of benchmarking. We print the amount of time 
        it takes for each step of the benchmarking. The main file also contains a dictionary p that contains
        parameters that can be adjusted to change the benchmarking process or change the SBT implementation.        

    utils.py
        Implementation of functions that are important for benchmarking (like reading in the files themselves,
        converting sequences to stuff insertable into the SBT). The file also contains additional 
        optional hash functions and similarity functions that can be set as a parameter to the benchmarking or SBT.    

    generate_test_data.py
        Generate completely random strings of 'ACGT' of custom length    

    test.py
        Random non-rigorous end to end tests for SBT

