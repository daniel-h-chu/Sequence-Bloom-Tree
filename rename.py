""" Rename simulated genomes once generated so that the main.py pipeline can use them as inputs """
import os
for i in range(500):
    if i % 10 == 0:
        print(i)
    name = 'fasta/sim' + str(i)
    os.rename('fasta/sim' + str(i) + '.simseq.genome.fa', name)
    with open(name, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(name, 'w') as fout:
        fout.writelines(data[1:])
