""" Rename simulated genomes once generated so that the main.py pipeline can use them as inputs """
import os
i = 0
for file in os.listdir('fasta'):
    if 'simseq.genome.fa' in file:
        if i % 10 == 0:
            print(i, "files renamed")
        new_name = 'fasta/sim' + str(i)
        i += 1
        try:
            os.remove(new_name)
        except FileNotFoundError:
            pass
        os.rename('fasta/' + file, new_name)
        with open(new_name, 'r') as fin:
            data = fin.read().splitlines(True)
        with open(new_name, 'w') as fout:
            fout.writelines(data[1:])
