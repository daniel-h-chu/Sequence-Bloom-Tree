import random

num_files = 500
str_length = 100000

for i in range(num_files):
    f = open('test_data/sequence' + str(i), 'w')
    f.write(''.join(random.choice("ACGT") for i in range(str_length)))
    f.close()
