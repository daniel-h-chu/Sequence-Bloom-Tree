""" Generates completely random strings as test data for the SBT """
import random

num_files = 500
str_length = 10000

for i in range(num_files):
    if i % 10 == 0:
        print(i + 1, "of", num_files, "strings generated")
    f = open('test_data/sequence' + str(i), 'w')
    f.write(''.join(random.choice("ACGT") for i in range(str_length)))
    f.close()
print("Done")