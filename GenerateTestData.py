import random

num_files = 100
str_length = 1000

for i in range(num_files):
    f = open('testData/string' + str(i), 'w')
    string = ''.join(random.choice("ACGT") for i in range(str_length))
    f.write(string)