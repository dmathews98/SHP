'''
Graphs temperature flucuation for MD
'''

import numpy as np
import matplotlib.pyplot as plt

def read_in(filename):
    filein = open(filename, 'r')
    num_lines = sum(1 for line in open(filename))

    temps = np.zeros(num_lines)
    line_num = 0
    for line in filein:
        val = line.split()
        temps[line_num] = float(val[5])
        line_num += 1

    filein.close()

    return temps

def main():
    temps3 = read_in('temps3')
    temps1 = read_in('temps1')
    temps10 = read_in('temps10')

    plt.figure(1)
    plt.plot(temps3, 'r', label = 'SMASS=3')
    plt.plot(temps1, 'b', label = 'SMASS=1')
    plt.plot(temps10, 'g', label = 'SMASS=10')
    plt.title('Temperature variation for varying SMASS value')
    plt.xlabel('Step')
    plt.ylabel('Temperature /K')
    plt.legend(loc='upper left')

    plt.show()

main()
