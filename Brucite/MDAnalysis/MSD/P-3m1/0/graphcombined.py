

import numpy as np
import matplotlib.pyplot as plt

def main():
    fileMg = 'p3m1_combmsd_Mg_p0.txt'
    fileO = 'p3m1_combmsd_O_p0.txt'
    fileH = 'p3m1_combmsd_H_p0.txt'

    Mg_MSD = []
    O_MSD = []
    H_MSD = []

    fMg = open(fileMg, 'r')
    for line in fMg:
        Mg_MSD.append(float(line))
    fMg.close()

    fO = open(fileO, 'r')
    for line in fO:
        O_MSD.append(float(line))
    fO.close()

    fH = open(fileH, 'r')
    for line in fH:
        H_MSD.append(float(line))
    fH.close()

    plt.plot(Mg_MSD, 'r')
    plt.plot(O_MSD, 'b')
    plt.plot(H_MSD, 'g')
    plt.xlabel('Step')
    plt.ylabel('MSD /Angstrom^2')
    plt.title('MSD over varying temperature')
    plt.legend(('Mg','O','H'), loc='upper left')
    plt.savefig('p3m1_combined_msd_p0.png')

main()
