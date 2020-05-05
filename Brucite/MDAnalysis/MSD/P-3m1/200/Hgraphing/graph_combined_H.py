'''
Graphs H distributions over varying pressures
'''

import numpy as np
import matplotlib.pyplot as plt

def main():
    t300 = []
    t600 = []
    t900 = []
    t1200 = []
    t1500 = []
    t1800 = []
    time = []

    filet300 = open('p3m1_msd_H_p200_t300.txt', 'r')
    for line in filet300:
        t300.append(float(line))
    filet300.close()

    filet600 = open('p3m1_msd_H_p200_t600.txt', 'r')
    for line in filet600:
        t600.append(float(line))
    filet600.close()

    filet900 = open('p3m1_msd_H_p200_t900.txt', 'r')
    for line in filet900:
        t900.append(float(line))
    filet900.close()

    filet1200 = open('p3m1_msd_H_p200_t1200.txt', 'r')
    for line in filet1200:
        t1200.append(float(line))
    filet1200.close()

    filet1500 = open('p3m1_msd_H_p200_t1500.txt', 'r')
    for line in filet1500:
        t1500.append(float(line))
    filet1500.close()

    filet1800 = open('p3m1_msd_H_p200_t1800.txt', 'r')
    for line in filet1800:
        t1800.append(float(line))
    filet1800.close()

    step = 10
    for i in range(len(t300)):
        t = i*step*(0.5e-3)
        time.append(t)

    plt.plot(time, t300)
    plt.plot(time, t600)
    plt.plot(time, t900)
    plt.plot(time, t1200)
    plt.plot(time, t1500)
    plt.plot(time, t1800)
    plt.xlabel('Time (ps)')
    plt.ylabel(r'MSD ($\AA^2$)')
    plt.title('MSDs for p-3m1 for Hydogen at different temperatures and 200kbar')
    plt.legend(('300K', '600K', '900K', '1200K', '1500K', '1800K'), loc='upper left')
    plt.savefig('p3m1_msd_H_p200.png')
    plt.show()

main()
