
import numpy as np
import matplotlib.pyplot as plt

def main():

    p1 = [0, 0, 0, 0, 0, 0]
    p2 = [5, 5, 5, 5, 5, 5]
    p3 = [20, 20, 20, 20, 20, 20]
    p4 = [30, 30, 30, 30, 30, 30]
    p5 = [40, 40, 40, 40, 40, 40]

    t = [300, 600, 900, 1200, 1500, 1800]
    tp = [0, 300, 600, 900, 1200, 1500, 1800]

    tns = [300, 600, 900, 1200, 1500]
    pns = [10, 10, 10, 10, 10]

    ts = [1800]
    ps = [10]

    plt.plot(p1, t, 'bo', p2, t, 'bo', pns, tns, 'bo', ps, ts, 'b^', p3, t, 'ro', p4, t, 'ro', p5, t, 'ro')
    axes = plt.gca()
    plt.xlabel('Pressure /GPa')
    plt.ylabel('Temperature /K')
    plt.ylim([0, 2000])
    plt.xlim([0, 45])
    axes.set_xticks([0, 5, 10, 20, 30, 40])
    axes.set_yticks(tp)
    plt.show()


main()
