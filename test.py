import ghofAtkinson
import matplotlib.pyplot as plt
import numpy as np

'''
TEST 1:

Recreate figure 9 in the paper
R = 100km, Vs = 760, M = 7 to 9, cascadia events
'''

def test1():
    M = [7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0]
    Vs30 = [760] * len(M)
    R = [100] * len(M)
    Faba = [0] * len(M)
    Cascadia = [1] * len(M)
    epistemic = 0
    periods = [5.0, 2.0, 1.0, 0.2]
    spectra, phi, tau = ghofAtkinson.spectra(M, R, Faba, Vs30, Cascadia, epistemic, periods)
    print spectra
    plt.subplot(2,2,1)
    plt.semilogy(M,[spec[0] for spec in spectra])
    plt.ylim([2,100])
    plt.yticks([2,10,20,100], [2,10,20,100])
    plt.subplot(2,2,2)
    plt.semilogy(M,[spec[1] for spec in spectra])
    plt.ylim([10,200])
    plt.yticks([10,20,30,100,200], [10,20,30,100,200])
    plt.subplot(2,2,3)
    plt.semilogy(M,[spec[2] for spec in spectra])
    plt.ylim([20,200])
    plt.yticks([20,30,40,100,200], [20,30,40,100,200])
    plt.subplot(2,2,4)
    plt.semilogy(M,[spec[3] for spec in spectra])
    plt.ylim([30,350])
    plt.yticks([30,100,200,300], [30,100,200,300])
    plt.show()

def test2():
    M = [7.0]
    Vs30 = [760]
    R = [100]
    Faba = [0]
    Cascadia = [1]
    epistemic = 0
    periods = [1.0]
    spectra, phi, tau = ghofAtkinson.spectra(M, R, Faba, Vs30, Cascadia, epistemic, periods)
    print spectra


if __name__ == '__main__':
    test1()
