import numpy as np
import matplotlib.pyplot as plt

def main():
    file = open("data/HomogenousHistory", "r")
    eta_list, a_list = [], []

    for line in file.readlines():
        eta, a = line.split()
        eta_list.append(float(eta))
        a_list.append(float(a))

    plt.plot(eta_list, a_list)
    plt.show()

if __name__ == "__main__":
    main()