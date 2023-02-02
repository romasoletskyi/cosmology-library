import numpy as np
import matplotlib.pyplot as plt

def main():
    file = open("data/HomogenousHistory", "r")
    eta_list, a_list, xe_list = [], [], []

    for line in file.readlines():
        eta, a, xe = line.split()
        eta_list.append(float(eta))
        a_list.append(float(a))
        xe_list.append(float(xe))

    eta = np.array(eta_list)
    a = np.array(a_list)
    xe = np.array(xe_list)
    z = 1 / a - 1

    plt.plot(-np.log(z), np.log(xe))
    plt.show()

if __name__ == "__main__":
    main()