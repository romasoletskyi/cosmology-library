import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
    file = open(sys.argv[1], "r")
    l_list, cl_list = [], []

    for line in file.readlines():
        l, cl = line.split()
        l_list.append(float(l))
        cl_list.append(float(cl))

    l = np.array(l_list)
    cl = np.array(cl_list)

    plt.plot(l, cl, 'k--')
    file.close()

    file = open(sys.argv[2], "r")
    l_list, cl_list, err_list = [], [], []

    for line in file.readlines()[1:]:
        l, cl, err = line.split()
        l_list.append(float(l))
        cl_list.append(float(cl))
        err_list.append(float(err))

    l_exp = np.array(l_list)
    cl_exp = np.array(cl_list)
    err_exp = np.array(err_list)
    file.close()

    plt.errorbar(l_exp, cl_exp, yerr=err_exp, fmt='ro', markersize=2)
    plt.show()

if __name__ == "__main__":
    main()