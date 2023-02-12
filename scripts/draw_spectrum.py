import numpy as np
import matplotlib.pyplot as plt

def main():
    file = open("data/spectrum", "r")
    l_list, cl_list = [], []

    for line in file.readlines():
        l, _, cl = line.split()
        l_list.append(float(l))
        cl_list.append(float(cl))

    l = np.array(l_list)
    cl = np.array(cl_list)

    plt.plot(l, cl)
    plt.show()

if __name__ == "__main__":
    main()