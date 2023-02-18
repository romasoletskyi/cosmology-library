import sys
from astropy.io import fits


def main():
    hdul = fits.open(sys.argv[1])
    l_list, spectrum, errors = hdul[2].data['ELL'], hdul[2].data['D_ELL'], hdul[2].data['ERR']
    hdul.close()

    with open(sys.argv[2], 'w') as file:
        print(len(spectrum), file=file)
        for l, cl, err in zip(l_list, spectrum, errors):
            print(l, cl, err, file=file)


if __name__ == "__main__":
    main()
