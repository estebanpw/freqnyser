import os
import sys
import getopt

import csv
import numpy as np
import pandas as pd

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def main(argv):
    input_path = ''
    # path of execution
    output_path = os.path.dirname(os.path.abspath(__file__))
    file = ''

    # Get all arguments given
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["inputfolder=", "outputfolder="])
    except getopt.GetoptError:
        print('multi2simple.py -i <inputfolder> -o <outputfolder>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('multi2simple.py -i <inputfolder> -o <outputfolder>')
            sys.exit()
        elif opt in ("-i", "--inputfolder"):
            if os.path.isabs(arg):
                input_path = arg
            else:
                path = os.path.join(os.getcwd(), arg)
                input_path = arg
        elif opt in ("-o", "--outputfolder"):
            if os.path.isabs(arg):
                output_path = arg
            else:
                output_path = os.path.join(os.getcwd(), arg)

    # columns = ['FILE', 'LENGTH', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    reader = pd.read_csv(os.path.join(input_path, os.listdir(input_path)[0]))
    columns = reader['kmer'].array
    columns = np.insert(columns, 0, 'FILE', axis=0)
    print(columns)

    df = pd.DataFrame(columns=columns)

    for filename in os.listdir(input_path):
        file = pd.read_csv(os.path.join(input_path, filename))
        values = file['rel_freq']
        #values = np.around(values, decimals=8)
        # shifteamos todos los indices 1 a la derecha
        values.index = values.index + 1
        # concatenamos el nombre del fasta en la posicion 0 de la serie
        values = pd.concat([pd.Series([filename]), values])
        # establecemos nuestro indice a los valores
        values.index = columns
        # metemos los valores leidos en el dataframe
        df = df.append(values, ignore_index=True)
    print(df)

    df.to_excel(os.path.join(output_path, 'tabla_freq_rel.xlsx'))
    df.to_csv(os.path.join(output_path, 'tabla_freq_rel.csv'), index=True, index_label='index')


if __name__ == "__main__":
    main(sys.argv[1:])


"""
        tabla[i] = values.to_numpy()
        if primerito:
            kmer = file['kmer']
            primerito = False

        i += 1

    print(type(tabla))

    #kmer.insert(0, 'FILE')
    df = pd.DataFrame(tabla, columns=kmer)
    df.to_excel(os.path.join(output_path, 'tabla.xlsx'))
    df.to_csv(os.path.join(output_path, 'tabla.csv'), index=True, index_label='index')

    '''
    print(x.shape)
    print(z.shape)
    print(z[0].shape)
    print(y.shape)
    '''

    y = np.arange(0, n_sequences, 1)
    #x = np.array(kmer)
    x = np.arange(0,n_kmers,1)
    z = np.array(tabla)
    print(type(z))
    xs, ys = np.meshgrid(x, y)

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(xs, ys, z, rstride=1, cstride=1, cmap='hot', antialiased=True)
    plt.show()
"""
