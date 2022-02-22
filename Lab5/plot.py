import numpy as np
import matplotlib.pyplot as plt
from io import StringIO 


def load_data(filename):
    with open(filename) as f:
        text_data = f.read()

    text_data = text_data.split('\n\n\n')

    data = []
    x = []
    y = []
    for d in text_data:
        if str(d):
            data.append(np.loadtxt(StringIO(str(d))))
    
    for i in range(len(data)):
        x.append(data[i][:,0])
        y.append(data[i][:,1])

    return x, y


def load_data(filename, params=2):
    with open(filename) as f:
        text_data = f.read()

    text_data = text_data.split('\n\n\n')

    data = []
    for d in text_data:
        if str(d):
            data.append(np.loadtxt(StringIO(str(d))))
    
    x = []
    y = []

    for i in range(len(data)):
        x.append(data[i][:,0])
        y.append(data[i][:,1])

    if params == 3:
        z = []
        for i in range(len(data)):
            z.append(data[i][:,2])
        
        return x, y, z
    return x, y


def plot_map(from_file, image_name):
    x, y, z = load_data(from_file, params=3)

    k = [16,8,4,2,1]
    labels = ["V(x,y), k = 16", "V(x,y), k = 8",
                "V(x,y), k = 4", "V(x,y), k = 2", "V(x,y), k = 1",]
    for i in range(len(list(zip(x,y)))):
        fig = plt.figure()
        x_0 = np.unique(x[i])
        y_0 = np.unique(y[i])
        values = np.reshape(z[i], (len(x_0), len(y_0))).transpose()
        #x_0 = np.arange(0,151)/10
        #y_0 = np.arange(0,101)/10

        x_0, y_0 = np.meshgrid(x_0, y_0)
        plt.pcolor(x_0, y_0, values, cmap='bwr', shading="auto")
        plt.title(labels[i])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.colorbar()
        plt.savefig(f'{image_name}_{k[i]}.png', dpi=200)



def line_plot(from_file, image_name, labels):
    x, y = load_data(from_file)
    
    fig = plt.figure()
    for i in range(len(list(zip(x,y)))):
        plt.plot(x[i], y[i], label=labels[i])

    l1 = plt.legend()
    plt.grid()
    plt.title("S(it)")
    plt.xlabel("Nr iteracji")
    plt.ylabel("S")
    plt.savefig(f'{image_name}.png', dpi=200)



plot_map("V.dat", "V")

labels = ["k = 16", "k = 8", "k = 4", "k = 2", "k = 1"]

line_plot("S.dat", "S",labels)
