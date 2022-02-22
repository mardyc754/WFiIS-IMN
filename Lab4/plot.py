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


def plot_map(from_file, plot_name, image_names):
    x, y, z = load_data(from_file, params=3)

    
    labels = [r"${{ {\omega} }}_G = 0.6$", r"${{ {\omega} }}_G = 1.0$"]
    for i in range(2):
        fig = plt.figure()
        x_0 = np.unique(x[i])
        y_0 = np.unique(y[i])
        values = np.reshape(z[i], (len(x_0), len(y_0))).transpose()

        x_0, y_0 = np.meshgrid(x_0, y_0)
        plt.pcolor(x_0, y_0, values, cmap='bwr', shading="auto")
        plt.title(f"{plot_name}, {labels[i]}")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.xlim(0,15)
        plt.ylim(0,10)
        plt.colorbar()
        plt.savefig(f'{image_names[i]}.png', dpi=200)


def line_plot(from_file, plot_name, image_name, labels, y_max):
    x, y = load_data(from_file)
    
    fig = plt.figure()
    for i in range(len(list(zip(x,y)))):
        plt.plot(x[i], y[i], label=labels[i])

    l1 = plt.legend()
    plt.grid()
    plt.xscale('log')
    plt.title(f"{plot_name} - S(it)")
    plt.xlabel("Nr iteracji")
    plt.ylabel("S")
    plt.xlim(1,5e4)
    plt.ylim(0,y_max)
    
    plt.savefig(f'{image_name}.png', dpi=200)


names = ["relaksacja_globalna_0_6", "relaksacja_globalna_1_0",
                "relaksacja_globalna_err_0_6", "relaksacja_globalna_err_1.0"]

plot_map("relaksacja_globalna_V_x_y.dat", "V(x,y)", names[:2])
plot_map("relaksacja_globalna_d_x_y.dat", r"${{ {\delta} }}(x,y)$", names[2:])


labels_global = [r"${{ {\omega} }}_G = 0.6$", r"${{ {\omega} }}_G = 1.0$"]

line_plot("relaksacja_globalna_S_t.dat", "Relaksacja globalna", "relaksacja_globalna_S", labels_global, 5000)

labels_local = [r"${{ {\omega} }}_L = 1.0$", r"${{ {\omega} }}_L = 1.4$",
                r"${{ {\omega} }}_L = 1.8$", r"${{ {\omega} }}_L = 1.9$"]

line_plot("relaksacja_lokalna_S_t.dat", "Relaksacja lokalna", "relaksacja_lokalna_S", labels_local, 4000)
