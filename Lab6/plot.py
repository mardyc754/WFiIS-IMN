import numpy as np
import matplotlib.pyplot as plt
from io import StringIO 


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


def plot_map(from_file, image_name, values, val_min, val_max, map_tiltes):
    x, y, z = load_data(from_file, params=3)

    
    for i in range(len(list(zip(x,y)))):
        fig = plt.figure()
        x_0 = np.unique(x[i])
        y_0 = np.unique(y[i])
        z_0 = np.reshape(z[i], (len(x_0), len(y_0)))

        x_0, y_0 = np.meshgrid(x_0, y_0)
        plt.pcolor(x_0, y_0, z_0, cmap='bwr', shading="auto", vmin=val_min, vmax=val_max)
        plt.title(map_tiltes[i])
        plt.xlabel("i")
        plt.ylabel("j")
        plt.colorbar()
        plt.savefig(f'{image_name}_{values[i]}.png', dpi=200)


titles_n_xy = [r"$n_{x} = n_{y} = 50$", r"$n_{x} = n_{y} = 100$", r"$n_{x} = n_{y} = 200$"]
plot_map("zad_5.dat", "nx_ny", [50, 100, 200], -10, 10, titles_n_xy)

titles_epsilon = [r"${{ {\varepsilon} }}_{{ 1 }} =  1, {{ {\varepsilon} }}_{{ 2 }} =  1$",
                  r"${{ {\varepsilon} }}_{{ 1 }} =  1, {{ {\varepsilon} }}_{{ 2 }} =  2$",
                  r"${{ {\varepsilon} }}_{{ 1 }} =  1, {{ {\varepsilon} }}_{{ 2 }} =  10$"]
plot_map("zad_6.dat", "epsilon", ["1_1", "1_2", "1_10"], -0.7, 0.7, titles_epsilon)
