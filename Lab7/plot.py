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


def color_map(from_file, image_name, title, values):
    x, y, z = load_data(from_file, params=3)

    for i in range(len(list(zip(x,y)))):
        fig = plt.figure(figsize=(10,4.5))
        x_0 = np.unique(x[i])
        y_0 = np.unique(y[i])
        
        z_0 = np.reshape(z[i], (len(x_0), len(y_0))).transpose()
        
        x_0, y_0 = np.meshgrid(x_0, y_0)
        plt.pcolor(x_0, y_0, z_0, cmap='cool', shading="auto")
        plt.title(f"{title}, Q = {values[i]}")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.colorbar()
        plt.savefig(f'{image_name}_{values[i]}.png', dpi=200)


def contour_map(from_file, image_name, title, values): 
    x, y, z = load_data(from_file, params=3)
    i_1 = 50
    j_1 = 55
    
    for i in range(len(list(zip(x,y)))):
        fig = plt.figure(figsize=(10,4.5))
        x_0 = np.unique(x[i])
        y_0 = np.unique(y[i])
        
        z_0 = np.reshape(z[i], (len(x_0), len(y_0))).transpose()

        z_0[0:j_1, 0:i_1] = np.nan # wyeliminowanie brzegu, aby nie zaburzał wyników
        x_0, y_0 = np.meshgrid(x_0, y_0)
        
        plt.contour(x_0, y_0, z_0, cmap='cool', levels=40, linewidths=0.7)
        plt.title(f"{title}, Q = {values[i]}")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.colorbar()
        plt.savefig(f'{image_name}_{values[i]}.png', dpi=200)


Q = [-1000, -4000, 4000]
contour_map("psi.dat", "psi", r"${{ \psi }}(x,y)$", Q)
contour_map("zeta.dat", "zeta", r"${{ \zeta }}(x,y)$", Q)
color_map("u.dat", "u", "u(x,y)", Q)
color_map("v.dat", "v", "v(x,y)", Q)





