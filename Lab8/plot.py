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



def color_map(from_file, image_name, values):
    x, y, z = load_data(from_file, params=3)
    
    titles = [r'$v^{x}(x,y)$', r'$v^{y}(x,y)$']
    for i in range(len(list(zip(x,y)))):
        fig = plt.figure(figsize=(20,4.5))
        x_0 = np.unique(x[i])
        y_0 = np.unique(y[i])        
        z_0 = np.reshape(z[i], (len(x_0), len(y_0))).transpose()
        
        x_0, y_0 = np.meshgrid(x_0, y_0)
        plt.pcolor(x_0, y_0, z_0, cmap='gnuplot', shading="auto")
        plt.title(titles[i])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.colorbar()
        plt.savefig(f'{image_name}_{values[i]}.png', dpi=200)


def line_plot(from_file, image_name, labels):
    x, y, z = load_data(from_file, params=3)
    
    fig = plt.figure()
    for i in range(len(list(zip(x,y)))):
        plt.plot(x[i], y[i], label=labels[i])
        plt.plot(x[i], z[i], label=labels[i+2])

    l1 = plt.legend(loc="upper right")
    plt.title(r"$c(t), x_{sr}(t)$")
    plt.xlabel("t")
    plt.ylabel(r"$c(t), x_{sr}(t)$")
    plt.grid()
    
    plt.savefig(f'{image_name}.png', dpi=200)



def colormap_group(filename, image_name, figure_title):
    fig = plt.figure( figsize=(20, 12))
    x, y, z = load_data(filename, params=3)
    
    it = [0, 2000,4000,6000,8000,10000]
    for i in range(len(list(zip(x,y)))):
        x_0 = np.unique(x[i])
        y_0 = np.unique(y[i])

        z_0 = np.reshape(z[i], (len(x_0), len(y_0))).transpose()
        
        axs = fig.add_subplot(3, 2, i+1)
        x_0, y_0 = np.meshgrid(x_0, y_0)
        im = axs.pcolor(x_0, y_0, z_0, cmap='gnuplot', shading="auto", vmin=0)
        plt.subplots_adjust(left=0.05,
                    bottom=0.1, 
                    right=1, 
                    top=0.9, 
                    wspace=0, 
                    hspace=0.5)
        axs.set(title=f"it = {it[i]}", xlabel="x", ylabel="y")
        fig.colorbar(im, ax=axs)

    fig.suptitle(figure_title, fontsize=16)    
    plt.savefig(f"{image_name}", dpi=200)



color_map("v.dat", "v", ["x", "y"])
line_plot('c_x_sr.dat', 'c_x_sr', [r'$x_{sr}(D=0)$', r'$x_{sr}(D=0.1)$', 'c(D=0)', 'c(D=0.1)'])

colormap_group("u_D_0.dat", "u_D_0.png", "u(x,y), D = 0")
colormap_group("u_D_01.dat", "u_D_01.png", "u(x,y), D = 0.1")




