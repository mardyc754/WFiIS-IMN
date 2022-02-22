import numpy as np
import matplotlib.pyplot as plt
from io import StringIO 


plt.rcParams.update({'font.size': 22})


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


def colormap_group(filename, image_name, figure_title):
    fig = plt.figure(figsize=(35,20))
    x, y, z = load_data(filename, params=3)
    
    it = [100,200,500,1000,2000]
    for i in range(len(list(zip(x,y)))):
        x_0 = np.unique(x[i])
        y_0 = np.unique(y[i])

        z_0 = np.reshape(z[i], (len(x_0), len(y_0))).transpose()
        
        axs = fig.add_subplot(2, 3, i+1)
        x_0, y_0 = np.meshgrid(x_0, y_0)
        im = axs.pcolor(x_0, y_0, z_0, cmap='gnuplot')
        plt.subplots_adjust(left=0.05,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.2, 
                    hspace=0.2)
        axs.set(title=f"it = {it[i]}", xlabel="x", ylabel="y")
        fig.colorbar(im, ax=axs)

    fig.suptitle(figure_title, fontsize=30)    
    plt.savefig(f"{image_name}", dpi=200)


colormap_group("T.dat", "T.png", "T(x,y)")
colormap_group("nabla_2_T.dat", "nabla_2_T.png", r"${{\nabla}}^{{2}}(x,y)$")
