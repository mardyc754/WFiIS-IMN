import matplotlib.pyplot as plt
import numpy as np


def plot(cell, x_axes, y_axes, x_title, y_title, plot_title):
    cell.plot(x_axes[0], y_axes[0], label=r"$tol = 10^{-2}$")
    cell.plot(x_axes[1], y_axes[1], label=r"$tol = 10^{-5}$")

    l1 = cell.legend(loc='upper right')
    cell.grid()

    cell.set(xlabel=x_title, ylabel=y_title, title=plot_title)



def plot_group(filename_start):
    fig, axs = plt.subplots(2, 2, figsize=(10,10))
    #plt.subplots_adjust(
    #                    wspace = 1,   # the amount of width reserved for blank space between subplots
    #                    hspace = 1)
    data_10_2 = np.loadtxt(f"{filename_start}_10_2.dat")
    data_10_5 = np.loadtxt(f"{filename_start}_10_5.dat")

    t_10_2 = data_10_2[:,0]
    t_10_5 = data_10_5[:,0]

    dt_10_2 = data_10_2[:,1]
    dt_10_5 = data_10_5[:,1]

    x_10_2 = data_10_2[:,2]
    x_10_5 = data_10_5[:,2]

    v_10_2 = data_10_2[:,3]
    v_10_5 = data_10_5[:,3]

    plot(axs[0,0], [t_10_2, t_10_5], [x_10_2, x_10_5], 
        "t", "x", "x(t)")
    plot(axs[0,1], [t_10_2, t_10_5], [v_10_2, v_10_5], 
        "t", "v", "v(t)")
    plot(axs[1,0], [t_10_2, t_10_5], [dt_10_2, dt_10_5], "t", 
            r"${ \Delta }t$", r"${ \Delta }t(t)$")
    plot(axs[1,1], [x_10_2, x_10_5], [v_10_2, v_10_5], 
        "x", "v", "v(x)")
    
    plt.savefig(f"{filename_start}.png")


plot_group("metoda_trapezow")
plot_group("rk2")