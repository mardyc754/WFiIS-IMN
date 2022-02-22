import numpy as np
import matplotlib.pyplot as plt


def plot_from_file(filename, plot_title, image_name):
    fig = plt.figure()

    data = np.loadtxt(filename)

    t = data[:,0]
    u = data[:,1]
    z = data[:,2]

    plt.plot(t, u, label="u(t)")
    plt.plot(t, z, label="z(t)")

    plt.xlim(1,100)
    plt.ylim(0, 500)

    l1 = plt.legend()
    plt.grid()

    plt.xlabel('Czas (t)')
    plt.ylabel('Populacja')
    plt.title(plot_title)

    plt.savefig(image_name, dpi=200)


plot_from_file("picard.dat", "Iteracja Picarda", "zad1.png")
plot_from_file("newton.dat", "Iteracja Newtona", "zad2.png")
plot_from_file("nrk2.dat", "Niejawna metoda RK2", "zad3.png")