import numpy as np
from io import StringIO 
import matplotlib.pyplot as plt


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


# metoda Eulera
fig = plt.figure()
t, y = load_data("euler.dat")

labels = ["t = 0.01", "t = 0.1", "t = 1.0"]
    
for i in range(len(t)):
    plt.plot(t[i], y[i], label=labels[i])

plt.plot(t[0], np.exp(-t[0]), 'k-', lw=2, label=r'$e^{{ {\lambda} t}}$')

plt.xlim(0,5)
plt.ylim(0,1)

l1 = plt.legend()
plt.grid()

plt.xlabel('t')
plt.ylabel('y')
plt.title('Euler')

plt.savefig("zad1a.eps")

# blad metody Eulera
fig = plt.figure()
t, err = load_data("euler_err.dat")

for i in range(len(t)):
    if(i>0):
        plt.plot(t[i][1:], err[i][1:], label=labels[i])
    else:
        plt.plot(t[i], err[i], label=labels[i])

plt.xlim(0,5)
plt.ylim(1e-5, 1)
plt.yscale('log')

l1 = plt.legend()
plt.grid()

plt.xlabel('t')
plt.ylabel(r'$\delta$')
plt.title('Euler')

plt.savefig("zad1b.eps")


# metoda RK2
fig = plt.figure()
t, y = load_data("rk2.dat")

for i in range(len(t)):
    plt.plot(t[i], y[i], label=labels[i])

plt.plot(t[0], np.exp(-t[0]), 'k-', lw=2, label=r'$e^{{ {\lambda} t}}$')
plt.xlim(0,5)
plt.ylim(0,1)

l1 = plt.legend()
plt.grid()

plt.xlabel('t')
plt.ylabel('y')
plt.title('RK2')

plt.savefig("zad2a.eps")


# blad metody RK2
fig = plt.figure()
t, err = load_data("rk2_err.dat")

for i in range(len(t)):
    if(i>0):
        plt.plot(t[i][1:], err[i][1:], label=labels[i])
    else:
        plt.plot(t[i], err[i], label=labels[i])


plt.xlim(0,5)
plt.ylim(1e-7, 1)
plt.yscale('log')

l1 = plt.legend()
plt.grid()

plt.xlabel('t')
plt.ylabel(r'$\delta$')
plt.title('RK2')

plt.savefig("zad2b.eps")


# metoda RK4
fig = plt.figure()
t, y = load_data("rk4.dat")

for i in range(len(t)):
    plt.plot(t[i], y[i], label=labels[i])

plt.plot(t[0], np.exp(-t[0]), 'k-', lw=2, label=r'$e^{{ {\lambda} t}}$')
plt.xlim(0,5)
plt.ylim(0,1)

l1 = plt.legend()
plt.grid()

plt.xlabel('t')
plt.ylabel('y')
plt.title('RK4')

plt.savefig("zad3a.eps")


# blad metody RK4
fig = plt.figure()
t, err = load_data("rk4_err.dat")

for i in range(len(t)):
    if(i>0):
        plt.plot(t[i][1:], err[i][1:], label=labels[i])
    else:
        plt.plot(t[i], err[i], label=labels[i])


plt.xlim(0,5)
plt.ylim(1e-13, 0.1)
plt.yscale('log')

l1 = plt.legend()
plt.grid()

plt.xlabel('t')
plt.ylabel(r'$\delta$')
plt.title('RK4')

plt.savefig("zad3b.eps")


labels = [r"${{ {\omega} }}_V = 0.5 {{ {\omega} }}_0$",
          r"${{ {\omega} }}_V = 0.8 {{ {\omega} }}_0$",
          r"${{ {\omega} }}_V = 1.0 {{ {\omega} }}_0$",
          r"${{ {\omega} }}_V = 1.2 {{ {\omega} }}_0$"]


# RRZ2 - wykresy Q(t)
fig = plt.figure()
t, Q = load_data("rlc_q.dat")

for i in range(len(t)):
    plt.plot(t[i], Q[i], label=labels[i])

l1 = plt.legend()
plt.grid()

plt.xlabel('t')
plt.ylabel('Q')
plt.title('z4_RLC_Q')

plt.savefig("zad4a.eps")


# RRZ2 - wykresy I(t)
fig = plt.figure()
t, I = load_data("rlc_i.dat")

for i in range(len(t)):
    plt.plot(t[i], I[i], label=labels[i])

l1 = plt.legend()
plt.grid()

plt.xlabel('t')
plt.ylabel('I')
plt.title('z4_RLC_I')

plt.savefig("zad4b.eps")