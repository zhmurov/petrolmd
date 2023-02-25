import csv, json

header = ['Time', 'Pxx', 'Pyy', 'Pzz']

n = 0
gamma = 0.0

with open('nvt.pressure.xvg') as datfile:
    for line in datfile:
        data = line.split()
        time = float(data[0])
        pxx = float(data[1])
        pyy = float(data[2])
        pzz = float(data[3])
        n = n + 1
        lx = 40.0
        gamma = gamma + 0.5*lx*(pxx - 0.5*(pyy + pzz))
    gamma = gamma/n
    print(gamma)