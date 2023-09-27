import csv, json, sys

header = ['Time', 'Pxx', 'Pyy', 'Pzz']

n = 0
gamma1 = 0.0
gamma2 = 0.0

with open(sys.argv[1]) as datfile:
    for line in datfile:
        data = line.split()
        time = float(data[0])
        pxx = float(data[1])
        pyy = float(data[2])
        pzz = float(data[3])
        gamma = float(data[4])
        n = n + 1
        lx = 40.0
        gamma1 = gamma1 + 0.5*lx*(pzz - 0.5*(pxx + pyy))
        gamma2 = gamma2 + gamma
    gamma1 = gamma1/n
    gamma2 = gamma2/n
    gamma2 = gamma2/2
    print(gamma1)
    print(gamma2)