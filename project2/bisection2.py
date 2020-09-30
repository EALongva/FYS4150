import numpy as np

n = 10
c = np.ones(n)*2
b = np.ones(n)*(-1)
beta = b**2
m1 = 0; #m1 = 0;
m2 = n-1; #m2 = 4;
eps1 = 10**(-10)
relfeh = 10**(-308) # ??
x = np.zeros(n)

beta[0]=b[0]=0.0;

# Set initial xmin and xmax - theorem???
xmin = c[n-1] - abs(b[n-1]);
xmax = c[n-1] + abs(b[n-1]);
for i in range(n-2,0,-1):
    h = abs(b[i]) + abs(b[i+1]);
    if c[i] + h > xmax:
        xmax = c[i] + h;
    if c[i] - h < xmin:
        xmin = c[i] - h;

wu = np.zeros(n);
x0 = xmax;
for i in range(m1,m2+1):
    x[i] = xmax;
    wu[i] = xmin;
z = 0;
for k in range(m2,m1-1,-1):
    xu = xmin;
    for i in range(k,m1-1,-1):
        if xu < wu[i]:
            xu = wu[i]
    if x0 > x[k]:
        x0 = x[k]
    while x0 - xu >  eps1 and z < 10**5:
        x1 = (xu + x0)/2;
        z += 1
        # Sturms sequence
        a = 0; q = 1; # P0 = 1
        for i in range(0,n):
            if q != 0:
                q = c[i] - x1 -  beta[i]/q;
            else:
                q = c[i] - x1 - beta[i]/relfeh
            if q < 0:
                a += 1;
        if a < k + 1: # redigert
            if a < m1 + 1: # redigert
                xu = wu[m1] = x1;
            else:
                xu = wu[a] = x1; # redigert
                if x[a-1] > x1: # redigert
                    x[a-1] = x1; # redigert
        else:
            x0 = x1;
    x[k] = (x0 + xu)/2

print(x)
print(z)
print(wu)
