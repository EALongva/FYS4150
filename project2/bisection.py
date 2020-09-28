import numpy as np
"""
n = 10
A = np.zeros((n,n))
d = np.ones(n)*2
a =
a = -1
A[0,0] = d[0];
A[n-1,n-1] = d[n-1];
A[0,1] = A[n-1,n-2] = a;
for i in range(1,n-1):
    A[i,i] = d[i]
    A[i,i+1] = a
    A[i,i-1] = a

for i in range(1,n-1):
    xmin = d[i] - abs(a)
    xmax = d[i] + abs(a)
    for i in range(1,n):
"""
n = 5
c = np.ones(n)*2
b = np.ones(n)*(-1)
#c = np.zeros(n)
#b = np.zeros(n)
#c[0] = 1; c[1] = 49;
#b[1] = 7;
beta = b**2
m1 = 0; #m1 = 0;
m2 = n-1; #m2 = 4;
eps1 = 10**(-10)
relfeh = 2**(-20)
x = np.zeros(int(m2-m1+1))

beta[0]=b[0]=0.0;

# Set initial xmin and xmax - theorem???
xmin = c[n-1] - abs(b[n-1]);
xmax = c[n-1] + abs(b[n-1]);
for i in range(0,n-1):
    h = abs(b[i]) + abs(b[i+1]);
    if c[i] + h > xmax:
        xmax = c[i] + h;
    if c[i] - h < xmin:
        xmin = c[i] - h;

# ????
if xmin + xmax > 0:
    eps2 = relfeh*xmax;
else:
    eps2 = relfeh*(-xmin);
# ?????

if eps1 <= 0:
    eps1 = eps2;
eps2 = 0.5*eps1 + 7*eps2 # ??????

wu = np.zeros(m2-m1+1);
x0 = xmax;
for i in range(m1,m2+1):
    x[i] = xmax;
    wu[i] = xmin;
z = 0;
for k in range(m1,m2,1):
    xu = xmin;
    for i in range(k,m1,-1):
        if xu < wu[i]:
            xu = wu[i]
    if x0 > x[k]:
        x0 = x[k]
    while x0 - xu > 2*relfeh*(abs(xu)+abs(x0)) + eps1 and z < 10**5:
        x1 = (xu + x0)/2;
        z += 1
        # Sturms sequence
        a = 0; q = 1;
        for i in range(0,n-1):
            if q != 0:
                q = c[i] - x1 -  beta[i]/q;
            else:
                q = c[i] - x1 - abs(b[i]/relfeh)
            if q < 0:
                a += 1;
        if a < k:
            if a < m1:
                xu = wu[m1] = x1;
            else:
                xu = wu[a+1] = x1;
                if x[a] > x1:
                    x[a] = x1;
        else:
            x0 = x1;
    x[k] = (x0 + xu)/2

print(x)
print(z)
