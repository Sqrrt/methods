#Для уравнений


#половинное деление
def f(x):
  return x**3+2*x**2+7*x-3
a,b=map(float,input().split())
x=0
e=float(input())
if f(a)/abs(f(a))==f(b)/abs(f(b)):
  print(’error’)
elif f(a)==0:
  print(’x=’,a)
elif f(b)==0:
  print(’x=’,b)
if b<a:
  a,b=b,a
while abs(b-a)>e:
  x=(b+a)/2
  if f(x)/abs(f(x))== -1:
    a=x
  else:
    b=x
print(’x=’, x)


#метод хорд

def f(x):
  return x**3-3*x=4
a,b=map(int,input().split())
e=float(input())
if f(a)==0:
  print(’x=’,a)
elif f(b)==0:
  print(’x=’,b)
if f(a)/abs(f(a))==f(b)/abs(f(b)):
  print(’error’)
else:
  if f(a)>0:
  while f(b)+e<0:
    b=f(a)*(b-a)/(f(a)-f(b))+a
    print(’x=’,b)
  else:
  while f(a)+e<0:
    a=f(b)*(a-b)/(f(b)-f(a))+b
  print(’x=’,a)

#метод секущих касательных
import math
def f(x):
  return x**3-3*x+4
a,b=map(float,input().split())
e=float(input())
x=0
if f(a)/abs(f(a))==f(b)/abs(f(b)):
  print(’error’)
elif f(a)==0:
  print(’x=’,a)
elif f(b)==0:
  print(’x=’,b)
while (f(b)/abs(f(b))!=f(a)/abs(f(a))):
  if a-b>0:
    b=(a+b)/2
  if b==0:
    b+=1
  else:
    a=(a+b)/2
     if a==0:
      a+=1
print(a,b)
while abs(f(b)-f(a))>e:
  x=(f(a)*b-f(b)*a)/(f(b)-f(a))
  a=b
  b=x
print(’x=’,x)

#метод Ньютона

import math
from math import *
def f(x):
  return x**3+2*x**2+7*x-3
from scipy.misc import derivative
def j(x):
  return derivative(f, x, dx=10(−10))
a,b=map(int,input().split())
e=float(input())
if f(a)/abs(f(a))==f(b)/abs(f(b)):
  print(’error’)
if f(a)==0:
  print(’x=’,a)
elif f(b)==0:
  print(’x=’,b)
d=max(a,b)
while (f(d)-e)/abs(f(d)-e)==1:
  d=d-f(d)/j(d)
print(’x=’, d)


#Для нелинейных систем


#метод простых итераций
import math as mh
import numpy as np
def f(x):
  f = np.zeros([n])
  f[0] = 2*x[0]**2 -x[0]*x[1] - 5*x[0]+1
  f[1] = x[0]+3*(mh.log10(x[0]))-x[1]**2
  return f
def g(x):
  g = np.zeros([n])
  g[0] = np.sqrt(0.5*(x[0]*x[1] + 5*x[0]-1))
  g[1] = np.sqrt(+x[0]+ 3*(mh.log10(x[0])))
  return g
n = 2
r = 100
b = float(input())
def simple(f, x, b):
  x1 = [0,0]
  f0 = f(x)
  for i in range(r):
    if (max(abs(x1[i]-x[i]) for i in range(n)) < b):
      return x, i
    x1 = x
    x = f(x)
    print(x)
    f0 = f(x)
  print('Iterations >= 100')
x0 = [3.5,2.2]
x, res = simple(g, x0, b)
print(x)
print(res)

#метод Зейделя

import math as mh
import numpy as np
def g1(x):
f = np.sqrt(0.5*(x[0]*x[1] + 5*x[0]-1))
  return f
def g2(x):
  f = np.sqrt(+x[0]+ 3*(mh.log10(x[0])))
  return f
def f(x):
  f = np.zeros([n])
  f[0] = 2*x[0]**2 -x[0]*x[1] - 5*x[0]+1
  f[1] = x[0]+3*(mh.log10(x[0]))+x[1]**2
  return f
def g(x):
  g = np.zeros([n])
  g[0] = np.sqrt(0.5*(x[0]*x[1] + 5*x[0]-1))
  g[1] = np.sqrt(+x[0]+ 3*(mh.log10(x[0])))
  return g
n = 2
r = 100
b = float(input())
def simple(f, x, b):
  x1 = [0,0]
  f0 = f(x)
  for i in range(r):
    if (max(abs(x1[i]-x[i]) for i in range(n)) < b):
    return x, i
    for j in range (n):
      if (j == 0):
        x[j] = g1(x)
      if(j == 1):
        x[j] = g2(x)
    x1 = x
    x = f(x)
    print(x)
    f0 = f(x)
  print('Iterations >= 100')
x0 = [3.5,2.2]
x, res = simple(g, x0, b)
print(x)

#метод Ньютона

import math as mh
import numpy as np
def f(x):
  f = np.zeros([n])
  f[0] = 2*x[0]**2 -x[0]*x[1] - 5*x[0]+1
  f[1] = x[0]+3*(mh.log10(x[0]))-x[1]**2
  return f
def g(x):
  g = np.zeros([n])
  g[0] = np.sqrt(0.5*(x[0]*x[1] + 5*x[0]-1))
  g[1] = np.sqrt(+x[0]+ 3*(mh.log10(x[0])))
  return g
n = 2
r = 100
b = float(input())
def simple(f, x, b):
  x1 = [0,0]
  f0 = f(x)
  for i in range(r):
    if (max(abs(x1[i]-x[i]) for i in range(n)) < b):
      return x, i
    x1 = x
    x = f(x)
    print(x)
    f0 = f(x)
  print('Iterations >= 100')
x0 = [3.5,2.2]
x, res = simple(g, x0, b)
print(x)
print(res)
