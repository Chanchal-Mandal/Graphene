import matplotlib.pyplot as plt
import math as mt
import numpy as np
from mpl_toolkits import mplot3d

X = np.array([])
Y= np.array([])
M = np.array([])
P = np.array([]) 


n =20 #No of hexagon in the row
m =20#No of rows
c0 = 1 # interlayer distance
s = (2*n+1)*(m+1)
print(s)

X= np.append(X,0)
Y= np.append(Y,0)
M = np.append(M,np.sqrt(3)/2)
P = np.append(P,1/2) 
Z= np.ones(s)*c0
#Z = np.append(Z,1)


for i in range (1,s) :
     if i in range (4*n+2 ,s,4*n+2) :
          print(i)
          print('xxx')
          x = X[i-1] - n*np.sqrt(3)
          y = Y[i-1] + 2
          m = x + np.sqrt(3)/2
          p = y + 1/2
          M = np.append(M,m)
          P = np.append(P,p)
          X = np.append(X,x)
          Y = np.append(Y,y)
         

     elif i%2 == 0 :
        print(i)
        x = X[i-1] + np.sqrt(3)/2
        y = Y[i-1] + 1/2
        m = x + np.sqrt(3)/2
        p = y + 1/2
        M = np.append(M,m)
        P = np.append(P,p)
        X = np.append(X,x)
        Y = np.append(Y,y)
        

     elif i in range (2*n+1 ,s,2*n+1) :
          print(i)
          x = X[i-1] - n*np.sqrt(3)
          y = Y[i-1] + 1
          m = x + np.sqrt(3)/2
          p = y + 1/2
          M = np.append(M,m)
          P = np.append(P,p)
          X = np.append(X,x)
          Y = np.append(Y,y)
          
        
     else :
          x = X [i-1] + np.sqrt(3)/2
          y=  Y[i-1] - 1/2
          m = x + np.sqrt(3)/2
          p = y + 1/2
          M = np.append(M,m)
          P = np.append(P,p)
          X= np.append(X,x)
          Y= np.append(Y,y)
          

ax = plt.axes(projection ='3d')
print(s)
print(X)
print(Y)

ax.scatter(X,Y,c="red")
ax.scatter(M, P, Z, c="blue")
plt.show()

