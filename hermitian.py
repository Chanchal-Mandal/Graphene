import numpy as np 
from numpy.linalg import eig
import csv
import matplotlib.pyplot as plt 

n = 4 # no of hexagon
m = 4# no of rows
s = (2*n +1)*(m+1) # Total no of sites 
t=1 # Hopping amplitude 

H = np.zeros([s, s], dtype= int)

# Code to find the eigen values of the hermitian system
for i in range(0, s-1) :      # Amplitude due to its preceeding and suceeding element except the boundary connection 
      if i in range  (2*n,s,2*n+1) :
         continue
      H[i, i+1] = 1 
      H[i+1, i] = 1 
      
                  
for j in range(2*n+1, s) :   # Amplitude due to vertical connections sites at even place 
            if  j % 2 != 0 :
             H[j,j-(2*n+1)] = 1
             
for k in range(0, s-(2*n+1)) :  # Amplitude due to vertical connections sites at odd place 
            if k%2 == 0 :
                H[k,k+(2*n+1)] = 1

w,v = np.linalg.eig(H)
z = np.array([w])
#z.T
z= np.sort(z)

with open('Matrix.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(H)
with open('eigenvalue.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(z)
    
"""rows,col = H.shape
index_values = np.arange(1, rows + 1)
H = np.insert(H, 0, index_values, axis=0)
H = np.insert(H, 0, index_values, axis=1)"""

w,v = np.linalg.eig(H)  
print(H)
print(w)

# Code for Density of states :

I = np.array([])
E = np.array([])

for l in np.arange (-4.0,4.0,0.01) :
    I =np.append(I,l)
    
"""for m in range (1,s) :
    x = (w[m])
    E= np.append(E,x)"""
    
plt.hist(w,bins=I)
plt.xlabel('Energy')
plt.ylabel('DOS')
plt.title("Density Of states")
plt.show()

