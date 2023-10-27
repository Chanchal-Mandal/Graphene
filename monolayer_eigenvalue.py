import numpy as np 
from numpy.linalg import eig
import csv
import matplotlib.pyplot as plt 

n = 31 # no of hexagon
m = 31 # no of rows
s = (2*n +1)*(m+1) # Total no of sites 
t=1 # Hopping amplitude 
print(s)
H = np.zeros([s, s], dtype= int)

# Code to find the eigen values of this hamiltonian system
for i in range(0, s-1) :      # Amplitude due to its preceeding and suceeding element except the boundary connection 
      if i in range (2*n+1,s,2*n+1) :
          continue
      H[i, i+1] = 1 
      if i in range (1,s,2*n+1) :
         continue
      H[i+1, i] = 1 
      
                
for j in range(2*n+1, s) :   # Amplitude due to vertical connections sites at even place 
            if  j % 2 != 0 :
             H[j,j-(2*n+1)] = 1
             
for k in range(0, s-(2*n+1)) :  # Amplitude due to vertical connections sites at odd place 
            if k%2 == 0 :
                H[k,k+(2*n+1)] = 1

print(H)
w,v = np.linalg.eig(H)
z = np.array([w])
#z.T
z= np.sort(z)

# To write the CSV file for Hamilatonian Matrix and its eigen value 
with open('Alayer.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(H)
with open('Aeigenvalue.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(z)

    
f = plt.figure()


plt.hist(w,bins=100)
plt.xlabel('Energy Eigenvalue')
plt.ylabel('DOS')
plt.title("Density Of states for Monolayer Graphene")
plt.show()

