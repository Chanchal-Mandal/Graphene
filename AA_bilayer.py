#CODE to plot the AA 2D strucutre of graphene and its Density of States 
import matplotlib.pyplot as plt
from numpy.linalg import eig
import numpy as np
import csv
from mpl_toolkits import mplot3d

X = np.array([])
Y= np.array([])
M = np.array([])
P = np.array([]) 


n = 20 #No of hexagon in the row
m = 20 #No of rows
c0 = 1 # interlayer distance
S = (2*n+1)*(m+1)
print(S)

X= np.append(X,0)
Y= np.append(Y,0)
M = np.append(M,np.sqrt(3)/2)
P = np.append(P,1/2) 
Z= np.ones(S)*c0

H = np.zeros([2*S, 2*S], dtype= int)


for i in range (1,S) :
     if i in range (4*n+2 ,S,4*n+2) :
          # print(i)
          # print('xxx')
          x = X[i-1] - n*np.sqrt(3)
          y = Y[i-1] + 2
          m = x + np.sqrt(3)/2
          p = y + 1/2
          M = np.append(M,m)
          P = np.append(P,p)
          X = np.append(X,x)
          Y = np.append(Y,y)
         

     elif i%2 == 0 :
        #print(i)
        x = X[i-1] + np.sqrt(3)/2
        y = Y[i-1] + 1/2
        m = x + np.sqrt(3)/2
        p = y + 1/2
        M = np.append(M,m)
        P = np.append(P,p)
        X = np.append(X,x)
        Y = np.append(Y,y)
        

     elif i in range (2*n+1 ,S,2*n+1) :
          #print(i)
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
          
p1 = np.array([])
p2 = np.array([])


t_0 = 1 # interlayer hopping
t_1 = 1 # interalayer hopping

# Hamiltonian matrix 
s = int(S/2)

for i in range(0, S) :      # Amplitude due to intra-layer AA hoping
     if i in range  (0,s) :
          H[i,s+i] = 1    
     else: 
          H[i,i-s] = 1
     
for i in range(S, 2*S) :      # Amplitude due to intra-layer AB hoping 
     if i in range  (S,S+s) :
          H[i,s+i] = 1
     else:
          H[i,i-s] = 1

for i in range (0,2*n+1,S):
     pass
     
for i in range (1,s):   # Amplitude of hopping due to AB inter-layer bond 
     # if i in range (2*n+1,5,2*n+1) :
     #      continue
     H[i, S-1+i] = 1 
     H[S-1+1,i]  = 1
     # if i in range (0,5,2*n+1) :
     #     continue
     H[i, S+i] = 1
     H[i, i] = 1
    
for i in range (s,S):  
     # if i in range (2*n+1,5,2*n+1) :
     #      continue
     H[S-1+i, i] = 1
     H[i, S-1+i] = 1 
     # if i in range (0,5,2*n+1) :
     #     continue
     H[i, S+i] = 1
     H[S+i, i] = 1
     

print(H)

w,v = np.linalg.eig(H)
z = (np.array([w]))
z= np.sort(z)
    

with open('AA.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(H)
with open('AAeigen.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(z)
    

plt.hist(w,bins = 100)
plt.xlabel('Energy Eigenvalue')
plt.ylabel('DOS')
plt.title("Density Of states for AA Bilayer Graphene")
#plt.show()

ax = plt.axes(projection ='3d')

# #plt.style.use('seaborn')
plt.scatter(X,Y,c="red")
plt.scatter(M, P,s=20)
plt.title("AB bilayer Graphene 2d structure")

# #plt.scatter(p1,p2,s=20,c="green")
plt.show()
