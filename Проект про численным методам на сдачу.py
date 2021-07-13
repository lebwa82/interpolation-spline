import numpy as np
from matplotlib import pyplot as plt


print("enter the number of points")
#a=int(input())
#x = np.linspace(1, 10, a)
#x = np.array([1,2,3, 4, 5,6,7,8,9,10,11,12,13,14,15]) #массив x
#y = np.array([10,9,12, 10, 11,54,1,45,4,45,2,23,4,12,1]) #массив y
#p = np.array([1000,1,0.0001,1,10,1,0.0001,100,0.0001,1,10,1,0.0001,1,1000]) #массив весов
print("enter the point first and then its weight")
#y = np.zeros(a)
x=np.array([0,1,2,3,4,5,6,7,8,9,10])
y=np.array ([1., 1.99495, 1.20545, 1.8703, 3.84104, 5.01165, 6.69188, 8.24921, 7.71055, 7.71338, 9.6174])
p=np.array([100,1,100,1,1,1,1,1,100,1,100])
#p = np.zeros(a)
#for i in range(a):
  #y[i],p[i]=map(float,input().split())
  #y[i] = np.sin(x[i])
  #p[i] = float(input())

m=len(x)-1 #количество точек

R = np.zeros((m+1, m+1))
for i in range(m+1):
  R[i][i] = 1/p[i] #диагональна€ матрица с весами

h = np.zeros(m) #длина интервала [xi,xi+1]
A = np.zeros((m+1,m+1))
H = np.zeros((m+1,m+1))
for i in range(m):
  h[i] = x[i+1] - x[i]

n = np.zeros(m)

v = np.zeros(m)
for i in range(m):
  v[i]=1/h[i]

w = np.zeros(m-1)
for i in range(m-1):
  w[i]=1/h[i]+1/h[i+1]


for i in range(2, m-1):
  A[i][i] = 2*(h[i-1]+h[i])
  A[i][i+1] = h[i]
  A[i][i-1] = h[i-1]
  
  
A[0][0] = 1
A[m][m] = 1
A[1][1] = 2*(h[0]+h[1])
A[1][2] = h[1]
A[m-1][m-1] =  2*(h[m-2]+h[m-1])
A[m-2][m-1] = h[m-2]


for i in range(1,m):
  H[i][i-1]=v[i-1]
  H[i][i]=-1*w[i-1]
  H[i][i+1]=v[i]
  
HT = H.transpose()
HR = 6*np.dot(H,R) #перемножает матрицы
HRHT = np.dot(HR, HT)
RES = A+HRHT
M = np.linalg.solve(RES, 6*np.dot(H,y))


RHT = np.dot(R, HT)
RHTM = np.dot(RHT,M)
z = y - RHTM 




#print(A, '\n', H, '\n', RHT, '\n', HT, '\n', HR, '\n', np.dot(H,y))
#print(m)
#print(z)

# Ћинейна€ зависимость
x_g = []
y_g = []
#np.linspace(0, 11, 50)
for i in range(m):
  v = 50
  t = np.linspace(0, 1, v)
  for j in range(v):
    k=z[i]*(1-t[j])+z[i+1]*t[j]-(h[i]**2)/6*t[j]*(1-t[j])*((2-t[j])*M[i]+(1+t[j])*M[i+1]) 
    y_g.append(k)
    x_g.append(h[i]*t[j]+x[i])


# ѕостроение графика
plt.title("«ависимости: ") # заголовок
plt.xlabel("")         # ось абсцисс
plt.ylabel("y1, y2")    # ось ординат
plt.grid()              # включение отображение сетки
plt.plot(x, y, x_g, y_g)  # построение графика

