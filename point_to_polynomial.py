from math import floor, ceil
import numpy as np
import matplotlib.pyplot as plt

def gauss(A):
    n = len(A)
 
    for i in range(0, n):
        # Search for maximum in this column
        maxEl = abs(A[i][i])
        maxRow = i
        for k in range(i+1, n):
            if abs(A[k][i]) > maxEl:
                maxEl = abs(A[k][i])
                maxRow = k
 
        # Swap maximum row with current row (column by column)
        for k in range(i, n+1):
            tmp = A[maxRow][k]
            A[maxRow][k] = A[i][k]
            A[i][k] = tmp
 
        # Make all rows below this one 0 in current column
        for k in range(i+1, n):
            c = -A[k][i]/A[i][i]
            for j in range(i, n+1):
                if i == j:
                    A[k][j] = 0
                else:
                    A[k][j] += c * A[i][j]
 
    # Solve equation Ax=b for an upper triangular matrix A
    x = [0 for i in range(n)]
    for i in range(n-1, -1, -1):
        x[i] = A[i][n]/A[i][i]
        for k in range(i-1, -1, -1):
            A[k][n] -= A[k][i] * x[i]
    return x

def point_to_polynomial(coordlist, degree):
    x=[]
    y=[]
    
    for point in coordlist:
        x.append(point[0])
        y.append(point[1])
    
    amount = len(x)
    exponentsx=[0]*(len(x)) 
    for i in range(0,len(x)):
        exponentsx[i] = x[i]
    
    diag = np.zeros((degree+1,degree+2))
    diag[0][0] = len(x)
    
    #loop for num of rows (excluding row 0)
    
    for i in range(1, degree+1):
        count = 0
        #loop for columns
        for j in range(0, degree+1):
            #making exponentsx hold proper powers of x values
            for k in range(0,len(exponentsx)):
                exponentsx[k] = x[k]**(count + i)
            diag[i][j] = sum(exponentsx)
            #making exponentsx hold proper values
            count += 1
    
    #making row 0 hold correct values
    for j in range(0, degree):
        #making exponentsx hold proper powers of x values
        for k in range(0,len(exponentsx)):
            exponentsx[k] = x[k]**(j+1)
        diag[0][j+1] = sum(exponentsx)
        #making exponentsx hold proper values
    
    #making column sum(xy) correct
    for i in range(0,degree+1):
    #compute correct sum
        temp = 0
        #if x is 0, python computes 0^0 as 1, but we still want to just add y
        for j in range(0,len(x)):
            if(x[j] == 0 and i == 0):
                temp += y[j]
            else:
                temp += x[j]**i * y[j]
        diag[i][degree+1] = temp
    
    mat = np.zeros((degree+1,degree+1))
    
    #get the square part of diag to find the determinant (if 0, matrix is singluar)
    for i in range(0,degree+1):
        for j in range(0, degree+1):
            mat[i][j] = diag[i][j]
        
    poly = gauss(diag) 
    v = degree
    #if not singluar, do row reduction to solve the matrix system
    if(np.linalg.det(mat) == 0):
        print("Singluar Matrix")
        print("")
        print("Options:")
        print("1) If the matrix is singular you may want to try changing the degree of your polynomial.")
        print("2) You can also try adding or taking away a point on your data plot")
    return poly[::-1]

def plot_polynomial(coefficients, points):
    xlim = [floor(points[0][0]), ceil(points[-1][0])]
    x = np.linspace(int(xlim[0]), int(xlim[1]), int(xlim[1] - xlim[0]) * 10)
    y = [np.polyval(coefficients, i) for i in x]
    plt.plot(x,y,'b-')
    for point in points:
        plt.plot(point[0], point[1], 'ro')
    plt.show()
    

def function_through_points(points, degree):
    points.sort(key=lambda x: x[0])
    coefficients = point_to_polynomial(points, degree)
    plot_polynomial(coefficients, points)
    print("In case you'd like to plot this in a graphing calculator of your own, the polynomial is as follows:")
    exp = len(coefficients) - 1
    polynomial = ""
    for co in coefficients:
        if co % 1 == 0:
            co = int(co)
        if co < 0 and exp != 0:
            polynomial += f" {co}x^{exp} "
            exp -= 1
        elif co < 0 and exp == 1:
            polynomial += f" {co}x "
            exp -= 1
        elif co < 0:
            polynomial += f" {co}"
            exp -= 1
        elif co == float(0):
            exp -= 1
        elif exp == 0:
            polynomial += f"+ {co} "
            exp -=1
        elif exp == 1:
            polynomial += f"+ {co}x "
            exp -= 1
        else:
            polynomial += f"+ {co}x^{exp} "
            exp -=1
    if polynomial[0] == "+":
        polynomial = polynomial[1::]
    polynomial = polynomial.replace("1x", "x")
    polynomial.strip()
    print(polynomial)
            