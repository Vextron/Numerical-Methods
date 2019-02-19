from sympy import symbols, Matrix, pi, cos, exp, sin

#Definir as variáveis a usar

fvars = symbols('x1, x2, x3')
x1, x2, x3 = fvars

#Definir sistema de equações

x = 3*x1 - cos(x2*x3) - 1/2
y = x1**2 - 81*(x2 + 0.1)**2 + sin(x3) + 1.06
z = exp(-1*x1*x2) + 20*x3 + (10*pi - 3)/3

#Definição do vetor de funções e da respetiva jacobiana

F = Matrix([x, y, z])
J = F.jacobian(fvars)

#Método de Newton

def newton(X, Tol, N, F, J):
    
    k = 0
        
    while k < N:
        
        Z = dict(zip(fvars, X))

        Jx = J.evalf(subs=Z)
        Fx = F.evalf(subs=Z)
        
        Y = -1*Jx.inv("LU")*Fx
        
        X += Y
     
        if Y.norm() < Tol:
            
            return X
        
        k += 1
        
    return "Maximum iterations reached"

X = Matrix([0.1, 0.1, -0.1])

newton(X, 9e-10, 6, F, J)