import matplotlib.pyplot as plt
from sympy import Symbol, sympify, diff, Subs


'''
Metodo Steffensen (SM)
Tomado del Artículo Científico:
    Steffensen type methods for solving nonlinear equations.
Autores:
    Alicia Cordero, José L. Hueso, Eulalia Martínez, Juan R. Torregrosa.
Pag.1-2
Ecuación (1)
Ejemplo : sne_fd_1("(x^2)-(E^x)-(3*x)+2",0.7,0.0000000001, 1)
'''
# CONSISTE EN
#   un método con convergencia cuadrática, libre de derivadas
#   y que utiliza dos evaluaciones por paso.
# ENTRADAS
#   f : función
#   x0 : valor inicial
#   tol : tolerancia
#   graf = parámetro para mostrar la gráfica
# SALIDAS
#   xaprox : aproximacián de x
#   iter : cantidad de iteraciones
#   graf : gráfica resultante
def sne_fd_1(f, x0, tol, graf=1):
    try:  # Se utiliza el Try para detecta un error de Sintaxis.
        x = Symbol('x')  # Se declara x como una variable independiente.
        ff = sympify(f)  # Se pasa de String a Ecuacion.
        itr = 0  # Se define una variable para contar las iteraciones.
        iterations = []  # Se define una lista para graficar las iteraciones.
        errors = []  # Se define una lista para graficar los errores.
        # Si f(Xk+1) es menor a la tolerancia, termina el metodo.
        while (abs(float(ff.subs(x, x0))) >= tol):
            fx = float(ff.subs(x, x0))  # Se calcula el valor de F(x)
            # Se calcula el Xk+1 dado por Steffensen's.
            x1 = x0 - (((fx ** 2) / float(ff.subs(x, x0 + fx) - fx)))
            # Se calcula f(Xk+1) para determinar el error.
            t = abs(float(ff.subs(x, x1)))
            x0 = x1  # Se actualiza el valor de Xk+1
            itr = itr + 1  # Se suma una iteracion
            # Se agregan las iteraciones a una lista para la grafica
            iterations.append(itr)
            # Se agregan los errores a una lista para la grafica
            errors.append(t)
        if graf == 1:
            plt.plot(iterations, errors)  # Llamada para generar la grafica.
            plt.title('Metodo de Steffensen')  # Titulo de la grafica
            plt.xlabel('iteraciones (k)')  # Nombre del eje x
            plt.ylabel('|f(xk)|')  # Nombre del eje y
            plt.grid(True)  # Despliege del grid
            plt.show()  # Despliege del grafico
        elif graf == 0:
            print("Sin grafica")  # Error
        else:
            print("Error: Entrada invalida (graf) : 1 | 0")  # Error
        print("xaprox: " + str(float(x0)) + "\niter: " + str(itr))  # Salida
    except:
        print('Error: Expresión no válida')  # Error: Syntax Error.


'''
Método Improved Ostrowski’s method free from derivatives (IODF)
Artículo: Steffensen type methods for solving nonlinear equations.
Autores: Alicia Cordero, José L. Hueso, Eulalia Martínez, Juan R. Torregrosa.
Página 3
Ecuación (5)
Ejemplo : sne_fd_2('(x^2)-(E^x)-(3*x)+2', 0.7, 0.000001,1)
'''
# CONSISTE EN
#   n método con convergencia 6, libre de derivadas.
# ENTRADAS
#   f : función
#   x0 : valor inicial
#   tol : tolerancia
#   graf = parámetro para mostrar la gráfica
# SALIDAS
#   xaprox : aproximacián de x
#   iter : cantidad de iteraciones
#   graf : gráfica resultante
def sne_fd_2(f, x0, tol, graf=1):
    try:  # Se utiliza el Try para detecta un error de Sintaxis.
        x = Symbol('x')  # Se declara x como una variable independiente.
        ff = sympify(f)  # Se pasa de String a Ecuacion.
        itr = 0  # Se define una variable para contar las iteraciones.
        iterations = []  # Se define una lista para graficar las iteraciones.
        errors = []  # Se define una lista para graficar los errores.
        # Si f(Xk+1) es menor a la tolerancia, termina el metodo
        while (abs(float(ff.subs(x, x0))) >= tol):
            fx = float(ff.subs(x, x0))  # Variable auxiliar para calcular F(x)
            # Variable auxiliar para simplificar el calculo de Yn
            aux1 = float(ff.subs(x, x0 + fx))
            # Variable auxiliar para simplificar el calculo de Yn
            aux2 = float(ff.subs(x, x0 - fx))
            # Se calcula el valor de Yn, requerido por el metodo.
            y = x0 - ((2 * (fx ** 2)) / (aux1 - aux2))
            fy = float(ff.subs(x, y))  # Variable auxiliar para calcular F(y)
            # Se calcula el valor de Zn, requerido por el metodo.
            z = y - (fy * ((y - x0) / ((2 * fy) - fx)))
            fz = float(ff.subs(x, z))  # Variable auxiliar para calcular F(z)
            # Se calcula el valor de Xn+1, requerido por el metodo.
            x1 = z - (fz * ((y - x0) / ((2 * fy) - fx)))
            # Se calcula f(Xk+1) para determinar el error.
            t = abs(float(ff.subs(x, x1)))
            x0 = x1  # Se actualiza el valor de Xk+1
            itr = itr + 1  # Se suma una iteracion para el resultado final.
            # Se agregan las iteraciones a una lista para la grafica
            iterations.append(itr)
            # Se agregan los errores a una lista para la grafica
            errors.append(t)
        if graf == 1:
            plt.plot(iterations, errors)  # Llamada para generar la grafica.
            # Titulo de la grafica
            plt.title('Metodo Improved Ostrowski’s method derivative free')
            plt.xlabel('iteraciones (k)')  # Nombre del eje x
            plt.ylabel('|f(xk)|')  # Nombre del eje y
            plt.grid(True)  # Despliege del grid
            plt.show()  # Despliege del grafico
        elif graf == 0:
            print("Sin grafica")  # Error
        else:
            print("Error: Entrada invalida (graf) : 1 | 0")  # Error
        print("xaprox: " + str(float(x0)) + "\niter: " + str(itr))  # Salida
    except:
        print('Error: Expresión no válida')  # Error: Syntax Error.


'''
Método Op4
Artículo: A class of Steffensen type methods with optimal order of convergence.
Autores:Alicia Cordero, Juan R. Torregrosa.
Páginas 3 y 4.
Ecuación (6)
Ejemplo : sne_fd_3('(x^2)-(E^x)-(3*x)+2', 0.7, 0.000001, 1)
'''
# CONSISTE EN
#   un método con convergencia 4, libre de derivadas.
# ENTRADAS
#   f : función
#   x0 : valor inicial
#   tol : tolerancia
#   graf = parámetro para mostrar la gráfica
# INTERMEDIOS
#   a,b,c,d pertenecen a los numeros Reales.
#   Se usan valores dados en el artículo.
# SALIDAS
#   xaprox : aproximacián de x
#   iter : cantidad de iteraciones
#   graf : gráfica resultante
def sne_fd_3(f, x0, tol, graf=1):
    try:  # Se utiliza el Try para detecta un error de Sintaxis.
        x = Symbol('x')  # Se declara x como una variable independiente.
        ff = sympify(f)  # Se pasa de String a Ecuacion.
        itr = 0  # Se define una variable para contar las iteraciones.
        # Parametros requeridos por el metodo op4 a,b,c tomados de los
        # ejemplos a=b=c=1.
        a = b = c = 1
        # Parametro requerido por el metodo op4 d tomados de los ejemplos d=0.
        d = 0
        iterations = []  # Se define una lista para graficar las iteraciones.
        errors = []  # Se define una lista para graficar los errores.
        # Si f(Xk+1) es menor a la tolerancia, termina el metodo
        while (abs(float(ff.subs(x, x0))) >= tol):
            # Se calcula el valor de Zk+1, requerido por el metodo.
            z = x0 + float(ff.subs(x, x0))
            w = float(ff.subs(x, x0))  # Variable auxiliar para calcular f(x)
            # Se calcula el valor de Yk+1, requerido por el metodo.
            y = x0 - ((w ** 2) / (float(ff.subs(x, z)) - w))
            # Variable auxiliar para calcular la primera parte del
            # denominador de Xk+1
            w1 = (a * float(ff.subs(x,
                                    y))) - (b * float(ff.subs(x,
                                                              z))) / (y - z)
            # Variable auxiliar para calcular la segunda parte del
            # denominador de Xk+1
            w2 = (c * float(ff.subs(x,
                                    y))) - (d * float(ff.subs(x,
                                                              x0))) / (y - x)
            # Se calcula el valor de Xk+1, requerido por el metodo.
            x1 = y - (float(ff.subs(x, y)) / (w1 + w2))
            # Se calcula f(Xk+1) para determinar el error.
            t = abs(float(ff.subs(x, x1)))
            x0 = x1  # Se actualiza el valor de Xk+1
            itr = itr + 1  # Se suma una iteracion para el resultado final.
            # Se agregan las iteraciones a una lista para la grafica
            iterations.append(itr)
            # Se agregan los errores a una lista para la grafica
            errors.append(t)
        if graf == 1:
            plt.plot(iterations, errors)  # Llamada para generar la grafica.
            plt.title('Metodode Steffensen')  # Titulo de la grafica
            plt.xlabel('iteraciones (k)')  # Nombre del eje x
            plt.ylabel('|f(xk)|')  # Nombre del eje y
            plt.grid(True)  # Despliege del grid
            plt.show()  # Despliege del grafico
        elif graf == 0:
            print("Sin grafica")  # Error
        else:
            print("Error: Entrada invalida (graf) : 1 | 0")  # Error
        print("xaprox: " + str(float(x1)) + "\niter: " + str(itr))  # Salida
    except:
        print('Error: Expresión no válida')  # Error: Syntax Error.

'''
Metodo del esquema de Steffensen
Applied Mathematics and Computation
A family of Kurchatov-type methods and its stability
Pag.1 Steffensen’s scheme [25]
Ecuación (2)
Ejemplo : sne_fd_4("x**2 - 18",3,0.0001,1)
'''
# CONSISTE EN
#   el método clásico de Newton libre de derivadas.
# ENTRADAS
#   f : función
#   x0 : valor inicial
#   tol : tolerancia
#   graf = parámetro para mostrar la gráfica
# SALIDAS
#   xaprox : aproximacián de x
#   iter : cantidad de iteraciones
#   graf : gráfica resultante
def sne_fd_4(f, x0, tol, graf=1):
    x = Symbol("x")  # Declaración de x como variable independiente
    try:
        ff = sympify(f)  # Función String -> Ecuación
        fd = diff(ff, x)  # Derivada de la función
        i = 0  # Cantidad de iteraciones
        eje_x = []  # Valores en el eje x para la grafica
        eje_y = []  # Valores en el eje y para la grafica
        xk = x0  # Sucesión de x (valor inicial)
        # Evaluacación de la condición de parada
        while (abs(ff.subs(x, xk)) >= tol):
            eje_x += [i]  # Nuevo punto al eje x
            eje_y += [abs(ff.subs(x, xk))]  # Nuevo punto al eje y
            wk = xk + ff.subs(x, xk)  # Sucesión de w
            # Sucesión de f(xk, wk)
            fxkwk = (ff.subs(x, xk) - ff.subs(x, wk)) / (xk - wk)
            xk = xk - ff.subs(x, xk) / fxkwk  # Sucesión de x
            xk = float(xk)  # xk decimal
            i += 1  # Incremento de iteraciones
        if graf == 1:
            plt.plot(eje_x, eje_y)  # Llamado para gerenar la grafica
            # Titulo de la grafica
            plt.title('Metodo del Esquema de Steffensen')
            plt.xlabel('iteraciones (k)')  # Nombre del eje x
            plt.ylabel('|f(xk)|')  # Nombre del eje y
            plt.grid(True)  # Despliege del grid
            plt.show()  # Despliege del grafico
        elif graf == 0:
            print("Sin grafica")  # Error
        else:
            print("Error: Entrada invalida (graf) : 1 | 0")  # Error
        print("xaprox: " + str(float(xk)) + "\niter: " + str(i))  # Salida
    except:
        print('Error: Expresión no valida')  # Error

'''
Metodo Parida and Gupta
Applied Mathematics Letters
A derivative free iterative method for finding multiple roots of
nonlinear equations
Pag.3 Parida and Gupta [11]
Ecuación (2)
Ejemplo : sne_fd_5("(1-x)**2",3,0.0001,1)
'''
# CONSISTE EN
#   una variacion al metodo de steffensen con convergencia cuadratica
#   funcional en casos especificos.
# ENTRADAS
#   f : función
#   x0 : valor inicial
#   tol : tolerancia
#   graf = parámetro para mostrar la gráfica
# SALIDAS
#   xaprox : aproximacián de x
#   iter : cantidad de iteraciones
#   graf : gráfica resultante
def sne_fd_5(f, x0, tol, graf=1):
    x = Symbol("x")  # Declaración de x como variable independiente
    try:
        ff = sympify(f)  # Función String -> Ecuación
    except:
        return print('Error: Expresión no valida')  # Error
    i = 0  # Cantidad de iteraciones
    eje_x = []  # Valores en el eje x para la grafica
    eje_y = []  # Valores en el eje y para la grafica
    xk = x0  # Sucesión de x (valor inicial)
    # Evaluacación de la condición de parada
    while (abs(ff.subs(x, xk)) >= tol):
        eje_x += [i]  # Nuevo punto al eje x
        eje_y += [abs(ff.subs(x, xk))]  # Nuevo punto al eje y
        # Calculo de la variable t para evitar 0 en el denominador
        t = ff.subs(x, (xk + ff.subs(x, xk))) - ff.subs(x, xk)
        if t < 0:
            t = -1  # Asignacion valor de t apropiado
            # Calculo sucesion x
            xk = xk - ((ff.subs(x, xk))**2 / (t*(ff.subs(x, xk))**2 +
                       ff.subs(x, (xk + ff.subs(x, xk))) - ff.subs(x, xk)))
        else:
            t = 1  # Asignacion valor de t apropiado
            # Calculo sucesion x
            xk = xk - ((ff.subs(x, xk))**2 / (t*(ff.subs(x, xk))**2 +
                       ff.subs(x, (xk + ff.subs(x, xk))) - ff.subs(x, xk)))
        i += 1  # Incremento de iteraciones
        # print(i, "      ", float(xk))
    if graf == 1:
        plt.plot(eje_x, eje_y)  # Llamado para gerenar la grafica
        plt.title('Metodo Parida and Gupta')  # Titulo de la grafica
        plt.xlabel('iteraciones (k)')  # Nombre del eje x
        plt.ylabel('|f(xk)|')  # Nombre del eje y
        plt.grid(True)  # Despliege del grid
        plt.show()  # Despliege del grafico
    elif graf == 0:
        print("Sin grafica")  # Error
    else:
        print("Error: Entrada invalida (graf) : 1 | 0")  # Error
    print("xaprox: " + str(float(xk)) + "\niter: " + str(i))  # Salida

'''
Metodo de Chun’s scheme (CM2).
Solving nonlinear problems by Ostrowski–Chun type parametric families
Pag.7
Ecuación (2)
Ejemplo : sne_fd_6("x**2 - 18",3,0.0001,1)
'''
# CONSISTE EN
#   el método de una familia paramétrica del esquemas originales de Chun sin
#   derivada, con a1 = 1 y b2 = 2
# ENTRADAS
#   f : función
#   x0 : valor inicial
#   tol : tolerancia
#   graf = parámetro para mostrar la gráfica
# SALIDAS
#   xaprox : aproximacián de x
#   iter : cantidad de iteraciones
#   graf : gráfica resultante
def sne_fd_6(f, x0, tol, graf=1):
    x = Symbol("x")  # Declaración de x como variable independiente
    try:
        ff = sympify(f)  # Función String -> Ecuación
        fd = diff(ff, x)  # Derivada de la función
        i = 0  # Cantidad de iteraciones
        eje_x = []  # Valores en el eje x para la grafica
        eje_y = []  # Valores en el eje y para la grafica
        xk = x0  # Sucesión de x (valor inicial)
        # Evaluacación de la condición de parada
        while (abs(ff.subs(x, xk)) >= tol):
            eje_x += [i]  # Nuevo punto al eje x
            eje_y += [abs(ff.subs(x, xk))]  # Nuevo punto al eje y
            zk = xk + ff.subs(x, xk)**2  # Sucesión de z
            # Sucesión de f(zk, xk)
            fzkxk = (ff.subs(x, zk) - ff.subs(x, xk)) / (zk - xk)
            yk = xk - ff.subs(x, xk) / fzkxk  # Sucesión de y
            # Sucesión de x
            xk = yk - (((ff.subs(x, xk) + 2*ff.subs(x, xk)) /
                       (ff.subs(x, xk))) * (ff.subs(x, yk) / fzkxk))
            xk = float(xk)  # xk decimal
            i += 1  # Incremento de iteraciones
        if graf == 1:
            plt.plot(eje_x, eje_y)  # Llamado para gerenar la grafica
            plt.title('Metodo de Chun (CM2)')  # Titulo de la grafica
            plt.xlabel('iteraciones (k)')  # Nombre del eje x
            plt.ylabel('|f(xk)|')  # Nombre del eje y
            plt.grid(True)  # Despliege del grid
            plt.show()  # Despliege del grafico
        elif graf == 0:
            print("Sin grafica")  # Error
        else:
            print("Error: Entrada invalida (graf) : 1 | 0")  # Error
        print("xaprox: " + str(float(xk)) + "\niter: " + str(i))  # Salida
    except:
        print('Error: Expresión no valida')  # Error
