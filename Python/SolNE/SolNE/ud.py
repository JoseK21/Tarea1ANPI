import matplotlib.pyplot as plt
from sympy import Symbol, sympify, diff, Subs

'''
Metodo Weerakoon and Fernando
Computers and Structures
Performance of cubic convergent methods for implementing nonlinear
constitutive models
Pag.2  M1 – Weerakoon and Fernando (2000) [21]
Ecuación (5)
Ejemplo : sne_ud_1("x**3 + 6*x**2 - 18",3,0.0001,1)
'''
# CONSISTE EN
#   un método multipunto con convergencia cúbica.
# ENTRADAS
#   f : función
#   x0 : valor inicial
#   tol : tolerancia
#   graf = parámetro para mostrar la gráfica
# SALIDAS
#   xaprox : aproximacián de x
#   iter : cantidad de iteraciones
#   graf : gráfica resultante
def sne_ud_1(f, x0, tol, graf=1):
    x = Symbol("x")  # Declaracion de x como variable independiente
    try:
        ff = sympify(f)  # Funcion String -> Ecuacion
        fd = diff(ff, x)  # Derivada de la funcion
        i = 0  # Cantidad de iteraciones
        eje_x = []  # Valores en el eje x para la grafica
        eje_y = []  # Valores en el eje y para la grafica
        xk = x0  # Sucesion de x (valor inicial)
        # Evaluacacion de la condicion de parada
        while (abs(ff.subs(x, xk)) >= tol):
            eje_x += [i]  # Nuevo punto al eje x
            eje_y += [abs(ff.subs(x, xk))]  # Nuevo punto al eje y
            xk = xk - (2 * ff.subs(x, xk) / (fd.subs(x, xk) + fd.subs(x, xk -
                       (ff.subs(x, xk) / fd.subs(x, xk)))))   # Sucesion de x
            i += 1  # Incremento de iteraciones
        if graf == 1:
            plt.plot(eje_x, eje_y)  # Llamado para gerenar la grafica
            plt.title('Metodo Weerakoon and Fernando')  # Titulo de la grafica
            plt.xlabel('iteraciones (k)')  # Nombre del eje x
            plt.ylabel('|f(xk)|')  # Nombre del eje y
            plt.grid(True)  # Despliege del grid
            plt.show()  # Despliege del grafico5
        elif graf == 0:
            print("Sin grafica")  # Error
        else:
            print("Error: Entrada invalida (graf) : 1 | 0")  # Error
        print("xaprox " + str(float(xk)) + "\niter: " + str(i))  # Salida
    except:
        print('Error: Expresion no valida')  # Error


'''
Metodo Frontini and Sormani, 0zban
Computers and Structures
Performance of cubic convergent methods for implementing nonlinear
constitutive models
Pag.2-3  M2 – Frontini and Sormani (2003) [22], 0zban (2004) [23]
Ecuación (7)
Ejemplo : sne_ud_2("x**3 + 6*x**2 - 18",3,0.0001,1)
'''
# CONSISTE EN
#   un método multipunto redescubrierto propuesto por Traub siguiendo el
#   procedimiento adoptado Kiran et al.
# ENTRADAS
#   f : función
#   x0 : valor inicial
#   tol : tolerancia
#   graf = parámetro para mostrar la gráfica
# SALIDAS
#   xaprox : aproximacián de x
#   iter : cantidad de iteraciones
#   graf : gráfica resultante
def sne_ud_2(f, x0, tol, graf=1):
    x = Symbol("x")  # Declaracion de x como variable independiente
    try:
        ff = sympify(f)  # Funcion String -> Ecuacion
        fd = diff(ff, x)  # Derivada de la funcion
        i = 0  # Cantidad de iteraciones
        eje_x = []  # Valores en el eje x para la grafica
        eje_y = []  # Valores en el eje y para la grafica
        xk = x0  # Sucesion de x (valor inicial)
        # Evaluacacion de la condicion de parada
        while (abs(ff.subs(x, xk)) >= tol):
            eje_x += [i]  # Nuevo punto al eje x
            eje_y += [abs(ff.subs(x, xk))]  # Nuevo punto al eje y
            # Sucesion de x
            xk = xk - (ff.subs(x, xk) /
                       fd.subs(x, xk - 0.5 * ff.subs(x, xk) / fd.subs(x, xk)))
            i += 1  # Incremento de iteraciones
        if graf == 1:
            plt.plot(eje_x, eje_y)  # Llamado para gerenar la grafica
            # Titulo de la grafica
            plt.title('Metodo Frontini and Sormani, 0zban')
            plt.xlabel('iteraciones (k)')  # Nombre del eje x
            plt.ylabel('|f(xk)|')  # Nombre del eje y
            plt.grid(True)  # Despliege del grid
            plt.show()  # Despliege del grafico
        elif graf == 0:
            print("Sin grafica")  # Error
        else:
            print("Error: Entrada invalida (graf) : 1 | 0")  # Error
        print("xaprox " + str(float(xk)) + "\niter: " + str(i))  # Salida
    except:
        print('Error: Expresion no valida')  # Error

'''
Metodo Halley
Servicio de Publicaciones, Universidad de La Rioja
EL METODO DE HALLEY: POSIBLEMENTE, EL METODO MAS REDESCUBIERTO DEL MUNDO
Pag.205
Ecuación (1)
Ejemplo : sne_ud_3("x**3 + 6*x**2 - 18",3,0.0001,1)
'''
# CONSISTE EN
#   un algoritmo de búsqueda de raíz utilizado para funciones de una variable
#   real con una segunda derivada continua.
# ENTRADAS
#   f : función
#   x0 : valor inicial
#   tol : tolerancia
#   graf = parámetro para mostrar la gráfica
# SALIDAS
#   xaprox : aproximacián de x
#   iter : cantidad de iteraciones
#   graf : gráfica resultante
def sne_ud_3(f, x0, tol, graf=1):

    x = Symbol("x")  # Declaración de x como variable independiente
    try:
        func = sympify(f)  # Función String -> Ecuación
    except ValueError:
        return "No se puede parsear la funcion."  # Error
    fp = diff(func, x)  # Derivada de la función
    itera = 0  # Cantidad de iteraciones
    xAx = []  # Valores en el eje x para la grafica
    yAx = []  # Valores en el eje y para la grafica
    # Evaluacación de la condición de parada
    while (abs(func.subs(x, x0)) >= tol):
        fe = func.subs(x, x0)  # funcion evaluada
        fp = diff(func, x)  # funcion prima
        fpe = fp.subs(x, x0)  # funcion prima evaluada
        fpp = diff(fp, x)  # funcion prima prima
        fppe = fpp.subs(x, x0)  # funcion prima prima evaluada

        xAx += [itera]  # Nuevo punto al eje x
        yAx += [abs(func.subs(x, x0))]  # Nuevo punto al eje y
        itera += 1  # Incremento de iteraciones

        if (2*((fpe)**2)-fe*fppe) == 0:  # Revision condicion 0/0
            return "Error de forma 0/0"  # Mensaje de error, caso 0/0

        x0 = x0 - (2*fe*fpe)/(2*((fpe)**2)-fe*fppe)  # Sucesion de x
    if (graf == 1):
        plt.plot(xAx, yAx)  # Llamado para gerenar la grafica
        plt.title("Metodo Halley")  # Titulo de la grafica
        plt.xlabel("iteraciones")  # Nombre del eje x
        plt.ylabel("|f(Xaprox)|")  # Nombre del eje y
        plt.grid(True)  # Despliege del grid
        plt.show()  # Despliege del grafico
    # Salida
    return print("xaprox " + str(float(x0)) + "\niter: " + str(itera))

'''
Metodo Darvishi and Barati
Computers and Structures
Performance of cubic convergent methods for implementing
nonlinear constitutive models
Pag.3 Darvishi and Barati (2007) [25]
Ecuación (11)
Ejemplo : sne_ud_4("x**3 + 6*x**2 - 18",3,0.0001,1)
'''
# CONSISTE EN
#   un método multipunto cúbico para evaluar raíces de un sistema
#   no lineal de ecuaciones.
# ENTRADAS
#   f : función
#   x0 : valor inicial
#   tol : tolerancia
#   graf = parámetro para mostrar la gráfica
# SALIDAS
#   xaprox : aproximacián de x
#   iter : cantidad de iteraciones
#   graf : gráfica resultante
def sne_ud_4(f, x0, tol, graf=1):
    x = Symbol("x")  # Declaracion de x como variable independiente
    try:
        ff = sympify(f)  # Funcion String -> Ecuacion
        fd = diff(ff, x)  # Derivada de la funcion
        i = 0  # Cantidad de iteraciones
        eje_x = []  # Valores en el eje x para la grafica
        eje_y = []  # Valores en el eje y para la grafica
        xk = x0  # Sucesion de x (valor inicial)
        # Evaluacacion de la condicion de parada
        while (abs(ff.subs(x, xk)) >= tol):
            eje_x += [i]  # Nuevo punto al eje x
            eje_y += [abs(ff.subs(x, xk))]  # Nuevo punto al eje y
            xk = xk - (ff.subs(x, xk) / fd.subs(x, xk) -
                       ff.subs(x, xk - ff.subs(x, xk) / fd.subs(x, xk)) /
                       fd.subs(x, xk))   # Sucesion de x
            i += 1  # Incremento de iteraciones
        if graf == 1:
            plt.plot(eje_x, eje_y)  # Llamado para gerenar la grafica
            plt.title('Metodo Darvishi and Barati')  # Titulo de la grafica
            plt.xlabel('iteraciones (k)')  # Nombre del eje x
            plt.ylabel('|f(xk)|')  # Nombre del eje y
            plt.grid(True)  # Despliege del grid
            plt.show()  # Despliege del grafico
        elif graf == 0:
            print("Sin grafica")  # Error
        else:
            print("Error: Entrada invalida (graf) : 1 | 0")  # Error
        print("xaprox " + str(float(xk)) + "\niter: " + str(i))  # Salida
    except:
        print('Error: Expresion no valida')  # Error

'''
Metodo Kou et al
Computers and Structures
Performance of cubic convergent methods for implementing
nonlinear constitutive models
Pag.3 M5 - Kou et al. (2006) [26]
Ecuación (13)
Ejemplo : sne_ud_5("x**3 + 6*x**2 - 18",3,0.0001,1)
'''
# CONSISTE EN
#   un método multipunto cúbico, caso especial del método multipunto
#   propuesto por Ostrowski
# ENTRADAS
#   f : función
#   x0 : valor inicial
#   tol : tolerancia
#   graf = parámetro para mostrar la gráfica
# SALIDAS
#   xaprox : aproximacián de x
#   iter : cantidad de iteraciones
#   graf : gráfica resultante
def sne_ud_5(f, x0, tol, graf=1):
    x = Symbol("x")  # Declaracion de x como variable independiente
    try:
        ff = sympify(f)  # Funcion String -> Ecuacion
        fd = diff(ff, x)  # Derivada de la funcion
        i = 0  # Cantidad de iteraciones
        eje_x = []  # Valores en el eje x para la grafica
        eje_y = []  # Valores en el eje y para la grafica
        xk = x0  # Sucesion de x (valor inicial)
        # Evaluacacion de la condicion de parada
        while (abs(ff.subs(x, xk)) >= tol):
            eje_x += [i]  # Nuevo punto al eje x
            eje_y += [abs(ff.subs(x, xk))]  # Nuevo punto al eje y
            # Sucesion de x
            xk = xk + (ff.subs(x, xk) /
                       fd.subs(x, xk) - ff.subs(x, xk + ff.subs(x, xk) /
                                                fd.subs(x, xk)) /
                       fd.subs(x, xk))
            i += 1  # Incremento de iteraciones
        if graf == 1:
            plt.plot(eje_x, eje_y)  # Llamado para gerenar la grafica
            plt.title('Metodo Kou et al')  # Titulo de la grafica
            plt.xlabel('iteraciones (k)')  # Nombre del eje x
            plt.ylabel('|f(xk)|')  # Nombre del eje y
            plt.grid(True)  # Despliege del grid
            plt.show()  # Despliege del grafico
        elif graf == 0:
            print("Sin grafica")  # Error
        else:
            print("Error: Entrada invalida (graf) : 1 | 0")  # Error
        print("xaprox " + str(float(xk)) + "\niter: " + str(i))  # Salida
    except:
        print('Error: Expresion no valida')  # Error

'''
Metodo Chun and Kim
Computers and Structures
Performance of cubic convergent methods for implementing
nonlinear constitutive models
Pag.3  M6 - Chun and Kim (2010) [27]
Ecuación (15)
Ejemplo : sne_ud_6("x**3 + 6*x**2 - 18",3,0.0001,1)
'''
# CONSISTE EN
#   un método multipunto cúbico utilizando análisis geométricos y utilizando
#   métodos de diferencias finitas para aproximar las segundas derivadas.
# ENTRADAS
#   f : función
#   x0 : valor inicial
#   tol : tolerancia
#   graf = parámetro para mostrar la gráfica
# SALIDAS
#   xaprox : aproximacián de x
#   iter : cantidad de iteraciones
#   graf : gráfica resultante
def sne_ud_6(f, x0, tol, graf=1):
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
            # Sucesión de x
            xk = xk - 0.5 * (3 - fd.subs(x, xk - ff.subs(x, xk) /
                                         fd.subs(x, xk)) /
                             fd.subs(x, xk)) * ff.subs(x, xk) / fd.subs(x, xk)
            i += 1  # Incremento de iteraciones
        if graf == 1:
            plt.plot(eje_x, eje_y)  # Llamado para gerenar la grafica
            plt.title('Metodo Chun & Kim')  # Titulo de la grafica
            plt.xlabel('iteraciones (k)')  # Nombre del eje x
            plt.ylabel('|f(xk)|')  # Nombre del eje y
            plt.grid(True)  # Despliege del grid
            plt.show()  # Despliege del grafico
        elif graf == 0:
            print("Sin grafica")  # Error
        else:
            print("Error: Entrada invalida (graf) : 1 | 0")  # Error
        print("xaprox " + str(float(xk)) + "\niter: " + str(i))  # Salida
    except:
        print('Error: Expresión no valida')  # Error
