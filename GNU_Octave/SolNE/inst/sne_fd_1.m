%{
Metodo Steffensen (SM)
Tomado del Artículo Científico: Steffensen type methods for solving nonlinear equations.
Autores: Alicia Cordero, José L. Hueso, Eulalia Martínez, Juan R. Torregrosa.
Pag.1-2
Ecuación (1)
Ejemplo : [xaprox, iter] = sne_fd_1("(x^2)-(E^x)-(3*x)+2",0.7,0.0000000001, 1)


CONSISTE EN
  un método con convergencia cuadrática, libre de derivadas y que utiliza dos evaluaciones por paso.
ENTRADAS
  f : función
  x0 : valor inicial
  tol : tolerancia
  graf = parámetro para mostrar la gráfica
SALIDAS
  xaprox : aproximacián de x
  iter : cantidad de iteraciones
  graf : gráfica resultante %}
%}
function [x0, itr] = sne_fd_1 (f, x0, tol, graf = 1)
  pkg load symbolic;
  try # Se utiliza try para atrapar errores de sintaxis
    syms x; # Se define x como variable global.
    fu = sym(f); # Variable auxiliar para la funcion F(X)
    ff = matlabFunction(fu); # Se pasa F(x) a simbolico.
    itr = 0; # Se define una variable para contar las iteraciones.
    iterations = []; # Se define una lista para graficar las iteraciones.
    errors = []; # Se define una lista para graficar los errores.

    while (abs(ff(x0)) >= tol) # Se itera hasta que se cumpla una condicion de parada.
        x1 = x0 - (((ff(x0)) ** 2) / (ff(x0 + ff(x0)) - ff(x0))); #Se calcula el Xk+1 dado por Steffensen's.
        t = abs(ff(x1)); # Se calcula f(Xk+1) para determinar el error.
        x0 = x1; # Se actualiza el valor de Xk+1
        itr = itr + 1; # Se suma una iteracion        
        iterations = [iterations, itr]; # Se llena la lista con las iteraciones.
        errors = [errors, t]; # Se llena la lista con los errores de cada iteracion
     endwhile
     
    if (graf==1)  #Se grafica si el valor de graf es 1.
       plot (iterations,errors); # Llamada para generar la grafica
       title("Metodo de Steffensen", "fontsize", 20); # Titulo del grafico.
       xlabel("iteraciones (k)", "fontsize", 20); # Etiqueta para eje X
       ylabel("|f(xk)|", "fontsize", 20); # Etiqueta para eje Y.
    else 
       if (graf != 0)
         printf('Error: Entrada invalida (graf) : 1 | 0\n'); #Error
       endif
    endif
    
   catch
        printf('Error: Expresión no valida\n'); #Error
   end_try_catch
endfunction

