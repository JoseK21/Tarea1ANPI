%{
Metodo Improved Ostrowski’s method free from derivatives (IODF)
Tomado del Artículo Científico: Steffensen type methods for solving nonlinear equations.
Autores: Alicia Cordero, José L. Hueso, Eulalia Martínez, Juan R. Torregrosa.
Pag.3
Ecuación (5)
Ejemplo : [xaprox, iter] = sne_fd_2("(x^2)-(E^x)-(3*x)+2",0.7,0.0000000001, 1)

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
   graf : gráfica resultante 
%}
function [x0, itr] = sne_fd_2 (f, x0, tol, graf = 1)
  pkg load symbolic;
  try # Se utiliza try para atrapar errores de sintaxis
    syms x; # Se define x como variable global.
    fu = sym(f); # Variable auxiliar para la funcion F(X)
    ff = matlabFunction(fu); # Se pasa F(x) a simbolico.
    itr = 0; # Se define una variable para contar las iteraciones.
    iterations = []; # Se define una lista para graficar las iteraciones.
    errors = []; # Se define una lista para graficar los errores.
    
    while (abs(ff(x0)) >= tol) # Se itera hasta que se cumpla una condicion de parada.
      
      fx = ff(x0); # Variable auxiliar para calcular F(x)
      aux1 = ff(x0 + fx); # Variable auxiliar para simplificar el calculo de Yn
      aux2 = ff(x0 - fx); # Variable auxiliar para simplificar el calculo de Yn
      y = x0 - ((2 * (fx ** 2)) / (aux1 - aux2)); # Se calcula el valor de Yn, requerido por el metodo.
      fy = ff(y); # Variable auxiliar para calcular F(y)
      z = y - (fy * ((y - x0) / ((2 * fy) - fx))); # Se calcula el valor de Zn, requerido por el metodo.
      fz = ff(z); # Variable auxiliar para calcular F(z)

      x1 = z - (fz * ((y - x0) / ((2 * fy) - fx))); # Se calcula el valor de Xn+1, requerido por el metodo.

      t = abs(ff(x1)); # Se calcula f(Xk+1) para determinar el error.
      x0 = x1; # Se actualiza el valor de Xk+1
      itr = itr + 1; # Se suma una iteracion para el resultado final
      iterations = [iterations, itr];  # Se llena la lista con las iteraciones.
      errors = [errors, t]; # Se llena la lista con los errores de cada iteracion
    endwhile
    
    if (graf==1)  #Se grafica si el valor de graf es 1.
       plot (iterations,errors); # Llamada para generar la grafica
       title("Metodo de Ostrowski", "fontsize", 20); # Titulo del grafico.
       xlabel("iteraciones (k)", "fontsize", 20); # Etiqueta para eje X
       ylabel("|f(xk)|", "fontsize", 20); # Etiqueta para eje Y.
    else 
       if (graf != 0)
         printf('Error: Entrada invalida (graf) : 1 | 0\n'); #Error
       endif
    endif
  catch
       printf('Error: Expresión no valida'); # Error.
  end_try_catch
endfunction
 
