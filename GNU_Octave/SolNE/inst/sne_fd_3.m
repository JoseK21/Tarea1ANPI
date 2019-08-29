%{
Metodo Op4
Tomado del Artículo Científico: A class of Steffensen type methods with optimal order of convergence.
Autores: Alicia Cordero, Juan R. Torregrosa.
Pag.3-4
Ecuación (6)
Ejemplo : [xaprox, iter] = sne_fd_3("(x^2)-(E^x)-(3*x)+2",0.7,0.0000000001, 1)

 CONSISTE EN
   un método con convergencia 4, libre de derivadas.
 ENTRADAS
   f : función
   x0 : valor inicial
   tol : tolerancia
   graf = parámetro para mostrar la gráfica
 INTERMEDIOS
   a,b,c,d pertenecen a los numeros Reales. Se usan valores dados en el artículo. 
 SALIDAS
   xaprox : aproximacián de x
   iter : cantidad de iteraciones
   graf : gráfica resultante
%}
function [x0, itr] = sne_fd_3 (f, x0, tol, graf = 1)
  pkg load symbolic;
  try # Se utiliza try para atrapar errores de sintaxis
    syms x; # Se define x como variable global.
    fu = sym(f); # Variable auxiliar para la funcion F(X)
    ff = matlabFunction(fu); # Se pasa F(x) a simbolico.
    itr = 0; # Se define una variable para contar las iteraciones.
    iterations = []; # Se define una lista para graficar las iteraciones.
    errors = []; # Se define una lista para graficar los errores.
    
    a = b = c = 1; # Parametros requeridos por el metodo op4 a,b,c tomados de los ejemplos a=b=c=1.
    d = 0; # Parametro requerido por el metodo op4 d tomados de los ejemplos d=0.
    
    while (abs(ff(x0)) >= tol)  # Se itera hasta que se cumpla una condicion de parada.
      z = x0 + ff(x0); # Se calcula el valor de Zk+1, requerido por el metodo.
      w = ff(x0); # Variable auxiliar para calcular f(x)

      y = x0 - ((w ** 2) / (ff(z) - w)); # Se calcula el valor de Yk+1, requerido por el metodo.
      w1 = (a * ff(y)) - (b * ff(z)) / (y - z); # Variable auxiliar para calcular la primera parte del denominador de Xk+1                                
      w2 = (c * ff(y)) - (d * ff(x0)) / (y - x0); # Variable auxiliar para calcular la segunda parte del denominador de Xk+1 

      x1 = y - (ff(y) / (w1 + w2)); # Se calcula el valor de Xk+1, requerido por el metodo.
      
      t = abs(ff(x1)); # Se calcula f(Xk+1) para determinar el error.
      x0 = x1; # Se actualiza el valor de Xk+1
      itr = itr + 1; # Se suma una iteracion
      iterations = [iterations, itr]; # Se llena la lista con las iteraciones.
      errors = [errors, t]; # Se llena la lista con los errores de cada iteracion
      
    endwhile
   
    if (graf==1)  #Se grafica si el valor de graf es 1.
       plot (iterations,errors); # Llamada para generar la grafica
       title("Metodo Op4", "fontsize", 20); # Titulo del grafico.
       xlabel("iteraciones (k)", "fontsize", 20); # Etiqueta para eje X
       ylabel("|f(xk)|", "fontsize", 20); # Etiqueta para eje Y.
    else 
       if (graf != 0)
         printf('Error: Entrada invalida (graf) : 1 | 0\n'); #Error
       endif
    endif
  catch
       printf('Error: Expresión no valida');  # Error
  end_try_catch
endfunction

