%{
Metodo de Parida and Gupta
Tomado del Artículo Científico: A derivative free iterative method for finding multiple roots of nonlinear equations
Pag.3
Ecuación (2)
Ejemplo : [xaprox, iter] = sne_fd_5("(1-x)^2",3,0.0001, 1)

 CONSISTE EN
   una variacion al metodo de Steffensen con convergencia cuadratica, funcional en casos especificos.
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
function [x0, itr] = sne_fd_5(f, x0, tol, graf = 1)
  pkg load symbolic;
  try # Se utiliza try para atrapar errores de sintaxis
    syms x; # Se define x como variable global.
    fu = sym(f); # Variable auxiliar para la funcion F(X)    
    ff = matlabFunction(fu); # Se pasa F(x) a simbolico.
    
        
    t_func = strcat('@(x)',f); % @x + función
    f = str2func(t_func); % función string a ecuación
    
    itr = 0; # Se define una variable para contar las iteraciones.
    iterations = []; # Se define una lista para graficar las iteraciones.
    errors = []; # Se define una lista para graficar los errores.

    while (abs(ff(x0)) >= tol)  # Se itera hasta que se cumpla una condicion de parada. 
      b = ff(x0 + ff(x0)) - ff(x0)  # Calculo de la variable t para evitar 0 en el denominador 
      if (b <0)
        x1 = x0 - (ff(x0))**2 / (-1*(ff(x0))**2 + ff(x0 + ff(x0)) - ff(x0))  # Calculo sucesion x  
      else
        x1 = x0 - (ff(x0))**2 / ((ff(x0))**2 + ff(x0 + ff(x0)) - ff(x0))  # Calculo sucesion x 
      endif 
      t = abs(ff(x0)); # Se calcula f(Xk+1) para determinar el error.
      x0 = x1; # Se actualiza el valor de Xk+1
      itr = itr + 1; # Se suma una iteracion
      iterations = [iterations, itr]; # Se llena la lista con las iteraciones.
      errors = [errors, t]; # Se llena la lista con los errores de cada iteracion
      
    endwhile
   
    if (graf==1)  #Se grafica si el valor de graf es 1.
       plot (iterations,errors); # Llamada para generar la grafica
       title("Metodo de Chun and Kim", "fontsize", 20); # Titulo del grafico.
       xlabel("iteraciones (k)", "fontsize", 20); # Etiqueta para eje X
       ylabel("|f(xk)|", "fontsize", 20); # Etiqueta para eje Y.
    else 
       if (graf != 0)
         printf('Error: Entrada invalida (graf) : 1 | 0\n'); #Error
       endif
    endif
  catch err
    warning(err.identifier, err.message); 
    #printf('Error: Expresión no valida');  # Error
  end_try_catch
endfunction




