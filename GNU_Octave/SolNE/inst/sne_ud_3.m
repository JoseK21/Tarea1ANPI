%{
Metodo de Halley
Tomado del Artículo Científico: Servicio de Publicaciones, Universidad de La Rioja
EL METODO DE HALLEY: POSIBLEMENTE, EL METODO MAS REDESCUBIERTO DEL MUNDO
Pag.205  
Ecuación (1)
Ejemplo : [xaprox, iter] = sne_ud_3("x**3 + 6*x**2 - 18",3,0.0001,1)

 CONSISTE EN
   un algoritmo de búsqueda de raíz utilizado para funciones de una variable real con una segunda derivada continua.
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
function [x0, itr] = sne_ud_3 (f, x0, tol, graf = 1)
  pkg load symbolic;
  try # Se utiliza try para atrapar errores de sintaxis
    syms x; # Se define x como variable global.
    fu = sym(f); # Variable auxiliar para la funcion F(X)    
    ff = matlabFunction(fu); # Se pasa F(x) a simbolico.
    
        
    t_func = strcat('@(x)',f); % @x + función
    f = str2func(t_func); % función string a ecuación
    fd = matlabFunction(diff(f, x)); % derivada de la función
    
    fdd = matlabFunction(diff(diff(f, x))); % segunda derivada de la función
    itr = 0; # Se define una variable para contar las iteraciones.
    iterations = []; # Se define una lista para graficar las iteraciones.
    errors = []; # Se define una lista para graficar los errores.

    while (abs(ff(x0)) >= tol)  # Se itera hasta que se cumpla una condicion de parada.   
      x1 = x0 - 2*ff(x0)*fd(x0) / ( 2*(fd(x0)**2 - ff(x0)*fdd(x0) ));  % sucesión de x      ;  # Sucesion de x
      
      t = abs(ff(x0)); # Se calcula f(Xk+1) para determinar el error.
      x0 = x1; # Se actualiza el valor de Xk+1
      itr = itr + 1; # Se suma una iteracion
      iterations = [iterations, itr]; # Se llena la lista con las iteraciones.
      errors = [errors, t]; # Se llena la lista con los errores de cada iteracion
      
    endwhile
   
    if (graf==1)  #Se grafica si el valor de graf es 1.
       plot (iterations,errors); # Llamada para generar la grafica
       title("Metodo de Halley", "fontsize", 20); # Titulo del grafico.
       xlabel("iteraciones (k)", "fontsize", 20); # Etiqueta para eje X
       ylabel("|f(xk)|", "fontsize", 20); # Etiqueta para eje Y.
    else 
       if (graf != 0)
         disp('Error: Entrada invalida (graf) : 1 | 0\n'); #Error
       endif
    endif
  catch err
    warning(err.identifier, err.message); 
    disp('Error: Expresión no valida');  # Error
  end_try_catch
endfunction




