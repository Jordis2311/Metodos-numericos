// Ejercicio 9 de la Práctica 7
clc()
clear()

// Función error del ejercicio
function w = err(x,n)
    // Entrada: x = valor real; n = grado del polinomio interpolante
    // Salida: w = error de interpolar por diferencias divididas en el punto x
    nodos = -5:10/n:5 // long(nodos) = n+1
    y = 1./(1.+nodos.^2)
    w = 1./(1.+x.^2) - horner(DD_Newton(nodos,y),x)
endfunction

// Método de Diferencias Divididas de Newton
function w = DD_Newton(x,y)
    // Entrada: x,y = vectores puntos de interpolación (x,y)
    // Salida: w = polinomio de diferencias divididas de Newton
    w = 0
    s = poly(0,"x")
    n = length(x)
    for j=n:-1:2
        w = w*(s-x(j-1)) + DD(x(1:j),y(1:j))*(s-x(j-1))
    end
    w = w + y(1)
endfunction

function y = j0(x)
    y = (1/%pi) * integrate('cos(x*sin(t))','t',0,%pi)
endfunction

// Diferencias divididas
function w = DD(x,y)
    // Entrada: 
    // Salida: 
    n = length(x)
    if n==1 then
        w = y(1)
    else
        w = (DD(x(2:n),y(2:n))-DD(x(1:n-1),y(1:n-1)))/(x(n)-x(1))
    end
endfunction

// - Ejercicio - //
x = -5:0.01:5
plot2d(x,[err(x,2)' err(x,4)' err(x,6)' err(x,10)' err(x,14)'],[2,3,4,5,6],leg="n=2@n=4@n=6@n=10@n=14")
a=gca()
a.x_location = "origin"
a.y_location = "origin"

plot(x,err(x,2))
