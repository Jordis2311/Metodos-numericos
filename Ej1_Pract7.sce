// Ejercicio 1 de la práctica 7
clc()
clear()

// Polinomio interpolador de Lagrange
function w = Lagrange(p,x,y)
        w = 0
        n = length(x)
    for i=1:n do
        w = w + L(p,i,x)*y(i)
    end
endfunction

// Función L_i(x) del polinomio interpolador de Lagrange
function w = L(p,i,x)
    w = 1
    n = length(x)
    for j=1:n do
        if j<>i then
            w = w*(p-x(j))/(x(i)-x(j))
        end
    end
endfunction

// 
function w = DD_Newton(p,x,y)
    w = 0
    n = length(x)
    for j=n:-1:2
        w = (w + DD(x(1:j),Y(1:j)))*(p-x(j-1))
    end
    w = w + y(1)
endfunction


// Diferencias divididas
function w = DD(x,y)
    n = length(x)
    if n==2 then
        w = (y(n)-y(1))/(x(n)-x(1))
    else
        w = (DD(x(2:n),y(2:n))-DD(x(1:n-1),y(1:n-1)))/(x(n)-x(1))
    end
endfunction
