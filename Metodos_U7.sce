// Método de Eliminación Gaussiana con pivoteo parcial
function [x,a] = gausselimPP(A,b)
[nA,mA] = size(A) 
[nb,mb] = size(b)
a = [A b]; // Matriz aumentada
n = nA;    // Tamaño de la matriz
// Eliminación progresiva con pivoteo parcial
for k=1:n-1
    kpivot = k; amax = abs(a(k,k));  //pivoteo
    for i=k+1:n
        if abs(a(i,k))>amax then
            kpivot = i; amax = a(i,k);
        end;
    end;
    temp = a(kpivot,:); a(kpivot,:) = a(k,:);
    a(k,:) = temp
    
    for i=k+1:n
        for j=k+1:n+1
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k)
        end;
        for j=1:k        // no hace falta para calcular la solución x
            a(i,j) = 0;  // no hace falta para calcular la solución x
        end              // no hace falta para calcular la solución x
    end
end
// Sustitución regresiva
x(n) = a(n,n+1)/a(n,n)
for i = n-1:-1:1
    sumk = 0
    for k=i+1:n
        sumk = sumk + a(i,k)*x(k)
    end
    x(i) = (a(i,n+1)-sumk)/a(i,i)
end
endfunction

// Lagrange(p,x,y) Metodo de interpolacion de lagrange
// p: punto que se quiere aproximar
// x,y: vectores de los Nodos de interpolacion (x,y)
// Ejemplo
// Lagrange(2,[1 3 5],[7 4 1]) = 5.5
function w = Lagrange(p,x,y)
        w = 0
        n = length(x)
    for i=1:n do
        w = w + L(p,i,x)*y(i)
    end
endfunction

// L(p,i,x)
// Función L_i(p) del polinomio interpolador de Lagrange
// p: punto a aproximar
// x: vector de los nodos de interpolacion
// i: indice de la funcion L
// retorna el resultado de la funcion L(p)
// Ejemplo
// L(3,2,[2 4 6]) = 0.75
function w = L(p,i,x)
    w = 1
    n = length(x)
    for j=1:n do
        if j<>i then
            w = w*(p-x(j))/(x(i)-x(j))
        end
    end
endfunction

// DD_Newton(p,x,y) Metodo de Diferencias divididas de Newton
// p: Punto que se quiere aproximar
// x,y: vectores de los nodos de interpolacion
// la imagen del punto p utilizando el polinomio interpolante
// resultante del metodo de diferencias divididas de newton
// en los puntos (x(i),y(i))
// Ejemplo
// c = DD(2.5,[2 3 4],[1 5 7])
// c = 3.25
function w = DD_Newton(p,x,y)
    w = 0
    n = length(x)
    for j=n:-1:2
        w = (w + DD(x(1:j),y(1:j)))*(p-x(j-1))
    end
    w = w + y(1)
endfunction

// Polinomio del Método de Diferencias Divididas de Newton
// Retorna el polinomio interpolante resultante
// de aplicar el método de Diferencias Divididas de Newton
// en los puntos(x(i),y(i))
// Ejemplo
// p = DD_Newton_pol([2 3 4],[2 1 5])
// p = 19 -13.5x +2.5x^2
function w = DD_Newton_pol(x,y)
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

// Diferencias divididas
// x,y: vectores de puntos de interpolación (x,y)
// retorna el resultado de las diferencias divididas
// en los vectores x e y
// Ejemplo:
// DD([2 3 4],[2 1 5]) = 2.5
function w = DD(x,y)
    n = length(x)
    if n==2 then
        w = (y(n)-y(1))/(x(n)-x(1))
    else
        w = (DD(x(2:n),y(2:n))-DD(x(1:n-1),y(1:n-1)))/(x(n)-x(1))
    end
endfunction

// MinCuad_pol(A,b)
// A: Matriz de aproximacion de funciones en base a x
// b: vector de imagenes de x
// Aproximación polinomial de mínimos cuadrados polinomial para matrices con rango completo
// retorna el polinomio de mínimos cuadrados junto con el vector de errores (eps = Ax-b)
// Ejemplo
// A = [1 1 1;1 3 9 ;1 2 4;1 4 16]
// b = [3 4 2 1]
// [p,err] = minCuad_pol(A,b)
// p = 1 +2.1x -0.5x^2
// err = [-0.4 -1.2 1.2 0.4]'
function [p,err] = MinCuad_pol(A,b)
     [w,a] = gausselimPP((A')*A,(A')*(b'))
     p = poly(w,"x","coeff")
     err = A*w-b'
endfunction

// A_mc(x,n)
// x: vector 1xn
// n: grado de la aproximación
// Dado Dado el vector x construye la matriz de Vandermonde que se utilizara
// para el metodo de minimos cuadrados
// Ejemplo
// A = A_mc([1 3 2,4],2)
// A = [1 1 1;1 3 9 ;1 2 4;1 4 16]
function A = A_mc(x,n)
    m = length(x)
    A = ones(m,1)
    for j=2:(n+1) do
        A = [A,(x').^(j-1)]
    end
endfunction

// raices_Chebyshev_intervaloç
// a,b: Intevalo donde se buscan las raices
// n: Cantidad de raices necesitadas
// retorna las raices de un polinomio de chebyshev
// adaptadas para que todas se encuentre en un intervalo [a,b]
// Ejemplo
// r = raices_Chebyshev_intervalo(1,2,4)
// r = [1.0380602 1.9619398 1.3086583 1.6913417]'
function x = raices_Chebyshev_intervalo(a,b,n)
    p = Chebyshev(n)
    t = roots(p)
    if (a<>-1)&&(b<>1) then
        for i = 1:n
            x(i) = ((b+a)+t(i)*(b-a))/2
        end
    else
        x = t
    end    
endfunction

// Chebyshev(n)
// n: grado del polinomio 
// retorna el polinomio monomico de chevyshev
// Ejemplo
// p = Chebyshev(3)
//`p = -3x +4x^3
function p = Chebyshev(n)
    if n == 0 then
        p = 1
    else 
        if n == 1 then
            p = poly([0 1],'x','coeff')
        else
            p2 = poly([0 2],'x','coeff')
            p = p2 * Chebyshev(n-1) - Chebyshev(n-2)
        end
    end
endfunction

// Chebyshev_mon(n)
// n: grado del polinomio 
// retorna el polinomio monomico de chevyshev
// Ejemplo
// p = Chebyshev_mon(3)
//`p = -0.75x + x^3
function p = Chebyshev_mon(n)
    p = Chebyshev(n)* 1/2**(n-1)
endfunction 

// DD_newtonP_chev(f,a,b,n)
// f: funcion
// a,b: Intervalos de interpolacion
// n: Grado del polinomio interpolante
// Dada una funcion, un intervalo y un nro natural n
// se retorna un polinomio de grado n que interpola
// a la funcion f en las n+1 raices del polinomio
// de chevyshev adaptadas a el intervalo dado
// Ejemplo:
// p = DD_newtonP_chev('cos(x)',0,%pi/2,3)
// p = 0.9984416 +0.0319396x -0.6049283x^2  +0.1142627x^3
function p = DD_newtonP_chev(f,a,b,n)
    deff("y = F(x)","y="+f);
    x = real(raices_Chebyshev_intervalo(a,b,n+1)) //Las raices del polinomio de chevishev siempre son reales, lo agrego para que scilab no muestre la parte imaginaria
    y = F(x)
    p = DD_Newton_pol(x,y)
endfunction

// graf_error_int(f,w,a,b)
// f = funcion; w = polinomio interpolante de f
// a,b = extremos de interpolacion
// c = cota de error
// Dados una funcion, su polinomio interpolante
// y los extremos de interpolacion se grafica el error de interpolacion
function graf_error_int(f,w,a,b,c)
    deff("y = F(x)","y="+f);
    z = (a-1):0.01:(b+1)
    deff("v = G(x)","v = F(x) - horner(w,z)")
    plot2d(z,G(z),rect = [(a-1), -c, (b+1), c],leg="f(x) - P(x)")
    a=gca()
    a.x_location = "origin"
    a.y_location = "origin"
endfunction
