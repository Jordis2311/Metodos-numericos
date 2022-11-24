function x = sumar3(a,b,c)
    x = a + b + c;
endfunction

function x = fact(n)
    f = 1;
    if n == 0 then x = 1;
    else x = prod(1:1:n);
    end
endfunction

function [a,b] = sumaresta(x,y)
    a = x + y;
    b = x - y;
endfunction

function valor = derivada(f,v,h)
    valor = (f(v+h) - f(v)) / h
endfunction

function x = dev_m(f,m,v,h)
    if m == 0 then
       x = f(v);
    else
       x = (dev_m(f,m-1,v+h,h) - dev_m(f,m-1,v,h)) / h
    end
endfunction

function r  = raices(p)
    c = coeff(p, 0);
    b = coeff(p, 1);
    a = coeff(p, 2);
    dis = (b**2 - 4*a*c);
    if(b < 0) then;
        r(1) = (-b + sqrt(dis))/(2*a);
        r(2) = (2*c)/(-b + sqrt(dis));
    else
        r(1) = (2*c)/(-b - sqrt(dis));
        r(2) = (-b - sqrt(dis))/(2*a); 
    end 
endfunction

function valor = horner2(pol,v)
    
    grado = degree(pol);
    
    b = coeff(p,grado);
    grado = grado - 1;
    q = 0;
        
    while grado > -1,
        q = q + b*(v**grado);
        b = coeff(p,grado) + v*b;
        grado = grado - 1;

    end
    
    valor(1) = b;
    valor(2) = q;
    
endfunction

function valor = tay(fun, v, n, a, h)
    
    valor = 0;
    
    for i = 0:n
        valor = valor + ( (1 / factorial(i)) * dev_m(fun,i,a,h) *(v - a)**i);
    end
    
endfunction

function salida = newton(fun,x0,tol,iter)
    deff("y = f(x)","y="+fun);
    i = 0;
    x1 = x0 - (f(x0) / numderivative(f,x0));
    while (abs(x1-x0) > tol) && (i<iter)
        i = i+1;
        x0 = x1;
        x1 = x0 - (f(x0) / numderivative(f,x0));
    end
    
    if (abs(x1-x0)<tol) then disp('Se alcanzo la tolerancia')
    end
    disp('Iteraciones realizadas: ');
    disp(i);
    salida = x1;
    
endfunction 

function salida = biseccion(fun,a0,b0,tol)
    deff("y = f(x)","y="+fun);
    i = 0;
    a = a0;
    b = b0;
    c = (a + b) / 2;
    
    if f(a)*f(b) >= 0 then
        disp("Intervalo no valido maestro");
        return
    end
    
    while (abs(a - b) > tol)
        i = i + 1;
        if (f(c)*f(a)) < 0 then
            b = c;
            c = (a + b) / 2; 
        end
        if (f(c)*f(b)) < 0 then
            a = c;
            c = (a + b) / 2;
        end
    end
    
    if (abs(f(c))<tol) then disp('Se alcanzo la tolerancia')
    end
    disp('Iteraciones realizadas: ');
    disp(i);
    salida = c;  
 
endfunction

function salida = secante(fun,x0,y0,tol,iter)
    deff("y = f(x)","y="+fun);
    i = 0;
    x = x0;
    y = y0;
    c = x - (f(x) * (x - y) / (f(x) - f(y)));
    while (abs(x-y) > tol) && (i < iter)
        i = i + 1;
        y = x;
        x = c;
        c = x - (f(x) * (x - y) / (f(x) - f(y)));
    end
    if (abs(x-y)<tol) then disp('Se alcanzo la tolerancia')
    end
    disp('Iteraciones realizadas: ');
    disp(i);
    salida = c;
endfunction

function x = resolver_MTS(m,b,n)
    x(n) = b(n) / m(n,n);
    
    for i = 1:(n-1)
        r = 0;
        for j = 0:i
           r = r + x(n-j) * m(n-i,n-j);
        end
        x(n-i) = (b(n-i) - r) / m(n-i,n-i);
    end
     
endfunction

function x = resolver_MTI(m,b,n)
    x(1) = b(1) / m(1,1);
    
    for i = 2:n
        r = 0;
        for j = 1:(i-1)
            r = r + x(j) * m(i,j);
        end
        x(i) = (b(i) - r) / m(i,i);
    end
     
endfunction

function x = gaussiana_boluda(m,b,n)
    [L,U,P] = lu(m);
    v = P * b;
    y = resolver_MTI(L,v,n);
    x = resolver_MTS(U,y,n);
    
endfunction

function [x,a] = gausselim(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana sin pivoteo.  

[nA,mA] = size(A) 
[nb,mb] = size(b)

pasos = 0;

if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b');
    abort;
end;

a = [A b]; // Matriz aumentada

// Eliminación progresiva
n = nA;
for k=1:n-1
    for i=k+1:n
        for j=k+1:n+1
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
            pasos = pasos + 1;
        end;
        for j=1:k        // no hace falta para calcular la solución x
            a(i,j) = 0;  // no hace falta para calcular la solución x
        end              // no hace falta para calcular la solución x
    end;
end;

// Sustitución regresiva
x(n) = a(n,n+1)/a(n,n);
for i = n-1:-1:1
    sumk = 0
    for k=i+1:n
        sumk = sumk + a(i,k)*x(k);
        pasos = pasos + 1;
    end;
    x(i) = (a(i,n+1)-sumk)/a(i,i);
    pasos = pasos + 1;
end;
disp('La cantidad de pasos es ');
disp(pasos);
endfunction

// GaussElim modificada para resolver varios sistemas (ej 3)
function [x,a] = gausselim_mod(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana sin pivoteo.  

[nA,mA] = size(A) 
[nb,mb] = size(b)

if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b');
    abort;
end;


a = [A b]; // Matriz aumentada

// Eliminación progresiva
n = nA;
for k=1:n-1
    for i=k+1:n
        for j=k+1:n+mb
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
        end;
        for j=1:k        // no hace falta para calcular la solución x
            a(i,j) = 0;  // no hace falta para calcular la solución x
        end              // no hace falta para calcular la solución x
    end;
end;

// Sustitución regresiva
for j = 1:n
    x(n,j) = a(n,n+j)/a(n,n);
    for i = n-1:-1:1
        sumk = 0
        for k=i+1:n
            sumk = sumk + a(i,k)*x(k,j);
        end;
        x(i,j) = (a(i,n+j)-sumk)/a(i,i);
    end
end


// Si decimos que b es la identidad, x va a ser la inversa de A
endfunction

function x = resolver_TRI(A,b)
    [nA,mA] = size(A);
    [nb,mb] = size(b);
    if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end 
    //Triangularizar
    a = [A b];
    for i = 2:nA-1
        a(i,i) = a(i,i) - a(i-1,i) * a(i,i-1) / a(i-1,i-1);
        a(i,nA+1) = a(i,nA+1) - a(i-1,nA+1) * a(i,i-1) / a(i-1,i-1);
        a(i,i-1) = 0
        
    end
    
    a(nA,nA) = a(nA,nA) - a(nA-1,nA) * a(nA,nA-1) / a(nA - 1,nA-1);
    a(nA,nA+1) = a(nA,nA+1) - a(nA-1,nA+1) * a(nA,nA-1) / a(nA - 1,nA-1);
    a(nA,nA-1) = 0;
    
    // Sustitución regresiva
    n = nA;
    x(n) = a(n,n+1)/a(n,n);
    for i = n-1:-1:1
        sumk = 0
        for k=i+1:n
            sumk = sumk + a(i,k)*x(k);
        end;
        x(i) = (a(i,n+1)-sumk)/a(i,i);
    end;
    
endfunction

function B = intercambiar_filas(A,i,j)
    [n,m] = size(A);
    aux = A(i,:);
    A(i,:) = A(j,:);
    A(j,:) = aux;
    B = A;
endfunction

function [L,U,P] = factorizar_PLU(A)
    [nA,mA] = size(A);
    if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    end
    n = nA;
    P = eye(n,n);
    L = eye(n,n);
    U = A;
    
    for k = 1:n-1
        i = k;
        for j = k+1:n
            //Seleccionar i ≥ k que maximiza |Uik|
            if abs(U(i,k)) < abs(U(j,k)) then
                i = j;
            end
        end
        if i <> k then
           //intercambio de filas
           for j = k:n
              aux = U(k,j);
              U(k,j) = U(i,j);
              U(i,j) = aux;
           end
           for j = 1:k-1
              aux = L(k,j);
              L(k,j) = L(i,j);
              L(i,j) = aux;
           end
           P = intercambiar_filas(P,i,k);
        end
        for j = k+1:n
            L(j,k) = U(j,k)/U(k,k);
            U(j,k:n) = U(j,k:n) - L(j,k)*U(k,k:n);
        end
    end    
    
endfunction

function [L,U] = doolit(A)
    [nA,mA] = size(A);
    if nA<>mA then
        error('La matriz A debe ser cuadrada');
        abort;
    end
    n = nA;
    L = zeros(n,n);
    U = zeros(n,n);
    
    
    for i = 1:n
        for k = i:n
            suma = 0;
            for j = 1:i
                suma = suma + (L(i,j) * U(j,k));
            end
            U(i,k) = A(i,k) - suma;
        end
        
        for k = i:n
            if (i == k) then
                L(i,i) = 1;
            else
                suma = 0;
                for j = 1:i
                    suma = suma + L(k,j) * U(j,i);
                end
                L(k,i) = (A(k,i) - suma) / U(i,i);
            end    
        end
    end
endfunction

function x = resolver_con_DL(A,b)
    [nA,mA] = size(A);
    if nA<>mA then
        error('La matriz A debe ser cuadrada');
        abort;
    end
    n = nA;
    [L,U] = doolit(A);
    y = resolver_MTI(L,b,n);
    x = resolver_MTS(U,y,n);
endfunction

function [U, ind] = Cholesky(A)
eps = 1.0e-8
n = size(A,1)
U = zeros(n,n)
for k = 1:n
    if k==1 then
            t = A(k,k)
    else 
            t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
    end

    if t <= eps
        printf("Matriz no definida positiva.\n")
        ind = 0
        return
    end
    U(k,k)= sqrt(t)
    for j = k+1:n
        if k==1 then 
                    U(k,j) = A(k,j)/U(k,k)
        else 
                    U(k,j) = ( A(k,j) - U(1:k-1,k)' * U(1:k-1,j) )/U(k,k)
        end
    end
end
ind = 1
endfunction

function x = resolver_chol(A,b)
    [nA,mA] = size(A);
    if nA<>mA then
        error('La matriz A debe ser cuadrada');
        abort;
    end
    n = nA;
    [R,i] = Cholesky(A);
    if i == 0 then
        error("La matriz no es definida positiva")
        abort;
    end
    disp(R);
    y = resolver_MTI(R',b,n);
    x = resolver_MTS(R,y,n);
endfunction

function res = presente(l,x)
    [n,m] = size(l)
    res = 0
    for i = 1:n
       if l(i) == x then
           res = 1
       end
    end
endfunction

// Permuta las filas de A hasta encontrar aquella en donde
// Ningun elemento de la diagonal sea nulo
// FUNCIONA PERO NO SE SI ES EXHAUSTIVA TODAVIA NO ENCONTRE UNA MATRIZ EN DONDE NO FUNCIONE
function [B,P] = reacomodar(A)
    [n,m] = size(A);
    if n<>m then
        error('La matriz A debe ser cuadrada');
        abort;
    end
    
    P = eye(n,n);
    pos = zeros(n,1);
    D = eye(A).*A
    while det(D) == 0
        for i = 1:n
            C = A(:,i); // Copiamos la columna i de A en C
            ban = 0;
            while (ban == 0) && (max(abs(C)) > 0) // Mientras su maximo de C no sea 0
                [m,p] = max(abs(C)); // m = maximo de C, p es la posicion de este maximo
                if A(i,p) == 0 then
                    // Si el intercambio de filas anula el elemento de la diagonal
                    // de la fila que queremos intercambiar entonces eliminamos ese elemento
                    C(p) = 0; 
                else
                    ban = 1;
                end
            end
            A = intercambiar_filas(A,i,p);
            P = intercambiar_filas(P,i,p);
        end
        D = eye(A).*A
    end
    B = A
endfunction

function [L,D,U] = separar_matriz_LDU(A)
    D = eye(A).*A
    n = size(A,1)
    L = zeros(n,n)
    
    for i = 1:n
        for j = i+1:n
            L(j,i) = A(j,i)
        end
    end
    U = A - L - D
endfunction

function [x,P] = jacobi(A,b,x0,eps,maxiter)
    [nA,mA] = size(A);
    if nA<>mA then
        error('La matriz A debe ser cuadrada');
        abort;
    end
    if det(A) == 0 then
        error('La matriz A es singular');
        abort;
    end
    
    n = nA
    D = eye(A).*A
    //Reacomodamos A para que no tenga elementos nulos
    //en la diagonal
    P = eye(n,n);
    if det(D) == 0 then
        [B,P] = reacomodar(A)
        A = B
        b = P*b
        D = eye(B).*B
    end

    I = eye(n,n);
    M = inv(D);  
    T = (I-M*A);   
    
    p = max(abs(spec(T)))
    if p > 1 then
        error("El metodo diverge")
        abort;
    end  
    
    //Primera iteracion
    x = x0

    z = M*b + T*x
    
    iter = 1
    while (norm(z - x) > eps) &&(iter < maxiter)
        x = z
        z = (M*b + T*x)
        iter = iter + 1
    end
    
    x = z
    disp(iter)
endfunction

function x = gauss_seidel(A,b,x0,eps,maxiter)
    tic()
    [n,m] = size(A);
    if n<>m then
        error('La matriz A debe ser cuadrada');
        abort;
    end
    if det(A) == 0 then
        error('La matriz A es singular');
        abort;
    end
    
    D = eye(A).*A
    //Reacomodamos A para que no tenga elementos nulos
    //en la diagonal
    P = eye(n,n);
    if det(D) == 0 then
        [B,P] = reacomodar(A)
        A = B
        b = P*b
        D = eye(B).*B
    end
    
    [L,D,U] = separar_matriz_LDU(A);
    
    N = L + D;
    M = inv(N);
    I = eye(n,n);
    T = (I-M*A);
    
    p = max(abs(spec(T)))
    if p > 1 then
        error("El metodo diverge");
        abort;
    end  
    
    //Primera iteracion
    x = x0;

    z = M*b + T*x;
    
    iter = 1;
    while (norm(z - x) > eps) &&(iter < maxiter)
        x = z;
        z = M*b + T*x; 
        iter = iter + 1;
    end
    
    x = z;
    disp(iter);
    t = toc();
    disp(t);
endfunction

function x = metodo_relajacion(A,b,x0,w,eps,maxiter)
    [n,m] = size(A);
    if n<>m then
        error('La matriz A debe ser cuadrada');
        abort;
    end
    if det(A) == 0 then
        error('La matriz A es singular');
        abort;
    end
    
    D = eye(A).*A
    //Reacomodamos A para que no tenga elementos nulos
    //en la diagonal
    P = eye(n,n);
    if det(D) == 0 then
        [B,P] = reacomodar(A)
        A = B
        b = P*b
        D = eye(B).*B
    end
    
    //Iteraciones
    iter = 1;
    x = zeros(n,1)
    //primera iteracion
    
    //i = 1
    g = A(1,2:n)*x0(2:n)
    x(1) = (1 - w)*x0(1) + (w / A(1,1))*(b(1) - g);
    
    // 1 < i < n
    for i = 2:n-1
        g1 = A(i,1:i-1)*x(1:i-1)
        g2 = A(i,i+1:n)*x0(i+1:n)
        x(i) = (1-w)*x0(i) + (w/A(i,i))*(b(i)-g1-g2);
    end
    
    // i = n 
    g = A(n,1:n-1)*x(1:n-1)
    x(n) = (1 - w)*x0(n) + (w / A(n,n))*(b(n) - g);
    
    z = x
    x = x0
    
    //Inician las iteraciones
    while (norm(z - x) > eps) && (iter < maxiter)
       x = z
       
       //i = 1
       g = A(1,2:n)*x(2:n)
       z(1) = (1 - w)*x(1) + (w / A(1,1))*(b(1) - g);
       
       // 1 < i < n
       for i = 2:n-1
            g1 = A(i,1:i-1)*z(1:i-1)
            g2 = A(i,i+1:n)*x(i+1:n)
            z(i) = (1 - w)*x(i) + (w / A(i,i))*(b(i) - g1 - g2);
       end
        
        //i = n
        s = 0
        g = A(n,1:n-1)*z(1:n-1)
        z(n) = (1 - w)*x(n) + (w / A(n,n))*(b(n) - g);
        
        iter = iter + 1;
    end
    disp(iter);
    x = z;
endfunction

function x = metodo_SOR(A,b,x0,eps,maxiter)
    //Este metodo funciona optimamente si A es Tridiagonal
    [n,m] = size(A);
    if n<>m then
        error('La matriz A debe ser cuadrada');
        abort;
    end
    D = eye(A).*A;
    I = eye(n,n);
    
    T = I - inv(D) * A;
    p = max(abs(spec(T)));
    
    w = 2 / (1 + sqrt(1 - p**2));
    
    //Iteraciones
    iter = 1;
    x = zeros(n,1)
    //primera iteracion
    
    //i = 1
    g = A(1,2:n)*x0(2:n)
    x(1) = (1 - w)*x0(1) + (w / A(1,1))*(b(1) - g);
    
    // 1 < i < n
    for i = 2:n-1
        g1 = A(i,1:i-1)*x(1:i-1)
        g2 = A(i,i+1:n)*x0(i+1:n)
        x(i) = (1-w)*x0(i) + (w/A(i,i))*(b(i)-g1-g2);
    end
    
    // i = n 
    g = A(n,1:n-1)*x(1:n-1)
    x(n) = (1 - w)*x0(n) + (w / A(n,n))*(b(n) - g);
    
    z = x
    x = x0
    
    //Inician las iteraciones
    while (norm(z - x) > eps) && (iter < maxiter)
       x = z
       
       //i = 1
       g = A(1,2:n)*x(2:n)
       z(1) = (1 - w)*x(1) + (w / A(1,1))*(b(1) - g);
       
       // 1 < i < n
       for i = 2:n-1
            g1 = A(i,1:i-1)*z(1:i-1)
            g2 = A(i,i+1:n)*x(i+1:n)
            z(i) = (1 - w)*x(i) + (w / A(i,i))*(b(i) - g1 - g2);
       end
        
        //i = n
        s = 0
        g = A(n,1:n-1)*z(1:n-1)
        z(n) = (1 - w)*x(n) + (w / A(n,n))*(b(n) - g);
        
        iter = iter + 1;
    end
    disp(iter);
    x = z;
endfunction

// UNIDAD 7 ----------------------------------------------------------

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

// DD_Newton(p,x,y) Calcula el resultado de aproximar
// la imagen del punto p utilizando el polinomio interpolante
// resultante del metodo de diferencias divididas de newton
// en los puntos (x(i),y(i))
function w = DD_Newton(p,x,y)
    w = 0
    n = length(x)
    for j=n:-1:2
        w = (w + DD(x(1:j),y(1:j)))*(p-x(j-1))
    end
    w = w + y(1)
endfunction

// Método de Diferencias Divididas de Newton
// Retorna el polinomio interpolante resultante
// de aplicar el método de Diferencias Divididas de Newton
// en los puntos(x(i),y(i))
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
function w = DD(x,y)
    // Entrada: x,y = vectores de puntos de interpolación (x,y)
    // Salida: w = diferencias divididas en los vectores x e y
    n = length(x)
    if n==2 then
        w = (y(n)-y(1))/(x(n)-x(1))
    else
        w = (DD(x(2:n),y(2:n))-DD(x(1:n-1),y(1:n-1)))/(x(n)-x(1))
    end
endfunction

// Error de interpolación
function w = err(p,x,cot)
    // Entrada: p = valor real, x = nodos de interpolación, cot = cota de |f^(n))|
    // Salida: w = error de interpolación en x = p
    n = length(x)
    w = cot/(factorial(n))
    for i=1:n do
        w = w*abs(p - x(i))
    end
endfunction

// Aproximación polinomial de mínimos cuadrados polinomial para matrices con rango completo
function [p,err] = MinCuad_pol(A,b)
    // Entrada: b = vectores 1xn
    // Salida: p = polinomio de mínimos cuadrados; err = vector de errores (eps = Ax-b)
     [w,a] = gausselimPP((A')*A,(A')*(b'))
     p = poly(w,"x","coeff")
     err = A*w-b'
endfunction


// Matriz del método de mínimo cuadrados polinomial
function A = A_mc(x,n)
    // Entrada: x,y = vectores 1xn; n = grado de la aproximación
    // Salida: A = matriz del método de mínimo cuadrados
    m = length(x)
    A = ones(m,1)
    for j=2:(n+1) do
        A = [A,(x').^(j-1)]
    end
endfunction

function p = transformar_pol(p,k)
    c = coeff(p)
    m = length(c)
    for i = 2:m
        c(i) = c(i) * k**(i-1)
    end
    p = poly(c,'x','coeff')
endfunction

// Chebyshev_intervalo
// Retorna un polinomio de grado n
// el cual todas sus raices se encuentran
function x = raices_Chebyshev_intervalo(a,b,n)
    p = Chebyshev_mon(n)
    t = roots(p)
    if (a<>-1)&&(b<>1) then
        for i = 1:n
            x(i) = ((b+a)+t(i)*(b-a))/2
        end
    else
        x = t
    end
    
endfunction

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

function p = Chebyshev_mon(n)
    p = Chebyshev(n)* 1/2**(n-1)
endfunction 

function p = DD_newtonP_chev(f,a,b,n)
    deff("y = F(x)","y="+f);
    x = raices_Chebyshev_intervalo(a,b,n+1)
    y = F(x)
    p = DD_Newton_pol(x,y)
endfunction

// graf_error_int(f,w,a,b)
// f = funcion; w = polinomio interpolante de f
// a,b = extremos de interpolacion
// n = Grado del polinomio interpolante
// Dados una funcion, su polinomio interpolante
// y los extremos de interpolacion se grafica el error de interpolacion
function graf_error_int(f,w,a,b,n)
    deff("y = F(x)","y="+f);
    z = (a-1):0.01:(b+1)
    c = 2**(-2*(n-1))
    deff("v = G(x)","v = F(x) - horner(w,z)")
    plot2d(z,G(z),rect = [(a-1), -c, (b+1), c],leg="f(x) - P(x)")
    a=gca()
    a.x_location = "origin"
    a.y_location = "origin"
endfunction

// Pr7 ej 1
// x = [0 0.2 0.4 0.6]
// y = [1.0 1.2214 1.4918 1.8221]
// w = Lagrange(1/3,x,y)
// disp("Resultado Lagrange")
// disp(w)
// w = DD_Newton(1/3,x,y) 
// disp("Resultado Newton")
// disp(w)


// Ej 7
//x = [0 0.15 0.31 0.5 0.6 0.75]
//y = [1 1.004 1.31 1.117 1.223 1.422]

//disp("(#) n=1.")
//A = A_mc(x,1)
//deter = det((A')*A)
//disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
//disp("El polinomio de mínimos cuadrados de grado 1 es:")
//[p1,err1] = MinCuad_pol(A,y)
//disp(p1)

//disp("(#) n=2.")
//A = A_mc(x,2)
//deter = det((A')*A)
//disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
//disp("El polinomio de mínimos cuadrados de grado 2 es:")
//[p2,err2] = MinCuad_pol(A,y)
//disp(p2)

//disp("(#) n=3.")
//A = A_mc(x,3)
//deter = det((A')*A)
//disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
//disp("El polinomio de mínimos cuadrados de grado 3 es:")
//[p3,err3] = MinCuad_pol(A,y)
//disp(p3)

//disp("(#) Analizamos los errores err = norm(Ax-y,2).")
//disp("Para la aproximación lineal: "+string(norm(err1,2)))
//disp("Para la aproximación cuadrática: "+string(norm(err2,2)))
//disp("Para la aproximación cúbica: "+string(norm(err3,2)))
//disp("Podemos decir que es mejor la aproximación cúbica en este caso.")

// Ejercicio 10

// UNIDAD 8 -------------------------------------------------


//Regla del trapecio
function ap = Regla_del_trapecio(f,a,b)
    deff("y = F(x)","y="+f);
    h = (b - a)
    ap = (F(a) + F(b)) * (h/2)
endfunction

function ap = RT_Compuesta(f,a,b,n)
    deff("y = F(x)","y="+f);
    ap = 0
    h = (b-a)/n
    ap = ap + (1/2)*(F(a) + F(b))
    for i = 1:n-1
        ap = ap + F(a+i*h)
    end
    ap = ap * h
endfunction

//Regla de Simpson
function ap = simpson(f,a,b)
    deff("y = F(x)","y="+f);
    h = (b-a)/2
    ap = (h/3)*(F(a) + 4*(F(a+h)) + F(b))
endfunction

function ap = simpson_compuesta(f,a,b,n)
    if(int(n/2) <> (n/2)) then
        disp("La cantidad de intervalos no es par")
        abort;
    end
    deff("y = F(x)","y="+f);
    h = (b-a)/n
    ap = 0
    ap = ap + F(a) + F(b)
    
    for i = 1:2:n-1
        ap = ap + 4*F(a+i*h)
    end
    for i = 2:2:n-2
        ap = ap + 2* F(a+i*h)
    end
    ap = ap*(h/3)
    
endfunction

