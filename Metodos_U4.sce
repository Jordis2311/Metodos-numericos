
function B = intercambiar_filas(A,i,j)
    [n,m] = size(A);
    aux = A(i,:);
    A(i,:) = A(j,:);
    A(j,:) = aux;
    B = A;
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

function x = gauss_Jordan(A,b)
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
    n = nA;    // Tamaño de la matriz
    
    // Eliminación progresiva con pivoteo parcial
    for k=1:n-1
        kpivot = k; amax = abs(a(k,k));  //pivoteo
        for i=k+1:n
            if abs(a(i,k))>amax then
                kpivot = i; amax = a(i,k);
            end;
        end;
        temp = a(kpivot,:); a(kpivot,:) = a(k,:); a(k,:) = temp;
        
        for i=k+1:n
            for j=k+1:n+1
                a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
            end;
            for j=1:k        // no hace falta para calcular la solución x
                a(i,j) = 0;  // no hace falta para calcular la solución x
            end              // no hace falta para calcular la solución x
        end;
    end;
    
    //Eliminacion regresiva (transformamos A en la matriz identidad)
    for k = n:-1:2
       a(k,:) = a(k,:) / a(k,k); //Divide toda la fila por el elemento de la diagonal
       
       for j = k-1:-1:1
           a(j,n+1) = a(j,n+1) - a(k,n+1) * a(j,k);
           a(j,k) = 0
       end
    end
    a(1,:) = a(1,:) / a(1,1);SS
    
    x(:) = a(:,n+1);
endfunction

function [x,a] = gausselimPP(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana con pivoteo parcial.

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
n = nA;    // Tamaño de la matriz

// Eliminación progresiva con pivoteo parcial
for k=1:n-1
    kpivot = k; amax = abs(a(k,k));  //pivoteo
    for i=k+1:n
        if abs(a(i,k))>amax then
            kpivot = i; amax = a(i,k);
        end;
    end;
    temp = a(kpivot,:); a(kpivot,:) = a(k,:); a(k,:) = temp;
    
    for i=k+1:n
        for j=k+1:n+1
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
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
    end;
    x(i) = (a(i,n+1)-sumk)/a(i,i);
end;
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

function [L,U,P] = factorizar_PLU(A)
    [nA,mA] = size(A);
    if nA<>mA then
        error('La matriz A debe ser cuadrada');
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

function [U,ind] = choleskyV2(A)
// Factorización de Cholesky.
// Trabaja únicamente con la parte triangular superior.
//
// ind = 1  si se obtuvo la factorización de Cholesky.
//     = 0  si A no es definida positiva
//
//******************
eps = 1.0e-8
//******************

n = size(A,1)
U = zeros(n,n)

t = A(1,1)
if t <= eps then
    printf('Matriz no definida positiva.\n')
    ind = 0
    return
end
U(1,1) = sqrt(t)
for j = 2:n
    U(1,j) = A(1,j)/U(1,1)
end
    
for k = 2:n
    t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
    if t <= eps then
        printf('Matriz no definida positiva.\n')
        ind = 0
        return
    end
    U(k,k) = sqrt(t)
    for j = k+1:n
        U(k,j) = ( A(k,j) - U(1:k-1,k)'*U(1:k-1,j) )/U(k,k)
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


function [Q, R] = QR(A)
    [mA,nA] = size(A);
    if mA<nA then
        error('Columnas no linealmente independientes');
        abort;
    end
    
    m = mA;
    n = nA;
    Q = zeros(m,n);
    R = zeros(n,n);
    
    for k = 1:n
        if k == 1 then
            v = norm(A(1:m,1));
            Q(1:m,1) = A(1:m,1) / v;
        else
            // Establecemos el vector a restar
            vect = zeros(1,m)';
            for i = 1:k-1
                vect = vect + ((A(1:m,k)'*Q(1:m,i)) * Q(1:m,i));
            end
            // Calculamos el v y la columna de Q correspondientes
            v = norm(A(1:m,k) - vect);
            Q(1:m,k) = (A(1:m,k) - vect) / v;
        end
        // Despues de esto tenemos el v y la columna de Q correspondientes a la iteracion
        // Obtengamos la fila correspondiente de R
        for i = k:n
            if i == k then
                R(k,i) = v;
            else
                R(k,i) = A(1:m,i)'*Q(1:m,k);    
            end
        end
    end
endfunction
