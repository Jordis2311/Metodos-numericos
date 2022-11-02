// UTILIDADES -------------------------------------------------------------------------------------------
function res = presente(l,x)
    [n,m] = size(l)
    res = 0
    for i = 1:n
       if l(i) == x then
           res = 1
       end
    end
endfunction

function B = intercambiar_filas(A,i,j)
    [n,m] = size(A);
    aux = A(i,:);
    A(i,:) = A(j,:);
    A(j,:) = aux;
    B = A;
endfunction

// Permuta las filas de A hasta encontrar aquella en donde
// Ningun elemento de la diagonal sea nulo
// FUNCIONA PERO NO SE SI ES EXHAUSTIVA TODAVIA NO ENCONTRE UNA MATRIZ EN DONDE NO FUNCIONE
function [B,P] = reacomodar(A)
    [n,m] = size(A)
    
    if det(A) == 0 then
        disp("La matriz es singular");
        abort;
    end
    
    P = eye(n,n) // matriz de permutacion
    
    D = diag(A);
    while prod(D) == 0
        for i = 1:n
           if A(i,i) == 0 then
               C = A(:,i)
               [m,p] = max(abs(C))
               m1 = m // El unico asegurado no nulo
               p1 = p
               ban = 0
               // Buscamos el maximo que no anule el elemento de la diagonal
               // De la fila que queremos intercambiar
               while (ban == 0) && (m > 0) 
                   if A(i,p) <> 0 then
                        aux = A(p,:)
                        A(p,:) = A(i,:)
                        A(i,:) = aux
                        
                        aux = P(p,:)
                        P(p,:) = P(i,:)
                        P(i,:) = aux
                        
                        ban = 1;
                   else
                        C(p) = 0
                        [m,p] = max(abs(C))
                   end
               end
               // Si no se encontro entonces intercambiamos filas con
               // El primer maximo que encontramos (Como A es no singular) no tiene fila nula
               // Por lo que al menos este es no nulo
               if ban == 0 then
                    aux = A(p1,:)
                    A(p1,:) = A(i,:)
                    A(i,:) = aux
                    
                    aux = P(p1,:)
                    P(p1,:) = P(i,:)
                    P(i,:) = aux
               end
           end
        end
        D = diag(A)
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


// METODOS -------------------------------------------------------------------------------------------
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
    D = diag(A);
    //Reacomodamos A para que no tenga elementos nulos
    //en la diagonal
    P = eye(n,n);
    if prod(D) == 0 then
        [B,P] = reacomodar(A)
        A = B
        b = P*b
    end
    
    D = eye(A).*A
    I = eye(n,n);
    M = inv(D);  
    T = (I-M*A);   
    
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
    
    D = diag(A);
    //Reacomodamos A para que no tenga elementos nulos
    //en la diagonal
    P = eye(n,n);
    if prod(D) == 0 then
        [B,P] = reacomodar(A)
        A = B
        b = P*b
    end
    
    [L,D,U] = separar_matriz_LDU(A);
    
    N = L + D;
    M = inv(N);
    I = eye(n,n);
    T = (I-M*A);
    
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
    
    D = diag(A);
    //Reacomodamos A para que no tenga elementos nulos
    //en la diagonal
    P = eye(n,n);
    if prod(D) == 0 then
        [B,P] = reacomodar(A)
        A = B
        b = P*b
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
