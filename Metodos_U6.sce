function p = pol_car(A)
   x = poly([0 1], "x", "coeff")
   n = size(A, 1)
   I = eye(n,n)
   p = det(A - x * I)
endfunction

function circ(r,x,y);
    xarc(x-r,y+r,r*2,r*2,0,360*64)
endfunction

function gers(A)
    n = size(A,1)
    c = diag(A)
    r = sum(abs(A), 'c') - abs(c)
    Mx = round(max(c + r) + 1)
    mx = round(min(c - r) - 1)
    My = round(max(r) + 1)
    my = -My
    rect = [mx,my,Mx,My]
    plot2d(real(spec(A)),imag(spec(A)),-1,"031","",rect)
    replot(rect)
    xgrid()
    for i=1:n
        circ(r(i), c(i), 0)
    end
endfunction

function gers2(A)
    [n,m] = size(A);
    centros = diag(A);
    radios = sum(abs(A),'c') - abs(centros) ;
    
    // buscamos calcular un rectángulo que contenga a todos los circulos
    // esquina inferior izquierda
    
    mx = round (min(centros - radios)-1);
    my = round (min(-radios)-1);
    
    // esquina superior derecha
    
    Mx = round(max(centros+radios)+1);
    My = round(max(radios)+1);
    
    rectangulo = [mx my Mx My];
    
    // dibujamos los autovalores
    plot2d(real(spec(A)),imag(spec(A)),-2,"031","",rectangulo)
    replot(rectangulo); // reeplaza al rect
    xgrid();
    
    for i=1:n
        circ(radios(i),centros(i),0)
    end
    
endfunction

function circGersValor_acotador(A)
    gers(A);
    gers(A');
endfunction

function [z,l] = metodo_potencia(A,x0,eps,miter)
    [n,m] = size(A)
    if n <> m then
        disp("La matriz debe ser cuadrada");
        abort
    end
    z = x0
    w = A*z;
    z_ant = z;
    z = w / norm(w);
    
    iter = 1;
    while (norm(z-z_ant) > eps) && (iter < miter)
        iter = iter + 1;
        w_ant = w;
        z_ant = z;
        w = A*z;
        z = w / norm(w);
    end
    [a,p] = max(abs(w_ant));
    if a == 0 then
        disp("La matriz es singular");
        l = 0;
    else
        l = w(p) / z_ant(p)
    end
    z = w / norm(w);
endfunction

