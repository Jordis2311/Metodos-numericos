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

function z=hornet(x,p)
    if  degree(p) <= 0 then
        z = 0
    else
        ai = coeff(p)
        n = degree(p)
        b(n+1) = ai(n+1)
        for i = n : -1 : 1
            b(i) = ai(i) + x*b(i+1);  
        end
        z = b(1)
    end
endfunction


function valor = tay(fun, v, n, a, h)
    
    valor = 0;
    
    for i = 0:n
        valor = valor + ( (1 / factorial(i)) * dev_m(fun,i,a,h) *(v - a)**i);
    end
    
endfunction
