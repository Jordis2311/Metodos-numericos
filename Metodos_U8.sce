// Regla_del_trapecio(f,a,b,n) Regla del Trapecio
// f: funcion a integrar
// a,b: Extremos de integracion
// Calcula una aproximacion de la integral de f en los intervalos
// [a,b] mediante la regla del trapecio
// Ejemplo
// Regla_del_trapecio('%e**x',0,1) = 1.8591409
function ap = Regla_del_trapecio(f,a,b)
    deff("y = F(x)","y="+f);
    h = (b - a)
    ap = (F(a) + F(b)) * (h/2)
endfunction

// RT_Compuesta(f,a,b,n) Regla del Trapecio compuesta
// f: funcion a integrar
// a,b: Extremos de integracion
// n: Cantidad de intervalos
// Calcula una aproximacion de la integral de f en los intervalos
// [a,b] mediante la regla del trapecio compuesto
// Ejemplo
// RT_Compuesta('%e**x',0,1,10) = 1.7197135
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

// simpson(f,a,b) Regla de Simpson
// f: funcion a integrar
// a,b: Extremos de integracion
// Calcula una aproximacion de la integral de f en los intervalos
// [a,b] mediante el metodo de simpson
// Ejemplo
// simpson('%e**x',0,1) = 1.7188612
function ap = simpson(f,a,b)
    deff("y = F(x)","y="+f);
    h = (b-a)/2
    ap = (h/3)*(F(a) + 4*(F(a+h)) + F(b))
endfunction

// simpson_compuesta(f,a,b,n) Regla de Simpson Compuesta
// f: funcion a integrar
// a,b: Extremos de integracion
// n: Cantidad de intervalos
// Calcula una aproximacion de la integral de f mediante
// el metodo de simpson compuesto
// Ejemplo
// simpson_compuesta('%e**x',0,1,10) = 1.7182828
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
