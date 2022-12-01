
// biseccion(fun,a0,b0,tol) Metodo de la biseccion
// fun: string de una funcion continua
// a0,b0: intervalo inicial donde queremos encontra la raiz
// tol: Tolerancia usada en el criterio de parada
// Calcula una estimacion raiz de la funcion en un intervalo 
// donde fun(a0) * fun(b0) < 0 utilizando el metodo de la biseccion
// Ejemplo
// biseccion('%e**x - 2',0,1,0.001) = 0.6923828
function salida = biseccion(fun,a0,b0,tol)
    deff("y = f(x)","y="+fun);
    iter = 1;
    a = a0;
    b = b0;
    c = (a + b) / 2;
    
    if f(a)*f(b) >= 0 then
        disp("Intervalo no valido");
        abort;
    end
    band = 1
    while (abs(b - c) > tol) && (band == 1)
        iter = iter + 1;
        if f(c) == 0 then
            band = 0
        else
            if (f(c)*f(a)) < 0 then
                b = c;
                c = (a + b) / 2; 
            else
                if (f(c)*f(b)) < 0 then
                    a = c;
                    c = (a + b) / 2;
                end
            end
        end
    end
    
    err = (1/2)**iter * abs(b0 - a0)
    disp('Se alcanzo la tolerancia con una cota de error de: '+ string(err))
    disp('Iteraciones realizadas: ');
    disp(iter);
    salida = c;  
 
endfunction

// newton(fun,x0,tol,iter) Metodo de Newton
// fun: string de una funcion continua
// x0: estimacion inicial de la raiz
// tol: Tolerancia usada en el criterio de parada
// iter: maximo de iteraciones que se realizaran
// Calcula una estimacion de la raiz de f con una tolerancia dada
// Utilizando el metodo de newton
// Ejemplo:
// newton(fun,x0,tol,iter) = 0.6931476
function salida = newton(fun,x0,tol,iter)
    deff("y = f(x)","y="+fun);
    i = 1;
    x1 = x0 - (f(x0) / numderivative(f,x0));
    while (abs(x1-x0) > tol) && (i<iter)
        i = i+1;
        x0 = x1;
        x1 = x0 - (f(x0) / numderivative(f,x0));
    end
    
    disp('Iteraciones realizadas: ');
    disp(i);
    salida = x1;
    
endfunction 

// secante(fun,x0,y0,tol,maxiter) Metodo de la secante
// fun: string de una funcion continua
// x0,y0: estimaciones iniciales de la raiz
// tol: Tolerancia usada en el criterio de parada
// maxiter: maximo de iteraciones que se realizaran
// Calcula una aproximacion a la raiz de la funcion
// utilizando el metodo de la secante
function salida = secante(fun,x0,y0,tol,maxiter)
    deff("y = f(x)","y="+fun);
    i = 1;
    x = x0;
    y = y0;
    c = x - (f(x) * (x - y) / (f(x) - f(y)));
    while (abs(x-y) > tol) && (i < maxiter)
        i = i + 1;
        y = x;
        x = c;
        c = x - f(x) * ((x - y) / (f(x) - f(y)));
    end
    disp('Iteraciones realizadas: '+ string(i));

    salida = c;
endfunction

function salida = secante2(fun, x0, x1, tol)
	deff("y=f(x)", "y=" + fun)
	while abs(x0 - x1) > tol
		x2 = x1 - f(x1)*((x1-x0)/(f(x1)-f(x0)))
		x0 = x1
		x1 = x2
	end
	salida = x0
endfunction

// falsa_posicion(fun,a1,b1,tol) Metodo de la Falsa posicion
// fun: string de una funcion continua
// a1,b1: Intervalo inicial tal que fun(a)*fun(b) < 0
// tol: Tolerancia usada en el criterio de parada
// Calcula una estimacion de la raiz de la funcion utilizando
// el metodo de la falsa posicion
// Ejemplo
// falsa_posicion('%e**x -2',0,1,0.001) = 0.6931472
function salida = falsa_posicion(fun,a1,b1,tol)
    deff("y=f(x)", "y=" + fun)
    if (f(a1)*f(b1)) > 0 then
        disp("intervalo no valido")
        abort;
    end
    b = b1
    a = a1
    c = b - f(b) *((b - a) / (f(b)-f(a)))
    iter = 1
    ban = 1
    while(abs(b-c) > tol) && ban
       if f(c) == 0 then
            ban = 0
       else if(f(a)*f(c) < 0) then
            b = c
            c = b - f(b) *((b - a) / (f(b)-f(a)))
            else if (f(b)*f(c) < 0) then
                    a = c
                    c = b - f(b) *((b - a) / (f(b)-f(a)))
                 end
            end
       end
    end
    
    salida = c
endfunction

// punto_fijo(g,x0,tol,maxiter) Metodo iterativo de punto fijo
// g: string de la funcion continua
// x0: Estimacion inicial del punto fijo
// tol: tolerancia utilizada para el criterio de parada
// maxiter: maximo de iteraciones realizables
// Calcula una estimacion de un punto fijo de la funcion
// utilizando el metodo de punto fijo
function x = punto_fijo(g,x0,tol,maxiter)
    deff("y=G(x)", "y=" + g)
    x1 = G(x0)
    iter = 1
    while (norm(x1-x0) > tol) && (iter < maxiter)
        x0 = x1
        x1 = G(x0)
        iter = iter + 1
    end
    
    x = x1
    
endfunction



function y=fx(x)
    f1 = (x(1)**2) + (x(1)*(x(2)**3))-9
    f2 = (3*(x(1)**2)*x(2))-4-(x(2)**3)
    y = [f1;f2]
endfunction

function res=fy(x0)
    x = x0(1)
    y = x0(2) 
    f1 = 1+(x**2)-y**2+%e**x*cos(y)
    f2 = 2*x*y+%e**x*sin(y)
    res = [f1; f2]
endfunction

function y=newton_mv(f, x0, it)
    deff("y = f(x)","y="+f);
    x1 = x0 - ((numderivative(f, x0)**-1)*f(x0))
    i = 0
    while(i < it)
        x0 = x1
        x1 = x0 - ((numderivative(f, x0)**-1)*f(x0))
        i = i + 1
    end
    y = x1
endfunction

