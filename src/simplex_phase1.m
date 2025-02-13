function [x,gamma,ind,ind2]=simplex_phase1(A,b,n,nI,nD,itemax,imp)
%        [x,gamma,ind,ind2]=simplex_phase1(A,b,n,nI,nD,itemax,imp)
%
%********************************************************************
% OBJETIVO:                                                         *
%    minimizar gamma, siendo x y gamma variables reales sujetas a:  *
%                      (2)  A(:,i)'* x = b(i), i=1,...,nI           *
%                      (3)  A(:,i)'* x - gamma <= b(i),             *
%                           i=nI+1,...,nI+nD                        *
%    dado x  verificando (2) y (3)                                  *
% AUTOR: PABLO GARCIA PEREZ                                         *
% FECHA: 11/05/2021                                                 *
% VERSION: PRELIMINAR                                               *
% *******************************************************************
%
% EN ENTRADA:
%
%  A      En cada columna de la matriz A colocar los coeficientes de una 
%         restriccion. Las primeras nI columnas corresponden a las 
%         restricciones de igualdad.   
%  b      Vector con los terminos independientes de las restricciones. 
%         Las primeras nI coordenadas corresponden a las restricciones 
%         de igualdad.
%  n      numero de variables
% nI      numero de restricciones de igualdad
% nD      numero de restricciones de desigualdad
%
%  itemax Numero maximo de iteraciones a realizar por el algoritmo.
%
%  imp    Nivel de impresion:  
%           0 no hay salida impresa de resultados
%           1 solo se presenta informacion del final del proceso
%           2 con informacion de las iteraciones
%
% EN SALIDA
%
% x       punto inicial admisible para el problema PL
% gamma   valor obtenido de gamma.
%         si es 0, el problema est� mal puesto
% ind     Indicador del exito o fracaso del proceso con los valores:
%         0: Se ha encontrado un punto admisible para el problema
%         1: Se ha encontrado una solucion de PL auxiliar
%
% ind2    0: Se ha encontrado una solucion de PL auxiliar
%         1: Se ha superado el numero maximo de iteraciones permitido
%        -1: Dificultades en la readaptacion de los factores QR
%         2: Problema no acotado
%        -2: Las restricciones de igualdad son linealmente depenientes
%         3: Demasiadas iteraciones consecutivas con degeneracion
% ----------------------------------------------------------------------

%Inicializamos algunas variables:
x=[]; gamma=[]; ind=[];
eps1=eps^0.75;

%Buscamos una solucion del sistema (2)
if nI == 0            %No hay restricciones de igualdad
    x=zeros(n,1);
else
    x=A(:,1:nI)'\b(1:nI);
end

%Ahora buscamos gamma tal que A(:,i)'*x-b(i) <= gamma
gamma=max(A(:,nI+1:nI+nD)'*x-b(nI+1:nI+nD));

if  gamma <= eps1   %Sabemos que x admisible y no llamamos a FaseII
    gamma=0;
    ind=0;
    ind2=0;
    
else
    ind=1;
    %Llamamos a fase II, 
    c=[zeros(n,1); 1];
    A= [A                          zeros(n,1) ;
        zeros(1,nI)  -ones(1,nD)   -1         ];  
    b=[b; 0];
    n=n+1;
    %nI se queda igual
    nD=nD+1;
    x=[x; gamma];
    %itermax las mantenemos igual
    %imp se manteiene igual

    [x,fx,lambda,m,itrab,iter,ind2]=simplex_phase2(c,A,b,n,nI,nD,x,itemax,imp);
    
    gamma=x(n); 
    x=x(1:n-1);
       
end

end
