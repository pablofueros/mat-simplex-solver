function [x,fx,lambda,m,itrab,iter,ind]=simplex_solver(c,A,b,n,nI,nD,x,itemax,imp)
%        [x,fx,lambda,m,itrab,iter,ind]=simplex_solver(c,A,b,n,nI,nD,x,itemax,imp)
%
%********************************************************************
% OBJETIVO:                                                         *
%                                                                   *
% AUTOR: PABLO GARCIA PEREZ                                         *
% FECHA: 20/05/2021                                                 *
% VERSION: PRELIMINAR                                               *
% *******************************************************************
%
% EN ENTRADA:
%  c      Vector columna n-dimensional con los coeficientes de la 
%         funcion objetivo.
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
%  x      Vector columna de dimension n. Contiene un punto ADMISIBLE.
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
%
% ----------------------------------------------------------------------

%Inicializamos algunas variables
eps1=eps^0.75;

%Diferenciamos las tres situaciones posibles:

if length(x) == 0    %No se suministra un punto inicial. Se llama a faseI:
    [x,gamma,ind,ind2]=simplex_phase1(A,b,n,nI,nD,itemax,imp);
    
    %Miramos si el punto calculado es admisible
    if ind == 0                     %x es admisible
        %Llamamos a FaseII partiendo de x
        [x,fx,lambda,m,itrab,iter,ind]=simplex_phase2(c,A,b,n,nI,nD,x,itemax,imp);
    else
        if ind2 == 0                %Si FaseII no ha tenido problemas
            %Miramos condicion optimal sobre gamma
            if abs(gamma) > eps1    %x no es admisible
                x=[];
                ind=1;
                if imp >= 1
                    fprintf('\n El Problema no est� bien puesto');
                end
            else                    %x es pto admisible
                ind=0;
                if imp >= 1
                    fprintf('\n Se ha encontrado un punto admisible: \n');
                    fprintf('\t \t \t %13.5e\n',x);
                end

                %Llamamos a FaseII partiendo de x
                [x,fx,lambda,m,itrab,iter,ind]=simplex_phase2(c,A,b,n,nI,nD,x,itemax,imp);
            end
        else                        %FaseII ha tenido algun problema
            ind=1;
            if imp >= 1
                    fprintf('\n No se ha podido encontrar un punto admisible');
            end
        end
    end
    
else        %Si se suministra punto inicial
    
    %Veamos si x es admisible para el problema
    r=A'*x-b;
    
    if (norm(r(1:nI)) <= eps1 && max(r(nI+1:nI+nD)) <= eps1)  %Si lo es
        [x,fx,lambda,m,itrab,iter,ind]=simplex_phase2(c,A,b,n,nI,nD,x,itemax,imp);
        
    else %Mismo proceso que si x=[];
        [x,gamma,ind,ind2]=simplex_phase1(A,b,n,nI,nD,itemax,imp);
    
        %Miramos si el punto calculado es admisible
        if ind == 0                     %x es admisible
            %Llamamos a FaseII partiendo de x
            [x,fx,lambda,m,itrab,iter,ind]=simplex_phase2(c,A,b,n,nI,nD,x,itemax,imp);
        else
            if ind2 == 0                %Si FaseII no ha tenido problemas
                %Miramos condicion optimal sobre gamma
                if abs(gamma) > eps1    %x no es admisible
                    x=[];
                    ind=1;
                    if imp >= 1
                        fprintf('\n El Problema no est� bien puesto');
                    end
                else                    %x es pto admisible
                    ind=0;
                    if imp >= 1
                        fprintf('\n Se ha encontrado un punto admisible: \n');
                        fprintf('\t \t \t %13.5e\n',x);
                    end

                    %Llamamos a FaseII partiendo de x
                    [x,fx,lambda,m,itrab,iter,ind]=simplex_phase2(c,A,b,n,nI,nD,x,itemax,imp);
                end
            else                        %FaseII ha tenido algun problema
                ind=1;
                if imp >= 1
                        fprintf('\n No se ha podido encontrar un punto admisible');
                end
            end
        end
    end
end

end
    