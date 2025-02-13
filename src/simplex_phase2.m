function [x,fx,lambda,m,itrab,iter,ind]=simplex_phase2(c,A,b,n,nI,nD,x,itemax,imp)
%        [x,fx,lambda,m,itrab,iter,ind]=simplex_phase2(c,A,b,n,nI,nD,x,itemax,imp)
%
%*****************************************************************
% OBJETIVO:                                             *
%    minimizar c^Tx, siendo x variables reales sujetas a:        *
%                      (2)  A(:,i)'* x = b(i), i=1,...,nI        *
%                      (3)  A(:,i)'* x <= b(i), i=nI+1,...,nI+nD *
%    dado x  verificando (2) y (3)                               *
% AUTOR: PABLO GARCIA PEREZ                                      *
% FECHA: 10/05/2021                                              *
% VERSION: PRELIMINAR                                            *
% ****************************************************************
%
% SE READAPTAN LOS FACTORES QR (llamando a qr_add_column.m y qr_remove_column.m)
% ***OJOOOOO**: LAS RESTRICCIONES DE IGUALDAD SERAN LINEALMENTE 
%               INDEPENDIENTES Y NO SE TRATAN LAS RESTRICCIONES 
%               DE COTA SOBRE LAS VARIABLES SEPARADAMENTE DE LAS
%               DE LAS RESTRICCIONES DE DESIGUALDAD GENERALES 
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
% EN SALIDA:
%  x      Ultimo iterante calculado.
%
%  fx     Valor de la funcion objetivo en x, si ind=0.
%
% lambda Vector de multiplicadores de Lagrange, si ind=0.
%
% m       Numero de restricciones en el conjunto de trabajo final.
%
% trab   Vector que indica, en las m primeras componentes, las restricciones 
%         del conjunto de trabajo en x, de la siguiente forma:
%         trab(j)=p indica que la j-esima restriccion de trabajo  
%         es la p-esima restriccion del problema.
%
% iter  Numero de iteraciones realizadas.
%
% ind     Indicador del exito o fracaso del proceso con los valores:
%              0: Se ha encontrado una solucion (global)
%              1: Se ha superado el numero maximo de iteraciones permitido
%             -1: Dificultades en la readaptacion de los factores QR
%              2: Problema no acotado
%             -2: Las restricciones de igualdad son linealmente depenientes
%              3: Demasiadas iteraciones consecutivas con degeneracion
%---------------------------------------------------------------------

%Se determinan algunas constantes y se inicializan algunas variables
 fx=[]; lambda=[]; m=0; itrab=zeros(n,1); iter=0; ind=[];
 eps1=eps^0.75;  %tolerancia para el test de parada

 idege=0;    %contabiliza las iteraciones CONSECUTIVAS con degeneracion
 idegemax=2*nD;    %Pone un tope a la degeneracion consecutiva

%Se comienza calculando la factorizacion qr y los vectores 
% iredes e itrab asociados al punto inicial ADMISIBLE
%SE SUPONE LA INDEPENDENCIA LINEAL DE LAS RESTRICCIONES DE IGUALDAD

%Comenzamos con las de igualdad que por definicion son activas y por
%hipotesis son linealmente independientes:

if nI>0     %puede que el problema no tenga de igualdad
   [q,r] = qr(A(:,1:nI));
   itrab(1:nI)=(1:nI)';
else        %nI=0, no hay restricciones de igualdad
   q=eye(n); r=[];
end

m=nI;   %numero de restricciones en el conjunto de trabajo

%Se a�aden al conjunto de trabajo las restricciones de desigualdad
%activas linealmente independientes

iredes=zeros(nD,1);     %ninguna de las restriccinoes de desigualdad esta
                        %en el conjunto de trabajo, de momento

%recorremos las restricciones de desigualdad y miro si son li.
i=1;

while(m<n && i<=nD)
    if abs(A(:,nI+i)'*x-b(nI+i)) <= eps1           %miramos si es activa
        [q,r,ind01]=qr_add_column(q,r,A(:,nI+i),n,m,0);
        if ind01==0             %restricciones linealmente independientes
            m=m+1;              %la a�adimos al conjunto de trabajo
            itrab(m)=nI+i;      
            iredes(i)=1;
        end
    end
    i=i+1;
end

%Comienzan las iteraciones

if imp==2
    fprintf('\n Comienzan las iteraciones.')
    fprintf('\n Restricciones en el conjunto de trabajo:')
    fprintf('%3d',itrab(1:m))    
end

for iter=1:itemax
     %Calculo del gradiente reducido (gp)
     if m>0 && m<n
         gp=q(:,m+1:n)'*c;
     
     %Teniendo en cuenta los casos extremos al programar
     elseif m==0     %q=eye(n)
         gp=c;
     else
         gp=0;       %m=n
     end
     
     %Paso (2):
     if (m<n && norm(gp)/(1+abs(c'*x))>eps1)     %el gradiente reducido 
                                                 %NO es "=0",
       %se calcula una direccion de descenso
       if m>0
           d=-q(:,m+1:n)*gp;
       else %m=0, luego Q=Z y Z*Z'=Id por ser ortogonal
           d=-c;
       end
       
     else    %Si gp"=0", mirar los candidatos a multiplicadores (m=n)
       lambda=-r(1:m,1:m)\(q(:,1:m)'*c);
       %Se determina si x es solucion del PL
       if m>nI
           [minlagdes,iCOLdel]=min(lambda(nI+1:m));
       else
           minlagdes=1000000;
       end
   
       if minlagdes >= -eps1
          fprintf('\n \n SE HA ENCONTRADO UNA SOLUCION (GLOBAL): \n');
          fprintf('\t \t \t %13.5e\n',x);
          fx=c'*x;
          fprintf('\n VALOR DE LA FUNCION: %13.5e',fx);
          fprintf('\n VECTOR DE MULTIPLICADORES:\n')
          fprintf('\t \t \t %13.5e\n',lambda(1:m));
          ind=0;
          return
       else 
           
         %Si no se detecta SOLUCION, se elimina una restriccion del 
         %conjunto de trabajo y se calcula una direccion de descenso
         
         iDESdel=itrab(nI+iCOLdel)-nI;  %restriccion de desigualdad
                                        %que se elimina
         
          if imp==2
             fprintf('\n Se elimina la restriccion %3d', nI+iCOLdel);
          end
      
          iredes(iDESdel)=0;
          itrab(nI+iCOLdel)=[];
          [q,r]=qr_remove_column(q,r,n,m,nI+iCOLdel);
          m=m-1;
          d=-sign(A(:,nI+iCOLdel)'*q(:,m+1))*q(:,m+1);
       end
     end
    
     %Se calcula el desplazamiento (o paso)
     ro=realmax;
      
     for i=1:nD
         wd= A(:,nI+i)'*d;
         if (iredes(i)==0 &&  wd > eps1)
             roi=(b(nI+i)-A(:,nI+i)'*x)/wd;
             if roi < ro
                 ro=roi;
                 indro=nI+i;
             end
         end
     end

     
     %Test sobre la acotacion del problema
     if ro == realmax
         if imp>=1
             fprintf('\n \n PROBLEMA NO ACOTADO: \n');
         end
         ind= 2;
         return
     end
     
     %Se estudia la degeneracion
     if abs(ro) <= eps1
         idege=idege+1;
     else
         idege=0;
     end
     
     if idege >= idegemax
         fprintf('\n \n Demasiadas iteraciones CONSECUTIVAS con degeneracion')
        ind= 3;
        return
     end
     
     
     %Se a�ade una restriccion al conjunto de trabajo
     
     [q,r,ind01]=qr_add_column(q,r,A(:,indro),n,m,0);
     if ind01==0     %restricciones linealmente independientes
                     %la a�adimos al conjunto de trabajo
         m=m+1;
         itrab(m)=indro;
         iredes(indro-nI)=1;
     end

     
     if imp==2
        fprintf('\n Se a�ade la restriccion itrab(p)=%3d',itrab(m));   
     end
     
     %Se readapta el aproximante a la solucion
     x=x+ro*d;
end
  
  if (imp>=1 && iter==itemax)
      fprintf('\n \n Se ha superado el maximo numero de iteraciones permitido')
  end
  
  ind=1;
  
end