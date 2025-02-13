function [q,r]=qr_remove_column(q,r,n,m,icol)
%        [q,r]=qr_remove_column(q,r,n,m,icol)
%***************************************************
%OBJETIVO: Se readapta la factorizacion QR de una  *
%matriz cuando a esta se le quita una columna.     *
%AUTORA: C. Pola, Universidad de Cantabria.        *
%***************************************************
%LISTA DE LLAMADA
%DE ENTRADA
%  q   : contiene a la matriz "Q" inicial
%  r   : contiene a la matriz "R" inicial
%  n   : numero de filas de la matriz a factorizar
%  m   : numero de columnas de la matriz "R" inicial
% icol : numero de la columna que se desea quitar
%DE SALIDA
%  q   :contiene a la nueva matriz "Q" 
%  r   :contiene a la nueva matriz "R" 
%

%
% Si la columna que se desea quitar ocupa el ultimo lugar
% en la matriz se termina el proceso
%
 if m==icol,
    r(:,[m])=[];
    return; 
 end
%
% Se modifican las ultimas columnas de R y la matriz Q
%
 for i=icol+1:m,
    if r(i,i)~=0,
       i1=i-1;i11=i+1;
       a=norm(r(i1:i,i),2);
       if a>eps,
 	       if r(i1,i)~=0, a=sign(r(i1,i))*a; end
	       r(i1:i,i)=r(i1:i,i)/a;
	       r(i1,i)=r(i1,i)+1;
	       s1=r(i1,i);s2=r(i,i);s=s2/s1;
	       for j=1:n,
	          t=-q(j,i1)-q(j,i)*s;
	          q(j,i1)=q(j,i1)+t*s1;
 	          q(j,i)=q(j,i)+t*s2;
          end
	       for j=i11:m,    
	          t=-r(i1,j)-r(i,j)*s;
   	       r(i1,j)=r(i1,j)+t*s1;
	          r(i,j)=r(i,j)+t*s2;
 	       end
   	    r(i1,i)=-a;
       end
    end
 end
%
% Se coloca en las m-1 primeras columnas de la matriz R esultante
%
 for j=icol+1:m, 
     j1=j-1;
     r(1:j1,j1)=r(1:j1,j);
 end
 r(:,[m])=[];
 return
