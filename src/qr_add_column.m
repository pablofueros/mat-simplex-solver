function [q,r,ind]=qr_add_column(q,r,x,n,m,ind)
%        [q,r,ind]=qr_add_column(q,r,x,n,m,ind)
%*********************************************************
%OBJETIVO: Se readapta la factorizacion QR de una matriz *
%          cuando se le anade una columna que ocupara el *
%          ultimo lugar en la matriz).                   *
%AUTORA: C. Pola, Universidad de Cantabria.              *
%*********************************************************
%
%LISTA DE LLAMADA
%DE ENTRADA
%  q   : contiene a la matriz "Q" inicial
%  r   : contiene a la matriz "R" inicial
%  x   : contiene a la columna que se desea anadir
%        No se usa si ind es distinto de cero.
%  n   : numero de filas de la matriz a factorizar
%  m   : numero de columnas de la matriz "R" inicial
%  ind : j si la nueva columna es el j-esimo vector de
%          la base canonica 
%       -j si la nueva columna es el j-esimo vector de
%          la base canonica cambiado de signo
%        0 si la nueva columna es un vector cualquiera
%DE SALIDA
%  q   :contiene a la nueva matriz "Q" 
%  r   :contiene a la nueva matriz "R" 
%  ind : -1 la nueva columna es linealmente dependiente de
%           las columnas que ya formaban la matriz
%         0 el proceso se ha realizado sin problemas
%

%
% Se calcula la columna de R correspondiente a la columna
% anadida
%
m1=m+1;nm=n-m;k=0;
if ind<0,   k=1;ind=-ind;end
if ind==0,
   if m>0,r(1:m,m1)=q(:,1:m)'*x;end
   w=q(:,m1:n)'*x;
else
   if m>0,r(1:m,m1)=q(ind,1:m)';end
   w=q(ind,m1:n)';
end
if k==1,
   if m>0,r(1:m,m1)=-r(1:m,m1);end
   w=-w;
end 
%
% Se averigua si la nueva columna es linealmente dependiente
%de las columnas que ya forman la matriz
%
rnorma=norm(w,2);eps1=eps^0.9;
if rnorma < eps1,
   ind=-1;
   return;
end
%
% Si la columna es linealmente independiente se procede a
%triangular la matriz
%
ind=0;
if m==n,
   r(m1,m1)=w(1);
   return;
else,
   if w(1)~=0, rnorma=sign(w(1))*rnorma; end
   w(1)=rnorma+w(1);
   s=sqrt(w(1)*rnorma);
   s=1.d0/s;
   w=s*w;
   for j=1:n,
       q(j,m1:n)=q(j,m1:n)-(q(j,m1:n)*w)*w';
   end
   r(m1,m1)=-rnorma;
end
return
