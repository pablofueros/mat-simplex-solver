%EJEMPLO 1 SIMPLEX PL GENERAL:
clear all
format short e
clc

%Datos del problema
c=[1 2 3 4]';
A=[[1 1 1 1; 1 0 1 -3]' -eye(4)];
b=[1 1/2 zeros(1,4)]';
n=4;
nI=2;
nD=4;
itemax=1000;
imp=2;

%Primero calculamos punto admisible con Fase 1:
%Llamada: [x,gamma, ind]=simplex_phase2(A,b,n,nI,nD,itemax,imp)

%[x,gamma,ind]=simplex_phase1(A,b,n,nI,nD,itemax,imp);

%Obtenemos punto admisible:
% x =
% 
%             0
%             0
%    8.7500e-01
%    1.2500e-01

%Calculamos solucion con FaseII
%[x,fx,lambda,m,itrab,iter,ind]=simplex_phase2(c,A,b,n,nI,nD,x,itemax,imp)

%[x,fx,lambda,m,itrab,iter,ind]=simplex_phase2(c,A,b,n,nI,nD,x,itemax,imp);

[x,fx,lambda,m,itrab,iter,ind]=simplex_solver(c,A,b,n,nI,nD,[],itemax,imp);