%EJEMPLO 2 SIMPLEX PL GENERAL:
clear all
format short e
clc

%Primero calculamos punto admisible con Fase 1:
%Llamada: [x,gamma, ind]=simplex_phase2(c,A,b,n,nI,nD,itemax,imp)

c=[-3/4 20 -1/2 6]';
A=[[1/4 -8 -1 9; 1/2 -12 -1/2 3; 0 0 1 0]' -eye(4)];
b=[0 0 1 zeros(1,4)]';
n=4;
nI=0;
nD=7;
itemax=1000;
imp=2;

%[x,gamma,ind]=simplex_phase1(A,b,n,nI,nD,itemax,imp);

%[x,fx,lambda,m,itrab,iter,ind]=simplex_phase2(c,A,b,n,nI,nD,x,itemax,imp);

[x,fx,lambda,m,itrab,iter,ind]=simplex_solver(c,A,b,n,nI,nD,[],itemax,imp);