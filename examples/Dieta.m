%EJEMPLO DIETA SIMPLEX PL GENERAL:
clear all
format short e
clc

%Datosn=11;
n=11;
nD=21;
nI=0;
c=-[5; 3.50; 6; 5.50; 4.50; 3; 5.5; 5; 5.25; 3.75; 5.75];
A=zeros(10,11);
A(1,1:4)=0.01*[25,-75,25,25];
A(2,1:4)=0.01*[-20 -20 80 -20];
A(3,5:8)=0.01*[25 -75 25 25];
A(4,5:8)=0.01*[-25 -25 -25 75];
A(5,9:11)=0.01*[-75 25 25];
A(6,9:11)=0.01*[25 -75 25];
A(7:10,1:8)=[eye(4) eye(4)];
A(7:10,9:11)=[1 0 0;0 1 0; 0 0 0;0 0 1];
A=[A' -eye(11)];
b=[zeros(1,6) 1000 1000 750 800 zeros(n,1)']';
x=[];
itemax=1000;
imp=2;

[x,fx,lambda,m,itrab,iter,ind]=simplex_solver(c,A,b,n,nI,nD,x,itemax,imp);

%[x,gamma,ind,ind2]=simplex_phase1(A,b,n,nI,nD,itemax,imp);