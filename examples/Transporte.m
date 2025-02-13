%EJEMPLO TRANSPORTE SIMPLEX PL GENERAL:
clear all
format short e
clc

%Datos
n=12;
nI=6;
nD=22;
A=[[eye(3) eye(3) zeros(3,6);
     zeros(3,6) eye(3) eye(3)]' [ones(1,3) zeros(1,9);...
    zeros(1,3) ones(1,3) zeros(1,6);...
    zeros(1,6) ones(1,3) zeros(1,3);...
    zeros(1,9) ones(1,3);...
    eye(3) zeros(3) eye(3) zeros(3);...
    zeros(3) eye(3) zeros(3) eye(3)]' -eye(12)];
b=[4;5;4;3;3;4;5;8;7;4;6;3;8;3;9;3;zeros(12,1)].*10^3;
c=10^3.*[6;14;7;10;8;15;6;14;7;10;8;15];
x=[];
itemax=1000;
imp=2;


[x,fx,lambda,m,itrab,iter,ind]=simplex_solver(c,A,b,n,nI,nD,x,itemax,imp);