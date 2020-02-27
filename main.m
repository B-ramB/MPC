clc;
clear all;
close all;

M1 = 2.75; %kg
k1 = 511.395; %N/m 
m1 = 5.067; %kg
c1 = 2.5593; %Ns/m
M2 = 2.75; %kg
k2 = 511.395; %N/m 
m2 = 5.067; %kg
c2 = 2.5593; %Ndm

A= [ 0 1 0 0 0 0 0 0
    -k1*(1/M1 + 1/m1) -c1/m1 0 0 0 0 0 0
    0 0 0 1 0 0 0 0
    k1/M1 0 0 0 0 0 0 0
    0 0 0 0 0 1 0 0
    0 0 0 0 -k2*(1/M2 + 1/m2) -c2/m2 0 0
    0 0 0 0 0 0 0 1
    k1/m1 c1/m1 0 0 k2/M2 0 0 0];

B= [0 -1/M1 0 1/M1 0 0 0 0
    0 -1/m1 0 0 0 -1/M2 0 (1/m1 + 1/M2)]';

C= [0 0 1 0 0 0 0 0
    0 0 0 0 0 0 1 0];
D= 0;

MIMO=ss(A,B,C,D);
step(MIMO);

MIMOdisc = c2d(MIMO,0.020,'zoh');

figure(2)
step(MIMOdisc,1);

