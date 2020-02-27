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

K = [-51.59 32.61 314.44 96.21 284.79 9.68 -33.57 15.15
   163.16 37.98 33.57 33.98 -433.18 -7.09 314.441 80.18];

MIMO=ss(A,B,C,D);
KMIMO=ss(A-B*K,B,C,D);
lsim(KMIMO,zeros(2,501),0:0.01:5,ones(8,1));

% MIMOdisc = c2d(MIMO,0.020,'zoh');
% 
% figure(2)he
% step(MIMOdisc,1);

