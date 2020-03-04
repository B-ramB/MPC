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
l1 = 2;
l2 = 2;

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

t= 5;

x0 = [0 0 0 0 0 0 0 0]';
u = [500*ones(1,(t*10+1))
     100*ones(1,(t*10+1))];


MIMO=ss(A,B,C,D);
KMIMO=ss(A-B*K,B,C,D);
s = lsim(KMIMO,u,0:0.1:t,x0);

figure(1)
plot(0:0.1:t,s(:,1),'b',0:0.1:t,s(:,2),'r')
legend('s1','s2');



for i = 1:length(s)

    b1 = [cos(s(i,1))*l1; sin(s(i,1))*l1];
    b2 = b1 + [cos(s(i,2))*l2; sin(s(i,2))*l2];
    
    beam1 = [zeros(2,1) b1];
    beam2 = [b1 b2];        
    
    figure(2)
    plot(0,0,'ko',b1(1),b1(2),'bo',beam1(1,:),beam1(2,:),'b',b2(1),b2(2),'ro',beam2(1,:),beam2(2,:),'r')
    axis([-5 5 -5 5])
    
    pause(0.1)
    
    
end

%%
MIMOdisc = c2d(MIMO,0.1,'zoh');

N = 10;

A = MIMOdisc.A;
B = MIMOdisc.B;
C = MIMOdisc.C;
D = MIMOdisc.D;

Q = 1;
R = 1;



T = zeros(8*N,8);
T(1:8,1:8) = eye(8);
S = zeros(8*N,2*N);
for i = 2:N+1
 T((8*i)-7:8*i,:) = A^(i-1);
 for k = 1:i-1
 S((8*i)-7:8*i,2*k-1:2*k) = A^(i-k-1)*B;
 end
end


QH = Q*eye((N+1)*8);


H = 0.5*(S'*QH*S + R*eye(N*2));
h = x0'*T'*QH*S;


cvx_begin

variable u(N*2,1)

minimize(u'*H*u + h*u)
subject to
    ones(N,1)*[0 0 1 0 0 0 0 0]*(T*x0 + S*u) <= ones(N,1)*pi/2;
    
cvx_end



