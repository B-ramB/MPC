clc;
clear all;
close all;

%variables
M1 = 2.75; %kg
k1 = 511.395; %N/m 
m1 = 5.067; %kg
c1 = 2.5593; %Ns/m
M2 = 2.75; %kg
k2 = 511.395; %N/m 
m2 = 5.067; %kg
c2 = 2.5593; %Ndm
l1 = 2; %m
l2 = 2; %m


% Continuous time system
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

%Controller volgens de originele paper

K = [-51.59 32.61 314.44 96.21 284.79 9.68 -33.57 15.15
   163.16 37.98 33.57 33.98 -433.18 -7.09 314.441 80.18];

t= 5; %simulatie tijd

x0 = [0 0 0 0 0 0 0 0]'; %x0
u = [500*ones(1,(t*10+1))
     500*ones(1,(t*10+1))]; % input over het aantal simulatie stappen 


MIMO=ss(A,B,C,D); %MIMO plant
KMIMO=ss(A-B*K,B,C,D); % MIMO plant met controller
s = lsim(KMIMO,u,0:0.1:t,x0); % simulatie op de step response

figure(1) % plot van de twee outputs
plot(0:0.1:t,s(:,1),'b',0:0.1:t,s(:,2),'r')
legend('s1','s2');



for i = 1:length(s) %code voor de beam simulatie plot

    b1 = [cos(s(i,1))*l1; sin(s(i,1))*l1]; %locatie van het uiteinde van beam 1
    b2 = b1 + [cos(s(i,2)+s(i,1))*l2; sin(s(i,2)+s(i,1))*l2]; %locatie van het uiteinde van beam 2
    
    beam1 = [zeros(2,1) b1]; %totale beam1 als lijn
    beam2 = [b1 b2];        % totale beam2 als lijn
    
    figure(2) %de simulatie plot
    plot(0,0,'ko',b1(1),b1(2),'bo',beam1(1,:),beam1(2,:),'b',b2(1),b2(2),'ro',beam2(1,:),beam2(2,:),'r')
    axis([-5 5 -5 5])
    
    pause(0.1)
    
    
end

%%
MIMOdisc = c2d(MIMO,0.1,'zoh'); %discretizeren van de plant

N = 100; %horizon input

A = MIMOdisc.A;
B = MIMOdisc.B;
C = MIMOdisc.C;
D = MIMOdisc.D;

Q = 1; %Waardes van de diagonaal van de Q (allemaal dezelfde)
R = 0; %waarde van R

x0 = [0 0 0 0 0 0 0 0]';

%Opmaken van de T en S matrix van de vergelijking x = Tx0 + Su
T = zeros(8*N,8);
T(1:8,1:8) = eye(8);
S = zeros(8*N,2*N);
for i = 2:N
 T((8*i)-7:8*i,:) = A^(i-1);
 for k = 1:i-1
 S((8*i)-7:8*i,2*k-1:2*k) = A^(i-k-1)*B;
 end
end


QH = Q*eye((N)*8); %Q zo opzetten dat je het kan vermenigvuldigen met de vector

xdest = [0 0 pi 0 0 0 pi 0]'; %waardes waar het naartoe moet bewegen.
udest = [0 0]'; %uiteindelijke waardes van u

% uiteindelijke waardes van x en u in vectoren zetten.
xref = zeros(8*N,1);
uref = zeros(2*N,1);
for l = 1:N
    xref(l*8-7:l*8,1) = xdest;
    uref(l*2-1:l*2,1) = udest;
end


%waardes van de cost function zo opzetten dat je u'*H*u + h*u kan zeggen
H = 0.5*(S'*QH*S + 2*R*eye(N*2));
h = x0'*T'*QH*S - xref'*QH*S - uref'*R*eye(N*2);

% optimalisatie beginnen
cvx_begin

variable u(2*N,1)

minimize(u'*H*u + h*u)
       
cvx_end

%optimale u waardes in een correcte vector zetten
uopt = zeros(2,N); 

for j = 1:N
 uopt(:,j) = u(2*j-1:2*j);
end

s = lsim(MIMOdisc,uopt,0:0.1:0.1*(N-1),x0); %simuleren van van y waardes

for i = 1:length(s) %code voor de beam simulatie plot

    b1 = [cos(s(i,1))*l1; sin(s(i,1))*l1]; %locatie van het uiteinde van beam 1
    b2 = b1 + [cos(s(i,2)+s(i,1))*l2; sin(s(i,2)+s(i,1))*l2]; %locatie van het uiteinde van beam 2
    
    beam1 = [zeros(2,1) b1]; %vector gelijk aan beam 1
    beam2 = [b1 b2];         %vector gelijk aan beam 2
    
    figure(2) %plot voor de beams simulatie
    plot(0,0,'ko',b1(1),b1(2),'bo',beam1(1,:),beam1(2,:),'b',b2(1),b2(2),'ro',beam2(1,:),beam2(2,:),'r')
    axis([-5 5 -5 5])
    
    pause(0.1)
    
    
end
