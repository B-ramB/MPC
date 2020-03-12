clc;
clear all;
close all;

set(0,'defaultTextInterpreter','latex'); % Set latex as default text interpreter

% variables
M1 = 2.75;      % kg
k1 = 511.395;   % N/m 
m1 = 5.067;     % kg
c1 = 2.5593;    % Ns/m
M2 = 2.75;      % kg
k2 = 511.395;   % N/m 
m2 = 5.067;     % kg
c2 = 2.5593;    % Ndm
l1 = 2;         % m
l2 = 2;         % m
    % put variables in a struct
    Parameters.M1 = M1; % kg
    Parameters.k1 = k1; % N/m 
    Parameters.m1 = m1; % kg
    Parameters.c1 = c1; % Ns/m
    Parameters.M2 = M2; % kg
    Parameters.k2 = k2; % N/m 
    Parameters.m2 = m2; % kg
    Parameters.c2 = c2; % Ndm
    Parameters.l1 = l1; % m
    Parameters.l2 = l2; % m

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

% Simulation settings:
T = 5;          % simulatie tijd
dt = 0.1;       % time step
t = 0:dt:T;     % time vector
xdest = [0 0 pi/2 0 0 0 pi/2 0]'; %waardes waar het naartoe moet bewegen.
udest = [0 0]'; %uiteindelijke waardes van u
    % Also put it in the parameters struct
    Parameters.T = T;
    Parameters.dt = dt;
    Parameters.t = t;

%% Controller volgens de originele paper (State feedback)

K = [-51.59 32.61 314.44 96.21 284.79 9.68 -33.57 15.15
   163.16 37.98 33.57 33.98 -433.18 -7.09 314.441 80.18];

x0 = [0 0 0 0 0 0 0 0]';            % x0
r = [xdest(3)*ones(1,(length(t)))
     xdest(7)*ones(1,(length(t)))];  % input over het aantal simulatie stappen 


MIMO=ss(A,B,C,D);               % MIMO plant
KMIMO=ss(A-B*K,B,C,D);          % MIMO plant met state feedback controller
KMIMO = KMIMO /dcgain(KMIMO);   % Controlled MIMO plant gedeeld door zijn dc gain. 
s = lsim(KMIMO,r,t,x0);         % simulatie op de step response

 PlotBeams(s,1,'State Feedback',Parameters) % (s, l1, l2, figure number, Controller)

%% MPC
MIMOdisc = c2d(MIMO,dt,'zoh'); %discretizeren van de plant

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

% uiteindelijke waardes van x en u in vectoren zetten.
xref = zeros(8*N,1);
uref = zeros(2*N,1);
for j = 1:N
    xref(j*8-7:j*8,1) = xdest;
    uref(j*2-1:j*2,1) = udest;
end


% waardes van de cost function zo opzetten dat je u'*H*u + h*u kan zeggen
H = 0.5*(S'*QH*S + 2*R*eye(N*2));
h = x0'*T'*QH*S - xref'*QH*S - uref'*R*eye(N*2);

% optimalisatie beginnen
cvx_begin

variable u(2*N,1)

minimize(u'*H*u + h*u)
       
cvx_end

% optimale u waardes in een correcte vector zetten
uopt = zeros(2,N); 
for j = 1:N
 uopt(:,j) = u(2*j-1:2*j);
end

s = lsim(MIMOdisc,uopt,0:dt:dt*(N-1),x0); %simuleren van y met optimale inputs

PlotBeams(s,2,'MPC',Parameters)







%% funtions

function [] = PlotBeams(s,FigureNumber,Controller,Parameters)
% Plots the beams using the computed output angles. 
    % FigureNumber  : Which figure should it be plotted in
    % title         : What should be the title of the figure
    % y             : output sequence to be plotted (s)

    if length(s) < length(Parameters.t)    % extend s to the simulation time. 
        s = [s;
            s(end,1)*ones(length(Parameters.t)-length(s),1), s(end,2)*ones(length(Parameters.t)-length(s),1)];
    end
    
    for i = 1:length(Parameters.t)  % code voor de beam simulatie plot

        b1 = [cos(s(i,1))*Parameters.l1; sin(s(i,1))*Parameters.l1];                    % locatie van het uiteinde van beam 1
        b2 = b1 + [cos(s(i,2)+s(i,1))*Parameters.l2; sin(s(i,2)+s(i,1))*Parameters.l2]; % locatie van het uiteinde van beam 2

        beam1 = [zeros(2,1) b1];    % vector gelijk aan beam 1
        beam2 = [b1 b2];            % vector gelijk aan beam 2

        figure(FigureNumber)        % plot voor de beams simulatie
        sgtitle(['Controller: ',Controller])
        subplot(1,2,1)
            plot(0,0,'ko',b1(1),b1(2),'bo',beam1(1,:),beam1(2,:),'b',b2(1),b2(2),'ro',beam2(1,:),beam2(2,:),'r')
            axis([-5 5 -5 5])
            title 'Two arms'  
        subplot(1,2,2)
            plot(Parameters.t(1:i),s(1:i,1),'b',Parameters.t(1:i),s(1:i,2),'r')
            axis([0 Parameters.T 0 2])
            legend('s1','s2');
            title '$s_1$ \& $s_2$'
            
            
        pause(Parameters.dt) 
    end

end
