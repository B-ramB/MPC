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
    a = 2;

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
xdest = [0 0 pi/3 0 0 0 pi/3 0]'; %waardes waar het naartoe moet bewegen.
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
[s,~,x] = lsim(KMIMO,r,t,x0);         % simulatie op de step response

PlotBeams(s,x,1,'State Feedback',Parameters,xdest) % (s, l1, l2, figure number, Controller)

%% MPC
MIMOdisc = c2d(MIMO,dt,'zoh'); %discretizeren van de plant

N = 20; %horizon input

A = MIMOdisc.A;
B = MIMOdisc.B;
C = MIMOdisc.C;
D = MIMOdisc.D;

Q = diag([1 1 1 1 1 1 1 1]); %Waardes van de diagonaal van de Q (allemaal dezelfde)
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

QH = zeros((N)*8); %Q zo opzetten dat je het kan vermenigvuldigen met de vector

% Constraints
F = [eye(8);
    -eye(8);
    1 0 1 0 1 0 1 0;
    -1 0 -1 0 -1 0 -1 0];
e = [pi;     % x1        upper
    pi;        % xdot1     upper
    pi;       % s1        upper
    pi;       % sdot1     upper
    pi;      % x2        upper
    pi;        % xdot2     upper
    pi;       % s2        upper
    pi;       % sdot2     upper
    pi;      % x1        lower
    pi;        % xdot1     lower
    pi;       % s1        lower
    pi;       % sdot1     lower
    pi;      % x2        lower
    pi;        % xdot2     lower
    pi;       % s2        lower
    pi;       % sdot2     lower
    2*pi;       % x1 + s1 + x2 + s2 < 
    2*pi];      % -x1 - s1 - x2 - s2 < 

% uiteindelijke waardes van x en u in vectoren zetten.
xref = zeros(8*N,1);
uref = zeros(2*N,1);
FLarge = zeros(size(F,1)*N,8*N);
eLarge = zeros(length(e)*N,1);
for j = 1:N
    xref(j*8-7:j*8,1) = xdest;
    uref(j*2-1:j*2,1) = udest;
    QH(j*8-7:j*8,j*8-7:j*8) = Q;
    FLarge(j*size(F,1)-(size(F,1)-1):j*size(F,1),j*size(F,2)-(size(F,2)-1):j*size(F,2)) = F;
    eLarge(j*length(e)-(length(e)-1):j*length(e),1)= e; 
end

[P,~,~] = idare(A,B,Q,R,[],[]);

QH(N*8-7:N*8,N*8-7:N*8) = P;


uopt = zeros(2,N); 

x = zeros(8,length(t)+1);
x(:,1) = x0;

% optimalisatie beginnen
for n = 1:length(t)
x0 = x(:,n);    
    
H = 0.5*(S'*QH*S + 2*R*eye(N*2));
h = x0'*T'*QH*S - xref'*QH*S - uref'*R*eye(N*2);    
    
cvx_begin quiet

variable u(2*N,1)

minimize(u'*H*u + h*u)
subject to
FLarge*(T*x0+S*u) <= eLarge;

cvx_end

x(:,n+1) = A*x(:,n) + B*u(1:2);

uopt(:,n) = u(1:2);
 
end

[s,~,x] = lsim(MIMOdisc,uopt,t,x(:,1));


PlotBeams(s,x,2,'MPC',Parameters,xdest)







%% funtions

function [] = PlotBeams(s,x,FigureNumber,Controller,Parameters,xdest)
% Plots the beams using the computed output angles. 
    % FigureNumber  : Which figure should it be plotted in
    % title         : What should be the title of the figure
    % y             : output sequence to be plotted (s)

    if length(s) < length(Parameters.t)    % extend s to the simulation time. 
        s = [s;
            s(end,1)*ones(length(Parameters.t)-length(s),1), s(end,2)*ones(length(Parameters.t)-length(s),1)];
    end
    
    for i = 1:length(Parameters.t)      % code voor de beam simulatie plot
        xm(i,1) = x(i,1)*10^(-3);       % from mm to m
        xm(i,5) = x(i,5)*10^(-3);       % from mm to m
        
        alpha1 = asin(xm(i,1)/(Parameters.l1/2));
        alpha2 = asin(xm(i,5)/(Parameters.l2/2));
        
        b11 = [cos(s(i,1))*Parameters.l1/2; sin(s(i,1))*Parameters.l1/2];                    % locatie van het uiteinde van beam 11
        b12 = b11 + [Parameters.l1/2*cos(s(i,1)+alpha1); Parameters.l1/2*sin(s(i,1)+alpha1)];

        b21 = b12 + [Parameters.l1/2*cos(s(i,1)+alpha1+s(i,2)); Parameters.l1/2*sin(s(i,1)+alpha1+s(i,2))];
        b22 = b21 + [Parameters.l1/2*cos(s(i,1)+alpha1+s(i,2)+alpha2); Parameters.l1/2*sin(s(i,1)+alpha1+s(i,2)+alpha2)];
        
        beam11 = [zeros(2,1) b11];
        beam12 = [b11 b12];
        beam21 = [b12 b21];
        beam22 = [b21 b22];
        
        b1 = [cos(s(i,1))*Parameters.l1; sin(s(i,1))*Parameters.l1];                    % locatie van het uiteinde van beam 1
        b2 = b1 + [cos(s(i,2)+s(i,1))*Parameters.l2; sin(s(i,2)+s(i,1))*Parameters.l2]; % locatie van het uiteinde van beam 2

        beam1 = [zeros(2,1) b1];    % vector gelijk aan beam 1
        beam2 = [b1 b2];            % vector gelijk aan beam 2
        
        
        figure(FigureNumber)        % plot voor de beams simulatie
        sgtitle(['Controller: ',Controller])
        subplot(1,2,1)
            plot(0,0,'ko')
            hold on
                % plot(b1(1),b1(2),'bo',b2(1),b2(2),'ro') % rondjess
                % plot(beam1(1,:),beam1(2,:),'b')     % beam1
                % plot(beam2(1,:),beam2(2,:),'r')     % beam2
                
                plot(0,0,'ko',b12(1),b12(2),'bo',b22(1),b22(2),'ro')   % rondjess
                plot(beam11(1,:),beam11(2,:),'b')   % b11
                plot(beam12(1,:),beam12(2,:),'b')   % b12
                plot(beam21(1,:),beam21(2,:),'r')   % b21
                plot(beam22(1,:),beam22(2,:),'r')   % b22
            hold off
            axis([-5 5 -5 5])
            title 'Two arms' 
            axis square
        subplot(1,2,2)
            plot(Parameters.t(1:i),s(1:i,1),'b',Parameters.t(1:i),s(1:i,2),'r')
            hold on
                plot(Parameters.t(1:i),x(1:i,1),'color',[0, 0.4470, 0.7410])
                plot(Parameters.t(1:i),x(1:i,5),'color',[0.8500, 0.3250, 0.0980]);
            hold off
            yline(xdest(3),'b--');
            yline(xdest(7),'r--');
            axis([0 Parameters.T -1 2])
            title '$s_1$ \& $s_2$'
            legend('s_1 [rad]','s_2 [rad]','x_1 [mm]','x_2 [mm]');    
        drawnow;
        pause(Parameters.dt) 
    end

end
