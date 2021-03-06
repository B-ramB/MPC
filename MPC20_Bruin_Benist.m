clc; close all; clear all;
%%%%%%%%%%%%%%%%%
%{
Code written by:
    Bram Benist  (4281926)
    Pim de Bruin (4545702)

For the course 'Model Prediceive control' SC42125 at the TU Delft 
06/04/2020

The code requirec the CVX optimisation software to run. 
%} 
%%%%%%%%%%%%%%%%%%



set(0,'defaultTextInterpreter','latex'); % Set latex as default text interpreter
%% Settings
Plots = true;           % true/false

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


MIMO=ss(A,B,C,D);               % MIMO plantco
KMIMO=ss(A-B*K,B,C,D);          % MIMO plant met state feedback controller
KMIMO = KMIMO /dcgain(KMIMO);   % Controlled MIMO plant gedeeld door zijn dc gain. 
disp('Simulating original state feedback controller')
[s,~,x] = lsim(KMIMO,r,t,x0);         % simulatie op de step response
if Plots == true
    PlotBeams(s,x,1,'State Feedback',Parameters,xdest) % (s, l1, l2, figure number, Controller)
end

for n = 1:length(s)
    if abs(s(n,1)-xdest(3))< 0.01*xdest(3) && abs(s(n,1)-xdest(7)) < 0.01*xdest(3)
     % disp(n)
    end
end

%% State Measurement MPC


MIMOdisc = c2d(MIMO,dt,'zoh'); %discretizeren van de plant

N = 13; %horizon input

A = MIMOdisc.A;
B = MIMOdisc.B;
C = MIMOdisc.C;
D = MIMOdisc.D;

Q = diag([20 50 500 50 20 50 500 50]); %Waardes van de diagonaal van de Q (allemaal dezelfde)
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
     1 0 1 0 1 0 1 0;
    -eye(8);
    -1 0 -1 0 -1 0 -1 0];

m = [pi;      % x1        lower
    1000;        % xdot1     lower
    pi;       % s1        lower
    2*pi;       % sdot1     lower
    pi;      % x2        lower
    1000;        % xdot2     lower
    pi;       % s2        lower
    2*pi;         % sdot2     lower
    2*pi];       % x1 + s1 + x2 + s2 < 

M = [pi;     % x1        upper
    1000;        % xdot1     upper
    pi;       % s1        upper
    2*pi;       % sdot1     upper
    pi;      % x2        upper
    1000;        % xdot2     upper
    pi;       % s2        upper
    2*pi;       % sdot2     upper   
    2*pi];      % -x1 - s1 - x2 - s2 < 

e = [M; m];

Fu = [eye(2); -eye(2)];
eu = [100; 100; 100; 100];

% uiteindelijke waardes van x en u in vectoren zetten.
xref = zeros(8*N,1);
uref = zeros(2*N,1);
FLarge = zeros(size(F,1)*N,8*N);
eLarge = zeros(length(e)*N,1);
FuLarge = zeros(size(Fu,1)*N,2*N);
euLarge = zeros(length(eu)*N,1);
for j = 1:N
    xref(j*8-7:j*8,1) = xdest;
    uref(j*2-1:j*2,1) = udest;
    QH(j*8-7:j*8,j*8-7:j*8) = Q;
    FLarge(j*size(F,1)-(size(F,1)-1):j*size(F,1),j*size(F,2)-(size(F,2)-1):j*size(F,2)) = F;
    eLarge(j*length(e)-(length(e)-1):j*length(e),1)= e; 
    FuLarge(j*size(Fu,1)-(size(Fu,1)-1):j*size(Fu,1),j*size(Fu,2)-(size(Fu,2)-1):j*size(Fu,2)) = Fu;
    euLarge(j*length(eu)-(length(eu)-1):j*length(eu),1)= eu; 
end

[P,K,~] = idare(A,B,Q,R,[],[]);

QH(N*8-7:N*8,N*8-7:N*8) = P;

uopt = zeros(2,N); 

x = zeros(8,length(t)+1);
x(:,1) = x0;

% Main loop

%% Asymptotic stability
disp('Checking asymptotic stability')
T_s = 0.1;
[Sas,D] = eig(P);         % diagonalise
Sas = normc(Sas);           % normalise S
D = sqrt(2*T_s*inv(D));        % get magnitudes of boundary xf

Points = [Sas*D -Sas*D]; 


Cas = combnk(1:16,8);
Vertices = zeros(8,size(Cas,1));
for i = 1:size(Cas,1)
    Vertices(:,i) = sum(Points(:,Cas(i,:)),2);
end
Upper = Vertices <= M(1:end-1);
Lower = Vertices >= -m(1:end-1);

uas = K*Vertices;

Check = sum([-100;-100] < uas < [100; 100],2);

if sum(sum(Upper,2))+sum(sum(Lower,2)) < 2*size(Vertices,1)*size(Vertices,2) || sum(Check) < 2*size(Vertices,2)
    disp('xf not in x or u not in U')
elseif sum(sum(Upper,2))+sum(sum(Lower,2)) == 2*size(Vertices,1)*size(Vertices,2) && sum(Check) == 2*size(Vertices,2)
    disp('xf in x! and u in U!')
end


%%
% optimalisatie beginnen

Vfplus = zeros(1,length(t));
Vf = zeros(1,length(t));
l = zeros(1,length(t));

disp('Simulating state measurement MPC')

for n = 1:length(t)
x0 = x(:,n);    
    
H = 0.5*(S'*QH*S + 2*R*eye(N*2));
h = x0'*T'*QH*S - xref'*QH*S - uref'*R*eye(N*2);    
    
cvx_begin quiet

variable u(2*N,1)

minimize(u'*H*u + h*u)
subject to
FLarge*(T*x0+S*u) <= eLarge;
FuLarge*u <= euLarge;
0.5*((T(N*8-7:8*N,:)*x0 + S(N*8-7:N*8,:)*u) - xdest)'*P*((T(N*8-7:8*N,:)*x0 + S(N*8-7:N*8,:)*u) - xdest) <= T_s;

cvx_end

x(:,n+1) = A*x(:,n) + B*u(1:2);
y(:,n) = C*x(:,n);

uopt(:,n) = u(1:2);
 



Vfplus(n) = (x(:,n+1)-xdest)'*P*(x(:,n+1)-xdest);
Vf(n) = (x(:,n)-xdest)'*P*(x(:,n)-xdest);
l(n) = 1/2*(x(:,n)-xdest)'*Q*(x(:,n) - xdest) + uopt(:,n)'*R*uopt(:,n);
end


% figure 
%     hold on
%     plot(t,Vfplus - Vf,'b')
%     plot(t,-l,'r')
%     xlabel 'time (s)'
%     ylabel 'Cost'
%     legend('$V_f(f(x,u) - V_f(x)$','-L','Interpreter','latex')



if Plots == true
    PlotBeams(y',x',2,'state-measurement MPC',Parameters,xdest)
end

for n = 1:length(y)
    if abs(y(1,n)-xdest(3))< 0.01*xdest(3) && abs(y(2,n)-xdest(7)) < 0.01*xdest(7)
     % disp(n)
    end
end



%% output measurement met disturbance matrices en observer

T = 10;          % simulatie tijd
Parameters.T = T;
dt = 0.1;       % time step
t = 0:dt:T;
Parameters.t = t;

Q = diag([20 50 500 40 20 50 500 40]);

Q2 = blkdiag( Q, zeros(2)); %Waardes van de diagonaal van de Q 


% Disturbance
NoiseLevel = 80;
Dist = 1*[0.001;            % constant disturbance on [s1; s2]
        0.001];
Bd = [1, 0
      0, 0
      0, 0
      0, 0
      0, 1
      0, 0
      0, 0
      0, 0];
Cd = zeros(2,2);

%Checks
OB = obsv(A,C);
rank(OB');
nnd = rank([(eye(8)-A), -Bd
              C,         Cd]);
C2 = [C Cd];
A2 = [A, Bd;
      zeros(2,8), eye(2)];
B2 = [B; zeros(2,2)];

L = place(A2',C2',[0.7, 0.75, 0.8, 0.85, 0.9, 0.7, 0.75, 0.8, 0.85, 0.9]'); % Place observer poles

% Opmaken van de T en S matrix van de vergelijking x = Tx0 + Su voor de
% disturbance rejection
T2 = zeros(10*N,10);
T2(1:10,1:10) = eye(10);
S2 = zeros(10*N,2*N);
for i = 2:N
 T2((10*i)-9:10*i,:) = A2^(i-1);
 for k = 1:i-1
 S2((10*i)-9:10*i,2*k-1:2*k) = A2^(i-k-1)*B2;
 end
end

% matrices voor the estimators and true values
ytrue = zeros(2,length(t));
xhat = zeros(10,length(t)); % LET OP var xhat bevat zowel xhat als dhat
xtrue = zeros(10,length(t));
y = zeros(2,length(t));

% Initial values
xtrue(:,1) = [zeros(8,1);       % Initialise with the disturbance
               Dist];
yr = C*xdest;                   % zelfde destination als vorige mpc controller


equal = [eye(8)-A, -B           % equality matrix
         C, zeros(2,2)];

W1 = eye(8);                    % weighting matrices voor de OTS
W2 = eye(2);

% Cost function and constraints
F2 = [F zeros(18,2)];           
F2Large = zeros(size(F,1)*N,10*N);

Q2H = zeros(N*10,N*10);
for j = 1:N     
    Q2H(j*10-9:j*10,j*10-9:j*10) = Q2;
    F2Large(j*size(F2,1)-(size(F2,1)-1):j*size(F2,1),j*size(F2,2)-(size(F2,2)-1):j*size(F2,2)) = F2;
end
Q2H(N*10-9:N*10,N*10-9:N*10) = [P, zeros(8,2); zeros(2,8), zeros(2,2)]; % add terminal cost matrix

% Initialisaties
xrtot = zeros(8,length(t));
urtot = zeros(2,length(t));
x2ref = zeros(8*N,1);
u2ref = zeros(2*N,1);

disp(['Simulating Output-measurement MPC',newline...
    'Noiselevel = ',num2str(NoiseLevel), newline...
    'Disturbance = [',num2str(Dist'),']']);

% Main loop
for n = 1:length(t)
    
% OTS
dhat = xhat(9:10,n);   %dhat aanroepen van xhat 
cvx_begin quiet

variable ur(2,1)
variable xr(8,1)
minimize(xr'*W1*xr + ur'*W2*ur)                % choose best target
subject to
equal*[xr; ur] == [Bd*dhat; (yr-Cd*dhat)];

cvx_end

for j = 1:N                                     % put in appropriate matrices
    x2ref(j*10-9:j*10,1) = [xr; 0; 0];
    u2ref(j*2-1:j*2,1) = ur;
end

% MPC
x0 = xhat(1:10,n);    %xhat als 0 zetten 
H = 0.5*(S2'*Q2H*S2 + 2*R*eye(N*2));
h = x0'*T2'*Q2H*S2 - x2ref'*Q2H*S2 - u2ref'*R*eye(N*2);    
    
cvx_begin quiet

variable u(2*N,1)
minimize(u'*H*u + h*u)              % optimization zoals vorige
subject to
F2Large*(T2*x0+S2*u) <= eLarge;     
FuLarge*u <= euLarge;
0.5*((T2(N*10-9:10*N,:)*x0 + S2(N*10-9:N*10,:)*u) - [xr; 0; 0])'*[P, zeros(8,2); zeros(2,8), zeros(2,2)]*((T2(N*10-9:10*N,:)*x0 + S2(N*10-9:N*10,:)*u)- [xr; 0; 0]) <= T_s;
cvx_end


% Simulate
xhat(:,n+1) = A2*xhat(:,n) + B2*u(1:2,1) + L'*(y(:,n) - C2*xhat(:,n));      % Observer
xtrue(:,n+1) = A2*xtrue(:,n) + B2*u(1:2,1) ;                                 % True system
ytrue(:,n+1) = C2*xtrue(:,n+1);
y(:,n+1) = awgn(ytrue(:,n+1),NoiseLevel);                                           % Measurement noise

% Save
urtot(:,n) = ur;
xrtot(:,n) = xr;
uopt(:,n) = u(1:2);
end

% look at the observer state convergence:
% err = vecnorm(xhat - xtrue);
% figure
%     hold on
%     plot(t,xhat(9:10,1:end-1)-xtrue(9:10,1:end-1))
%     legend('$\hat{d}_1 - dist_1$','$\hat{d_2} - dist_2$','Interpreter','latex');


if Plots == true
    PlotBeams(y',xtrue(1:8,:)',3,'Output-measurement MPC',Parameters,xdest)
end

for n = 1:length(y)
    if abs(ytrue(1,n)-xdest(3))< 0.01*xdest(3) && abs(ytrue(2,n)-xdest(7)) < 0.01*xdest(7)
     % disp(n)
    end
end

disp('Done')

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
            title '$$s_1$$ \& $$s_2$$'
            legend('s_1 [rad]','s_2 [rad]','x_1 [mm]','x_2 [mm]');    
        drawnow;
        % pause(Parameters.dt) 
        
        
        
    end

%     figure(FigureNumber+10) 
%                 plot(Parameters.t(1:i),s(1:i,1),'b',Parameters.t(1:i),s(1:i,2),'r')
%             hold on
%                 plot(Parameters.t(1:i),x(1:i,1),'color',[0, 0.4470, 0.7410])
%                 plot(Parameters.t(1:i),x(1:i,5),'color',[0.8500, 0.3250, 0.0980]);
%             hold off
%             yline(xdest(3),'b--');
%             yline(xdest(7),'r--');
%             axis([0 Parameters.T -1 2])
%             title '$$s_1$$ \& $$s_2$$'
%             legend('s_1 [rad]','s_2 [rad]','x_1 [mm]','x_2 [mm]'); 
    grid on
end
