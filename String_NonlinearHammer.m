%-----------------------------------------------------
% PBMMI Assignment 4 BTB2
% ----- Tension modulation nonlinearity and coupling with ideal hammer
%
% Author: B200247 
% Date: 18/3/23
%
% (took me around 14s to run, and sounds 30s following the order:
% nonlinear, linear and 4 sets of couplings with hammer, thanks again for
% your time and patience!)  
%
% (Besides directly run, may also try turn on 'plot_on' in line 52)
% 
% This script includes 2 parts:
%
% A) > Tension modulation nonlinearity with modification from scheme 18 in
%      the reference - kept stiffness, frequency independent and
%      independent loss and tension modulation nonlinearity while excluded
%      fretboard and fingerstopping, parameters are also made similar with 
%      those in the paper while rho there equals rho*A in this script  
%    > Plotted a comparison between nonlinear and linear string output, 
%      shows similar change with figure 8.10 and 8.12 in NSS. By turning on
%      'plot_on' where you can also see the difference
%    > long elapsed time is mostly caused by the updates of A and C in loop
%      for non-linearity, Tf here is smaller than T60 for this trade off
%      issue
% B) > Coupled the linear stiff string with a very basic hammer (the one in
%      Chapter 4!) which, therefore doesn't really exist, but I had a lot
%      of fun with it. It had no width and only acts with one selected point
%    > plotted selected outputs and force/mass outputs uder 4 different
%      settings - may need to maximize the figure window for better view,
%      waveforms are all different  
%    > function starting from line 293
% 
% Also tried with another way of defining Dxx and Dxxxx (line 116) without
% mannually changing boundary values, as introduced in reference and it
% carries the theoretical meaning better  
%
% reference: Numerical Modeling of String/Barrier Collisions: The Fretboard. Bilbao, Stefan & Torin, Alberto. (2014). 
%-----------------------------------------------------

clc
clear all
close all


%% initial parameters

% flags -----------------------------------------------------
itype = 1;                  % type of input: 1: pluck, 2: strike
         
% in-loop plotting on (1) or off (0)
plot_on = 1;                % for nonlinearity
ploth1 = 0;                 % for hammer1
ploth2 = 0;                 % for hammer2
ploth3 = 0;                 % for hammer3
ploth4 = 0;                 % for hammer4

% parameters ------------------------------------------------
r = 0.00043;                % string radius (m)
E = 2e11;                   % Young's modulus (Pa) for steel
rho = 9000;                 % density (kg/m^3)
L = 0.65;                   % length (m)
loss = [100, 10; 1000, 8];  % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
SR = 44100;                 % sample rate (Hz)
Tf = 5;                     % duration of simulation (s)

%/////////////////////////////////////////////////////////////////////////%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Non-linearity %%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% more parameters
T = 60;                     % tension (N)
xi = 0.52;                  % coordinate of excitation (normalised, 0-1)
famp = 1.5;                 % peak amplitude of excitation (N)
dur = 0.001;                % duration of excitation (s)
exc_st = 0.001;             % start time of excitation (s)

xo = 0.1;                   % coordinate of output (normalised, 0-1)


%% derived parameters

A = pi*r^2;                 % string cross-sectional area
I = 0.25*pi*r^4;            % string moment of inertia

c = sqrt(T/(rho*A));        % wave speed
K = sqrt(E*I/(rho*A));      % stiffness constant 

k = 1/SR;                   % time step
Nf = floor(SR*Tf);          % number of time steps

%loss parameters
zeta1 = (-c^2+sqrt(c^4+4*K^2*(2*pi*loss(1,1))^2))/(2*K^2);
zeta2 = (-c^2+sqrt(c^4+4*K^2*(2*pi*loss(2,1))^2))/(2*K^2);
sig = 6*log(10)*(-zeta2/loss(1,2)+zeta1/loss(2,2))/(zeta1-zeta2);
sig1 = 6*log(10)*(1/loss(1,2)-1/loss(2,2))/(zeta1-zeta2);

%% grid spacing
hmin = sqrt(0.5*(c^2*k^2+4*sig1*k+sqrt((c^2*k^2+4*sig1*k)^2+16*K^2*k^2)));
N = floor(L/hmin);
h = L/N;


%% input
t = [0:Nf-1]'/SR;
ind = sign(max(-(t-exc_st).*(t-exc_st-dur),0));
f = 0.5*famp*ind.*(1-cos(itype*pi*(t-exc_st)/dur));


%% matrices for update

% input/output position in grid point index
xi_m = floor(N*xi);
xo_m = floor(N*xo);     

% matrices of spatial derivatives
e = ones(N,1);
Dx = spdiags([-1*e e], -1:0,N,N-1)/h;
Dxx = -Dx'*Dx;
Dxxxx = Dxx*Dxx;

% matrices for state update
% for nonlinearity
B = 2*speye(N-1)+c^2*k^2*Dxx-K^2*k^2*Dxxxx+2*sig1*k*Dxx;
% for linear comparison
B0 = 1/(1+sig*k)*(2*speye(N-1)+c^2*k^2*Dxx-K^2*k^2*Dxxxx+2*sig1*k*Dxx);
C0 = 1/(1+sig*k)*((-1+sig*k)*speye(N-1)-2*sig1*k*Dxx);

% input force vector
J = sparse(N-1,1);
J(xi_m) = 1/(1+sig*k)*k^2/(rho*A*h);
   
% output vector
cout = sparse(1,N-1);              
cout(xo_m) = 1;


%% initialise scheme variables

% for nonlinearity
u2 = zeros(N-1,1);              % state
u1 = u2;                        % state
u = u2;                         % state

y = zeros(Nf,1);                % output

% for linear comparison
u20 = zeros(N-1,1);             % state
u10 = u20;                      % state
u0 = u20;                       % state

y0 = zeros(Nf,1);               % output

xax = (1:N-1)'*h;               % x-axis for plotting

% Initialize the VideoWriter object
videoObj = VideoWriter('output_video.mp4', 'MPEG-4');
videoObj.FrameRate = 10; % Set the frame rate, e.g., 10 frames per second
open(videoObj);

%% main loop

for n=1:Nf
    
    % column vector dependant on u
    a = k/2*sqrt(E*h/(rho*L))*Dxx*u1;

    % A and C for state update: A*u = B*u1-C*u2+J*f(n)
    C = (sig*k-1)*speye(N-1)-a*a'-2*sig1*k*Dxx;
    % inverse of A with ShermanMorrison-Woodbury formula
    A_inv = 1/(1+sig*k)*(speye(N-1)-a*a'/(1+sig*k+a'*a)); 

    % update state, and insert current value of f.
    u = A_inv*(B*u1 + C*u2 + J*f(n));
    u0 = B0*u10 + C0*u20 + J*f(n);
    
    % read output
    y(n) = cout*u;
    y0(n) = cout*u0;
    
    % plot
    if(plot_on==1)
        % draw current state
        figure(1)
        plot(xax, u, 'k');
        xlabel('x (m)')
        ylabel('u (m)')
        axis([0 L -0.001 0.001])
        hold on
        plot(xax, u0, 'b');
        xlabel('x (m)')
        ylabel('u (m)')
        axis([0 L -0.001 0.001])
        hold off
        legend('nonlinear string','linear string')
        drawnow

        % Capture the current frame and write it to the video object
        currentFrame = getframe(gcf);
        writeVideo(videoObj, currentFrame);
        
    end
    
    % shift state
    
    u2 = u1;
    u1 = u;

    u20 = u10;
    u10 = u0;
    
end
close(videoObj);

%/////////////////////////////////////////////////////////////////////////%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hammer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% variables for the hammer
xH0 = -0.001; vH0 = 2;    % initial conditions of hammer
xH00 = 0.001; vH00 = 2;   % another initial conditions of hammer

xi1 = 0.52;               % coordinate of pluck excitation (normalised, 0-1)
xH1 = 0.3;                % coordinate of hammer excitation (normalised, 0-1)

exc_st1 = 0.3;            % start time of pluck excitation (s)
exc_st2 = 0.0;            % start time of pluck excitation (s)

xo1 = 0.3;                % coordinate of output (normalised, 0-1)

famp1 = 2;                % peak amplitude of excitation (N)
% dur = 0.001;            % duration of excitation (s) (already included before)

%% calculate output
[yh1,fh1] = hammerstring(ploth1,itype,SR,Tf,xi1,xH1,xo1,xH0,vH0,exc_st1,famp1);
[yh2,fh2] = hammerstring(ploth2,itype,SR,Tf,xi1,xH1,xo1,xH00,vH00,exc_st1,famp1);
[yh3,fh3] = hammerstring(ploth3,itype,SR,Tf,xi1,xH1,xo1,xH0,vH0,exc_st2,famp1);
[yh4,fh4] = hammerstring(ploth4,itype,SR,Tf,xi1,xH1,xo1,xH00,vH00,exc_st2,famp1);

%% play sound

soundsc([y;y0;yh1;yh2;yh3;yh4],SR);

%% plot spectrum

figure(2)
subplot(4,2,1)
plot([0:Nf-1]*k, yh1, 'k'); 
title('Position of Target Mass with negative initial hammer position and delayed pluck'); 
xlabel('t'); ylabel('y');
axis([0 0.5 -0.005 0.005])
subplot(4,2,2)
plot([0:Nf-1]*k, fh1, 'k'); 
title('Hammer Force/Mass with negative initial hammer position and delayed pluck'); 
xlabel('t'); ylabel('F');
axis([0 0.05 0 3000])
subplot(4,2,3)
plot([0:Nf-1]*k, yh2, 'k'); 
title('Position of Target Mass with positive initial hammer position and delayed pluck'); 
xlabel('t'); ylabel('y');
axis([0 0.5 -0.005 0.005])
subplot(4,2,4)
plot([0:Nf-1]*k, fh2, 'k'); 
title('Hammer Force/Mass with positive initial hammer position and delayed pluck'); 
xlabel('t'); ylabel('F');
axis([0 0.05 0 3000])
subplot(4,2,5)
plot([0:Nf-1]*k, yh3, 'k'); 
title('Position of Target Mass with negative initial hammer position and simultaneous pluck'); 
xlabel('t'); ylabel('y');
axis([0 0.5 -0.005 0.005])
subplot(4,2,6)
plot([0:Nf-1]*k, fh3, 'k'); 
title('Hammer Force/Mass with negative initial hammer position and simultaneous pluck'); 
xlabel('t'); ylabel('F');
axis([0 0.05 0 3000])
subplot(4,2,7)
plot([0:Nf-1]*k, yh4, 'k'); 
title('Position of Target Mass with positive initial hammer position and simultaneous pluck'); 
xlabel('t'); ylabel('y');
axis([0 0.5 -0.005 0.005])
subplot(4,2,8)
plot([0:Nf-1]*k, fh4, 'k'); 
title('Hammer Force/Mass with positive initial hammer position and simultaneous pluck'); 
xlabel('t'); ylabel('F');
axis([0 0.05 0 3000])


figure(3)
subplot(2,1,1)
yfft = 10*log10(abs(fft(y)));
plot([0:Nf-1]'*(SR/Nf), yfft, 'k')
xlim([0 SR/2])
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Transform of output (tension nonlinear string)')
subplot(2,1,2)
yfft0 = 10*log10(abs(fft(y0)));
plot([0:Nf-1]'*(SR/Nf), yfft0, 'k')
xlim([0 SR/2])
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Transform of output (linear string)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function for stiff string with frequency independent and dependent loss
function [y,fh] = hammerstring(plot_on,itype,SR,Tf,xi,xH,xo,xH0,vH0,exc_st,famp)
% this function works with a hammer string simulation with also loss and
% pluck/strike force, see input/output parameter definition as in the main
% script

% physical string parameters ----------------------------------------------
T = 240;                    % tension (N)
r = 0.00043;                % string radius (m)
E = 2e11;                   % Young's modulus (Pa) for steel
rho = 9000;                 % density (kg/m^3)
L = 0.65;                   % length (m)
loss = [100, 10; 1000, 8];  % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
dur = 0.001;                % duration of pluck/strike excitation (s)

% parameters for hammer ---------------------------------------------------

MR = 100;                   % parameter for force in function, MH/(rho*A*L)
wH = 1000;                  % stiffness parameter for hammer
alpha = 2;                  % hammer stiffness nonlinearity exponent

% derived parameters ------------------------------------------------------

A = pi*r^2;                 % string cross-sectional area
I = 0.25*pi*r^4;            % string moment of inertia

c = sqrt(T/(rho*A));        % wave speed
K = sqrt(E*I/(rho*A));      % stiffness constant 

k = 1/SR;                   % time step
Nf = floor(SR*Tf);          % number of time steps

%loss parameters
zeta1 = (-c^2+sqrt(c^4+4*K^2*(2*pi*loss(1,1))^2))/(2*K^2);
zeta2 = (-c^2+sqrt(c^4+4*K^2*(2*pi*loss(2,1))^2))/(2*K^2);
sig = 6*log(10)*(-zeta2/loss(1,2)+zeta1/loss(2,2))/(zeta1-zeta2);
sig1 = 6*log(10)*(1/loss(1,2)-1/loss(2,2))/(zeta1-zeta2);

% grid spacing ------------------------------------------------------------
hmin = sqrt(0.5*(c^2*k^2+sqrt(c^4*k^4+16*K^2*k^2)));
N = floor(L/hmin);
h = L/N;

% input -------------------------------------------------------------------
t = [0:Nf-1]'/SR;           % time vector
ind = sign(max(-(t-exc_st).*(t-exc_st-dur),0));     % select time period
f = 0.5*famp*ind.*(1-cos(itype*pi*(t-exc_st)/dur)); % calculate force

% matrices for update -----------------------------------------------------
% u = B*u1-C*u2+J*f(n)

% input/output position in grid point index
xi_m = floor(N*xi);
xo_m = floor(N*xo);     
xH_m = floor(N*xH);     

% matrices of spatial derivatives
e = ones(N,1);
Dx = spdiags([-1*e e], -1:0,N,N-1)/h;
Dxx = -Dx'*Dx;
Dxxxx = Dxx*Dxx;

% matrices for state update
B = 1/(1+sig*k)*(2*speye(N-1)+c^2*k^2*Dxx-K^2*k^2*Dxxxx+2*sig1*k*Dxx);
C = 1/(1+sig*k)*((-1+sig*k)*speye(N-1)-2*sig1*k*Dxx);

% input force vector
J = sparse(N-1,1);
J(xi_m) = 1/(1+sig*k)*k^2/(rho*A*h);

% hammer input vector (selects position)
d = sparse(1,N-1);              
d(xH_m) = 1/(1+sig*k);

% output vector
c = sparse(1,N-1);              
c(xo_m) = 1;


% initialise scheme variables ---------------------------------------------

% hammer state
uH2 = xH0; uH1 = xH0+k*vH0;
% string state
u2 = zeros(N-1,1); u1 = u2; u = u2; 

fh = zeros(Nf,1);               % hammer force
y = zeros(Nf,1);                % output
xax = (1:N-1)'*h;               % x-axis for plotting


% main loop ---------------------------------------------------------------

for n=1:Nf
    
    % calculate hammer force
    fh(n) = wH^(1+alpha)*(0.5*((uH1-u1(xH_m))+abs((uH1-u1(xH_m)))))^alpha;
    
    % update hammer state
    uH = 2*uH1-uH2-k^2*fh(n);
    uH2 = uH1; uH1 = uH;

    % update string state
    u = B*u1 + C*u2 + J*f(n) + d'*MR*k^2*fh(n);
    
    % read output
    y(n) = c*u;
    
    % plot
    if(plot_on==1)
        % draw current state
        figure(1)
        plot(xax, u, 'k');
        xlabel('x (m)')
        ylabel('u (m)')
        axis([0 0.65 -0.005 0.005])
        drawnow
    end
    
    % shift state
    u2 = u1;
    u1 = u;  
end

end


