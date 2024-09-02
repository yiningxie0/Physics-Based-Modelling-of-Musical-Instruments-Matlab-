%-------------------------------------------------
% PBMMI Assignment 3 BTB1 Nonlinear Moog VCF 
% - with Trapezoidal Integration and Newton-Raphson algorithm
%
% Author: B200247 
% Date: 24/2/23
%
% This script calculates Nonlinear Moog VCF with
%
% Trapezoidal Integration (function, starting line 210)
% Newton-Raphson algorithm (function, starting line 227)
%
%
% And plots 
% 1 a) a graph to show system state variables
% 1 b) a comparison of outputs in db between nonlinear and linear case with
%      0.1 and 1 amplitude impulse input - showing that smaller amplitude
%      makes result closer with the linear one   
% 
% 2 graphs with sinusoidal input under different amplitudes (0.001,1,100),
%   2a) similar to linear response in low amplitude
%   2b) harmonic distortion under moderate amplitude
%   2c) aliasing under high amplitude
%   out of sync between frequency peaks and over all trend is influenced by
%   feedback coeff: r closer to 1 results in smoother transition as on slides
%
% Not sure if writing this way runs correctly on every machine, but if it
% does not, you may try the latest Matlab online version which should be sure to run
%-------------------------------------------------
clear all; clc; close all;  

%% input parameters

SR = 44100;                  % sample rate [Hz]
Tf = 0.2;                    % total simulation time [s]
f0 = 120;                    % resonant filter frequency [Hz]
r = 0.7;                     % feedback coeff 0<=r<=1

% changed parameters to present aliasing for figure 2
Tf2 = 0.5;                  
f1 = 1000;                   

%% derived parameters

k = 1/SR;                    % time step [s]
om0 = 2*pi*f0;               % angular resonant frequency [rad/s]
om1 = 2*pi*f1;               % angular resonant frequency [rad/s]
Nf = floor(Tf*SR);           % total number of samples
Nf2 = floor(Tf2*SR);         % total number of samples


%% stability check

% stable with 0<=r<=1


%% initialize

% parameters for state-space form system
I = eye(4);
A = om0*[-1 0 0 -4*r; 1 -1 0 0; 0 1 -1 0; 0 0 1 -1];
b = om0*[1;0;0;0];
c = [0 0 0 1];

% for figure 1 -----------------------------
% input - 0-pad to include u(0) = 0 for update of x(1)
u1 = 0.1*[0;1;zeros(Nf-1,1)];
u0 = [0;1;zeros(Nf-1,1)];

% state variables
x1 = zeros(4,1);
x0 = zeros(4,1);

% overall state variable vector in case we need to obverse its change
% take the u0 input setting as example
xy0 = zeros(4,Nf);

% output
y1 = zeros(Nf,1);
y0 = zeros(Nf,1);

tvec = [0:Nf-1]'*k;         % time vector for plots
fvec = [0:Nf-1]'*SR/Nf;     % frequency vector for plots

Hc = zeros(Nf,1);           % exact transfer function in the linear case


% for figure 2 -----------------------------
% sinusoid input
ua0 = [0;0.001*sin(2*pi*1000*[0:Nf2-1]'/SR)];
ua1 = [0;10*sin(2*pi*1000*[0:Nf2-1]'/SR)];
ua2 = [0;100*sin(2*pi*1000*[0:Nf2-1]'/SR)];

% state variables
xa0 = zeros(4,1);
xa1 = zeros(4,1);
xa2 = zeros(4,1);

% output
ya1 = zeros(Nf2,1);
ya0 = zeros(Nf2,1);
ya2 = zeros(Nf2,1);

tvec2 = [0:Nf2-1]'*k;       % time vector for plots
fvec2 = [0:Nf2-1]'*SR/Nf2;  % frequency vector for plots

%% main calculation

% ----------------------- loop for figure 1 -------------------------------
tic
for n = 1:Nf

    % calculate system functions from Trapezoid rule
    % g: numerical function by Trapezoid rule
    % J: Jacobian of g
    [J1,g1] = myTrapz(x1,u1(n+1),u1(n),k,om0,r);
    [J0,g0] = myTrapz(x0,u0(n+1),u0(n),k,om0,r);

    % inplement Newton Raphson method
    x1 = NR(J1,g1);
    x0 = NR(J0,g0);

    % store state
    xy0(:,n)=x0;
    
    % write output
    y1(n) = c*x1;
    y0(n) = c*x0;

    % transfer function for linear VCF
    Hc(n) = c*((1j*2*pi*fvec(n)*I-A)\b); 

end

simTime = toc;

% calculate outputs in frequency domain
Yt1 = fft(y1);
Yt0 = fft(y0);
Yc1 = 0.1*Hc;  % analytical output for linear VCF with 0.1 amp impulse input
Yc0 = Hc;      % analytical output for linear VCF with 1 amp impulse input

% ----------------------- loop for figure 2 -------------------------------
for n = 1:Nf2

    % calculate system functions from Trapezoid rule
    % g: numerical function by Trapezoid rule
    % J: Jacobian of g
    [Ja0,ga0] = myTrapz(xa0,ua0(n+1),ua0(n),k,om1,r);
    [Ja1,ga1] = myTrapz(xa1,ua1(n+1),ua1(n),k,om1,r);
    [Ja2,ga2] = myTrapz(xa2,ua2(n+1),ua2(n),k,om1,r);

    % inplement Newton Raphson method
    xa1 = NR(Ja1,ga1);
    xa0 = NR(Ja0,ga0);
    xa2 = NR(Ja2,ga2);

    % write output
    ya1(n) = c*xa1;
    ya0(n) = c*xa0;
    ya2(n) = c*xa2;
end


%% plots

% plots to show system state and outputs under impulse input with different
% amplitudes  

figure
subplot(2,1,1)
t = [0:Nf-1]'*k;
plot(t,xy0);
xlabel('Time (s)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('State Variables, $f_0 = 120Hz$, $r = 0.7$, impulse input', 'Interpreter', 'latex')
legend('x1','x2','x3','x4', 'Interpreter', 'latex')

subplot(2,1,2)
loglog(fvec(1:ceil(Nf/2)), abs(Yt0(1:ceil(Nf/2))));
hold on
loglog(fvec(1:ceil(Nf/2)), abs(Yc0(1:ceil(Nf/2)))); 
hold on
loglog(fvec(1:ceil(Nf/2)), abs(Yt1(1:ceil(Nf/2))));
hold on
loglog(fvec(1:ceil(Nf/2)), abs(Yc1(1:ceil(Nf/2))));
xlabel('$f (Hz)$', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('Outputs, $f_0 = 120Hz$, $r = 0.7$, impulse input', 'Interpreter', 'latex')
legend('Non-linear VCF with 1 amp impulse','Linear VCF with 1 amp impulse','Non-linear with 0.1 amp impulse','Linear VCF with 0.1 amp impulse', 'Interpreter', 'latex')
xlim([0 1e4])
ylim([1e-8 100])

% plots to show the aliasing
figure
subplot(3,1,1)
semilogy(fvec2,abs(fft(ya0)))
xlim([0 ceil(SR/2)])
xlabel('$f (Hz)$', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('Moog VCF response, $f_0 = 1000Hz$, $r = 0.7$, sinusoid input with 0.001 amplitude', 'Interpreter', 'latex')

subplot(3,1,2)
semilogy(fvec2,abs(fft(ya1)))
xlim([0 ceil(SR/2)])
xlabel('$f (Hz)$', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('Moog VCF response, $f_0 = 1000Hz$, $r = 0.7$, sinusoid input with 10 amplitude', 'Interpreter', 'latex')

subplot(3,1,3)
semilogy(fvec2,abs(fft(ya2)))
xlim([0 ceil(SR/2)])
xlabel('$f (Hz)$', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('Moog VCF response, $f_0 = 1000Hz$, $r = 0.7$, sinusoid input with 100 amplitude', 'Interpreter', 'latex')

%% function for Trapezoid rule
function [J,g] = myTrapz(x,u1,u2,k,om0,r)
% function to store variable functions involved in Trapezoid rule
% inputs state vectors, system inputs, time step, angular resonant frequency, resonant coeff
% outputs the numerical function by Trapezoid rule and its Jacobian

% system function
f = @(x,u) om0*[-tanh(x(1))+tanh(u-4*r*x(4)); -tanh(x(2))+tanh(x(1)); -tanh(x(3))+tanh(x(2)); -tanh(x(4))+tanh(x(3))];
% Trapezoid rule
g = @(w) w-x-k/2*(f(w,u1)+f(x,u2));
% Jacobian of g
J = @(w) eye(4)-k/2*om0*[-sech(x(1))^2 0 0 -4*r*sech(u1-4*r*x(4)).^2;
                            sech(x(1)).^2 -sech(x(2)).^2 0 0;
                            0 sech(x(2)).^2 -sech(x(3)).^2 0; 
                            0 0 sech(x(3)).^2 -sech(x(4)).^2];
end

%% function for Newton Raphson method
function w = NR(J,g)
% function to implement Newton Raphson method to find roots
% inputs function to calculate and its Jacobian
% outputs the root

% % test code to estimate roots
% numPoints = 1000;
% xax = [1;1;1;1]*linspace(-5, 5, 1000);
% plot(xax, g(xax))
% grid on

% Newton Raphson method 
w = 0*ones(4,1);          % initial guess
maxIter = 10;             % max number of guesses
iter = 0;                 % iteration counter
tol = 1e-10;
step = ones(4,1);

while (iter < maxIter) && (abs(min(step)) > tol) 
   step = J(w)\g(w);      % calculate step size
   w = w - step;          % step update
   iter = iter + 1;       % iter update
end
end
