%-------------------------------------------------
% PBMMI Assignment 2 BTB2 Nonlinear Oscillators 
%
% Author: B200247 
% Date: 11/2/23
%
% This script includes the following progressively:
% 
% a) Basic nonlinear cubic oscillator following scheme (4.9b) on NSS,
%    unconditionally stable with numerical energy calculation 
%
% b) The mixed linear/cubic oscillator with numerical energy calculation
%    and a naive stability check      
%
% c) Non-drive Duffing Oscillator with finite difference time domain method
%
% d) Non-drive Duffing Oscillator with Runge–Kutta 4th-order method (as
%    read that it has relatively better accuracy) 
%
% e) A drum simulation with granular-like sound by Duffing Oscillator
%
% Outputs all the above sounds one by one
%
% reference: "Numerical Sound Synthesis: Finite Difference Schemes and Simulation in Musical Acoustics”. Stefan Bilbao, Wiley (2009).
%            https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Derivation_of_the_Runge%E2%80%93Kutta_fourth-order_method
%-------------------------------------------------

clear all; clc; close all;  


%% basic parameters

SR = 44100;                  % sample rate [Hz]
Tf = 1;                      % total duration [s]
f0 = 200;                    % reference frequency [Hz]
u0 = 1.0;                    % initial oscillator position
v0 = 0.0;                    % initial oscillator velocity
w1 = 20*pi;                  % parameter for cubic term

% parameters for duffing oscillator, here I follow the form of:
%
% d^u/dt^2 = -sigma*du/dt-alpha*u-beta*u^3
%
% though some can be shared with the other 2 - alpha corresponds to wo^2
% and beta to w1^4, in this equation they aren't necessarily positive,
% hence I set them independently though still with the corresponding value
% as described.
% by changing them smaller as more often seen, you might need to modify
% time step k to smaller values as well

alpha = (2*pi*f0)^2;         % spring force parameter
beta = 1e8;                  % nonlinearity parameter
sigma = 6*log(10);           % damping parameter

% sigma may also be calculated with sigma = 6*log(10)/T60, here uses the 
% 20log10 db definition

%% derived parameters

k = 1/SR;                    % time step [s]
w0 = 2*pi*f0;                % reference angular frequency [rad/s]
Nf = floor(Tf*SR);           % total number of steps


%% initialize

% for nonlinear cubic oscillator
uc_2 = u0;                   % initial displacement
uc_1 = u0+k*v0;              % second displacement
outc = zeros(Nf,1);          % output vector

Hc = zeros(Nf,1);            % total energy
Hc_k = zeros(Nf,1);          % kinetic energy
Hc_p = zeros(Nf,1);          % potential energy

% for mixed linear/cubic oscillator
um_2 = u0;                   % initial displacement
um_1 = u0+k*v0;              % second displacement
outm = zeros(Nf,1);          % output vector

Hm = zeros(Nf,1);            % total energy
Hm_k = zeros(Nf,1);          % kinetic energy
Hm_p = zeros(Nf,1);          % potential energy

% for Duffing Oscillator
ud_2 = u0;                   % initial displacement
ud_1 = u0+k*v0;              % second displacement
outd = zeros(Nf,1);          % output vector

% for Duffing Oscillator with Runge–Kutta fourth-order method
uD = zeros(Nf,1);            % displacement vector
vD = zeros(Nf,1);            % velocity vector

uD(1)=u0;                    % initial displacement
vD(1)=v0;                    % initial velocity

f = @(vD)vD;                                % du/dt
g = @(uD,vD)-sigma*vD-alpha*uD-beta*uD^3;   % right-hand side of EOM


%% ///////////////// a) Basic nonlinear cubic oscillator //////////////////

for n=1:Nf
    outc(n) = uc_2;          % allocate output
    
    % Energy Calculation
    Hc_k(n) = 0.5/k^2*(uc_1-uc_2)^2; 
    Hc_p(n) = w1^4/4*uc_1^2*uc_2^2;
    Hc(n) = Hc_k(n)+Hc_p(n);

    u = (2/(1+k^2*w1^4*uc_1^2/2))*uc_1-uc_2;      % calculate difference scheme
    uc_2 = uc_1;                                  % update state of difference scheme
    uc_1 = u;                                     % update state of difference scheme    
end

% soundsc(outc,SR);

%% ///////////////// b) The mixed linear/cubic oscillator /////////////////

for n=1:Nf
    outm(n) = um_2;          % allocate output
    
    % Energy Calculation
    Hm_k(n) = 0.5/k^2*(um_1-um_2)^2; 
    Hm_p(n) = w1^4/4*um_1^2*um_2^2+w0^2/2*um_1*um_2; 
    Hm(n) = Hm_k(n)+Hm_p(n);
    
    % Stability Check - though more efficient to pre-calculate the conditions and check outsie
    assert(Hm(n)>=0,"Error: unstable");

    u = (2*(2-k^2*w0^2)/(2+k^2*w1^4*um_1^2))*um_1-um_2; % calculate difference scheme
    um_2 = um_1;                                  % update state of difference scheme
    um_1 = u;                                     % update state of difference scheme    
end

% soundsc(outm,SR);

%% ///////////////// c) Duffing Oscillator with FDTD //////////////////////

for n=1:Nf
    outd(n) = ud_2;          % allocate output
    
    u = 2*(2-k^2*alpha)/(2+k*sigma+k^2*beta*ud_1^2)*ud_1-(2-k*sigma+k^2*beta*ud_1^2)/(2+k*sigma+k^2*beta*ud_1^2)*ud_2; % calculate difference scheme
    ud_2 = ud_1;                                  % update state of difference scheme
    ud_1 = u;                                     % update state of difference scheme    
end

% soundsc(outd,SR);

%% ////// d) Duffing Oscillator with Runge–Kutta fourth-order method //////

for n=1:Nf-1

    % calculate increments
    k1 = f(vD(n));
    l1 = g(uD(n),vD(n));

    k2 = f(vD(n)+k*l1/2);
    l2 = g(uD(n)+k*k1/2,vD(n)+k*l1/2);

    k3 = f(vD(n)+k*l2/2);
    l3 = g(uD(n)+k*k2/2,vD(n)+k*l2/2);

    k4 = f(vD(n)+k*l3);
    l4 = g(uD(n)+k*k3,vD(n)+k*l3);

    % update uD and vD
    uD(n+1) = uD(n)+k*(k1+2*k2+2*k3+k4)/6;
    vD(n+1) = vD(n)+k*(l1+2*l2+2*l3+l4)/6;
end

% soundsc(uD,SR);

%% //////// e) make sound with the Duffing Oscillator in d) ///////////////

% generate a sound
outD = uD.*(sin(2*pi*uD*SR));

% go through a flanger
M = 60;                 % depth of flanger 
outD(M:Nf)=outD(M:Nf)+0.5*outD(1:Nf-M+1);

outD = outD./max(abs(outD));
sound(outD,SR);
audiowrite('drum.wav',outD,SR);

%% make sound

% make a sequence from the above sound
outT = [outc;zeros(SR*0.5,1);outm;zeros(SR*0.5,1);outd;zeros(SR*0.5,1);uD;zeros(SR*0.5,1);outD];

soundsc(outT,SR);


%% plot
% calculate time and frequency vector
t = [0:Nf-1]'*k;
f = SR*[0:Nf-1]/Nf;

% plot the outputs and energy of a) cubic and b) mixed oscillator
figure
subplot(3,2,1)
plot(t,outc);
xlabel('t');
ylabel('u');
title({'Nonlinear Cubic Oscillator in Time Domain';[]});

subplot(3,2,3)
loglog(f,abs(fft(outc)));
title({'Nonlinear Cubic Oscillator in Frequency Domain';'with Runge–Kutta 4th-order method'});
xlabel('frequency (Hz)');
ylabel('transform');
xlim([0 SR/2]);

subplot(3,2,5)
plot(0:Nf-1,(Hc-Hc(1))/Hc(1));
xlabel('t');
ylabel('error');
title({'Error in Energy for Nonlinear Cubic Oscillator';[]});

subplot(3,2,2)
plot(t,outm);
xlabel('t');
ylabel('u');
title({'Mixed Linear/Cubic Oscillator in Time Domain';[]});

subplot(3,2,4)
loglog(f,abs(fft(outm)));
title({'Mixed Linear/Cubic Oscillator in Frequency Domain';[]});
xlabel('frequency (Hz)');
ylabel('transform');
xlim([0 SR/2]);

subplot(3,2,6)
plot(0:Nf-1,(Hm-Hm(1))/Hm(1));
xlabel('t');
ylabel('error');
title({'Error in Energy for Mixed Linear/Cubic Oscillator';[]});

% plot the outputs and phase portraits of Duffing Oscillator in c) and d)
figure
subplot(3,2,1)
plot(t,outd);
xlabel('t');
ylabel('u');
title({'Duffing Oscillator with FDTD in Time Domain';[]});

subplot(3,2,3)
loglog(f,abs(fft(outd)));
title({'Duffing Oscillator with FDTD in Frequency Domain';[]});
xlabel('frequency (Hz)');
ylabel('transform');
xlim([0 SR/2]);

subplot(3,2,5)

% calculate velocity
vd = zeros(length(outd),1);
vd(2:length(outd)) = 1/k*(outd(2:length(outd)) - outd(1:length(outd)-1));

plot(outd,vd)
xlabel('Displacement')
ylabel('Velocity')
title({'Phase Portrait of Duffing Oscillator with FDTD';[]});

subplot(3,2,2)
plot(t,uD);
xlabel('t');
ylabel('u');
title({'Duffing Oscillator with Runge–Kutta 4th-order method';'in Time Domain'});

subplot(3,2,4)
loglog(f,abs(fft(uD)));
title({'Duffing Oscillator with Runge–Kutta 4th-order method';'in Frequency Domain'});
xlabel('frequency (Hz)');
ylabel('transform');
xlim([0 SR/2]);

subplot(3,2,6)
plot(uD, vD)
xlabel('Displacement')
ylabel('Velocity')
title({'Phase Portrait of Duffing Oscillator';'with Runge–Kutta 4th-order method'});

% plot modulized output
figure
subplot(3,1,1)
plot(t,outD);
xlabel('t');
ylabel('u');
title({'A drum simulation with granular-like sound in time domain'});

subplot(3,1,2)
loglog(f,abs(fft(outD)));
title({'A drum simulation with granular-like sound in frequency domain'});
xlabel('frequency (Hz)');
ylabel('transform');
xlim([0 SR/2]);

subplot(3,1,3)
% calculate velocity
vv = zeros(length(outD),1);
vv(2:length(outD)) = 1/k*(outD(2:length(outD)) - outD(1:length(outD)-1));

plot(outD, vv)
xlabel('Displacement')
ylabel('Velocity')
title({'A drum simulation with granular-like sound phase portrait'});

