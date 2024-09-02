% ------------------------------------------------
%                       PBMMI Assignment 6 BTB1 
%     Generalised state-space model, Asymmetry and Damped Newton's
%   - applied to Diode Clipper, Tube Screamer and Dallas Rangemaster.
% 
%                           Author: B200247 
%                            Date: 14/4/23
%
% This code is generalised to accept any set of state space arrays A − N,
% and any nonlinear function f(v) - at least proven with Diode Clipper,
% Tube Screamer and Dallas Rangemaster. See flag 'circuit'. 
%
% Parameter arrays for discrete time and state update methods follow [1],
% circuit parameters from [1] and [2], DC initialization method for Dallas
% Rangemaster follows [2] 
%
% In detail, following the above methods which look slightly different from
% that on slide: 
%
% - Modelled asymmetric clipping
% - Modelled asymmetric Diode Clipper
% - Modelled asymmetric Tube Screamer 
% - Modelled Dallas Rangemaster.
% - Included Dapped Newton's Method in Newton Raphson function starting
%   from line 230
% - Generalised to allow all above (and other possible systems) to work in
%   one same main loop, including the same function for Newton's Method
% - Similar with basic, plots input and output in time domain, and output
%   signal in frequency domain. 
%
% Instruction ------------------------------------------------------------
% input = 0: coded sinusoid, 1: non-sinusoidal audio
% circuit = 1: Diode Clipper, 2: Tube Screamer, 3: Dallas Rangemaster
%
% Reference --------------------------------------------------------------
% [1] Carson, Alistair. (2020). Aliasing Reduction in Virtual Analogue Modelling.
% [2] Porter, Mac. (2019). Virtual Analog Modeling of Guitar Effects Circuits.
% ------------------------------------------------
clear all; clc; close all;  

%% input type
% 0: sinusoidal, 1: read audio
input = 0;
% circuit type: 1: Diode Clipper, 2: Tube Screamer, 3: Dallas Rangemaster
circuit = 1;
%% input and simulation parameters
if input == 1

    % --------- non-sinusoidal audio ------------
    [u,SR] = audioread('A6_PG_B200247_Input.wav');  % read audio, SR: sample rate [Hz] u: sample 
    [Nf,n] = size(u);                      % get length and number of channels of x
    if n ~=1
        u = sum(u,2)/n;                    % Average channels to mono
    end
    % would like to adjust amp
    amp = 1;
    u = amp*u;
    % time vector [s]
    tvec = (0:Nf-1)'/SR; 
else
    % ----------- sinusoidal input --------------
    SR = 88.2e3;                           % sample rate [Hz]
    dur = 0.3;                             % duration of sim [s]
    Nf = round(SR*dur);                    % number of samples in sim

    f0 = 200;                              % frequency [Hz]
    amp = 2;                               % amplitude of input
    tvec = (0:Nf-1)'/SR;                   % time vector [s]
    u = amp*sin(2*pi*f0*tvec);             % input vector
end


%% Physical Parameters
switch circuit
    case 3
        Is = 2.029e-6;                     % saturation current
        Bf = 108.1;                        % forward current gain
        Br = 8.602;                        % reverse current gain
        % resistance
        R1 = 470e3;R2 = 68e3;R3 = 10e3;
        R4 = 3.9e3;R5 = 1e6;
        % capacitance
        C1 = 4.7e-9;C2 = 47e-6;C3 = 10e-9;
    case 2
        Is = 2.52e-9;                      % saturation current
        % resistance
        R1 = 10e3; R2 = 551e3; R3 = 4.7e3;
        % capacitance
        C1 = 1e-6; C2 = 51e-12; C3 = 47e-9; 
    case 1
        Is = 2.52e-9;                      % saturation current
        r = 1e3;                           % resistance
        c = 33e-9;                         % capacitance
end

Vt = 25.83e-3;                             % thermal voltage (V)
Ni = 1.752;                                % ideality factor


%% K-method parameters
k = 1/SR; % time step

% parameter arrays
switch circuit
    case 3
        A = [-(R1+R2)/(C1*R1*R2) 0 0; 0 -1/(C2*R4) 0; 0 0 -1/(C3*(R3+R5))];
        B = [(R1+R2)/(C1*R1*R2) -1/(C1*R1); 0 0; 0 1/(C3*(R3+R5))];
        C = [-1/C1 0; -1/C2 -1/C2; 0 R3/(C3*(R3+R5))];
        D = [0 0 -R5/(R3+R5)];
        E = [0 R5/(R3+R5)];
        F = [0 R3*R5/(R3+R5)];
        G = [1 1 0; 1 0 R3/(R3+R5)];
        H = [-1 0; -1 R5/(R3+R5)];
        K = [0 0; 0 R3*R5/(R3+R5)];
        I = eye(size(A));
        L = Is*[1/Bf 1/Br; 1 -(Br+1)/Br];
    case 2
        A = -[1/(R1*C1), 0, 0; 1/(R3*C2), 1/(R2*C2), 1/(R3*C2); 1/(R3*C3), 0, 1/(R3*C3)]; 
        B = [1/(R1*C1); 1/(R3*C2); 1/(R3*C3)]; 
        C = [0; -1/C2; 0]; 
        D = [-1, 1, 0]; E = 1; F = 0; 
        G = [0, 1, 0]; H = 0; K = 0;
        I = eye(size(A));
    case 1
        A = -1/(r*c);
        B = 1/(r*c);
        C = -1/c;
        D = 1; E = 0; F = 0;
        G = 1; H = 0; K = 0; I = eye(size(A));
end

% discrete time model
Z = (2*SR*I-A)\I;
A_ = (2*SR*I+A)*Z;
B_ = 2*Z*B;
C_ = 2*Z*C;
D_ = 2*SR*D*Z;
E_ = E+D*Z*B;
F_ = F+D*Z*C;
G_ = 2*SR*G*Z;
H_ = H+G*Z*B;
K_ = K+G*Z*C;

%% functions
% i: nonlinear function
% g: implicit function to solve v
% J: Jacobian of g, for Newton's solver
switch circuit
    case 3
        i = @(v)L*(exp(v/Vt)-1);
        g = @(p,v) p+K_*L*(exp(v/Vt)-1)-v; % 
        J = @(v) K_/Vt*L.*exp([v v]'/Vt)-eye(length(v));
    case 2
        i = @(v)Is*(exp(v/(2*Vt*Ni))-exp(-v/(Vt*Ni)));
        g = @(p,v) K_*Is*(exp(v/(2*Vt*Ni))-exp(-v/(Vt*Ni)))+p-v; 
        J = @(v) K_*(Is/(Vt*Ni))*(0.5*exp(v/(2*Vt*Ni))+exp(-v/(Vt*Ni)))-1;
    case 1
        i = @(v)Is*(exp(v/(2*Vt*Ni))-exp(-v/(Vt*Ni)));
        g = @(p,v) K_*Is*(exp(v/(2*Vt*Ni))-exp(-v/(Vt*Ni)))+p-v; 
        J = @(v) K_*(Is/(Vt*Ni))*(0.5*exp(v/(2*Vt*Ni))+exp(-v/(Vt*Ni)))-1;
end

%% initialize
% initialize output
y = zeros(Nf,1);
% initialize input
switch circuit
    case 3
        % input
        in = [u -9*ones(length(u),1)]';    
        % initialize for steady state
        % DC matrices
        Adc = [-R2/(R1+R2); 0; 1]; 
        Bdc = [-R1*R2/(R1+R2) 0; -R4 -R4; 0 R3]; 
        Cdc = [-R2/(R1+R2); R1/(R1+R2)]; 
        Ddc = [((-R1-R4)*R2-R1*R4)/(R1+R2) (-R1*R4-R2*R4)/(R1+R2);-R1*R2/(R1+R2) (R1*R3+R2*R3)/(R1+R2)]; 
        % initialize for DC analysis
        u = -9; 
        v = [0;0]; 
        f = @(v)Cdc*u+Ddc*L*(exp(v/Vt)-1)-v;
        j = @(v)Ddc/Vt*L.*exp([v v]'/Vt)-eye(length(v));
        % steady state voltage
        v = NR(j,f,v);
        % initialized x
        x0 = Adc*u+Bdc*L*(exp(v/Vt)-1);
        % translate x to xc (as in [1])
        x = 1/(2*SR)*((2*SR*I+A)*x0+B*[0;-9]+C*L*(exp(v/Vt)-1));
    case 2
        in = u';
        x = [0;0;0];v = 0;
    case 1
        in = u';
        x = 0;v = 0;
end

%% main loop
for n = 1:Nf
    p = G_*x+H_*in(:,n);             % calculate p
    gnew = @(v)g(p,v);               % replace p variable with current value p
    v = NR(J,gnew,v);                % root search for v
    y(n) = D_*x+E_*in(:,n)+F_*i(v);  % output 
    x = A_*x+B_*in(:,n)+C_*i(v);     % state update
end


%% make sound
soundsc(y,SR);

%% plot
% Input/Output signal in Time Domain
subplot(2,1,1)
plot(tvec,in(1,1:end),tvec,y);
xlabel('Time(s)');                         
ylabel('Amplitude');  
xlim([0.03 0.05])
legend('input','output')
title('Input/Output signal in Time Domain');   

% Output signal in Frequency Domain
subplot(2,1,2)
f = (0:ceil(Nf/2)-1)/Nf*SR;              % frequency vector
Y = fft(y,Nf);                           % y after fft
Ydb = 20*log10(abs(Y(1:ceil(Nf/2),:)));
plot(f,Ydb);
xlabel('Frequency(Hz)');                         
ylabel('dB');             
title('Output signal in Frequency Domain');   


%% Damped Newton's method
function v = NR(J,g,v0)
% Inputs functions to solve and its Jacobian, and initial guess value
% Newton Raphson method 
v = v0;                   % initial guess
maxIter = 50;             % max number of guesses
maxDIter = 10;            % max number of guesses in damping
iter = 0;                 % iteration counter
tol = 1e-9;               % tolerance
step = ones(length(v),1); % initial step
f = g(v);                 % function to solve
df = J(v);                % Jacobian
f0 = f;                   % previous f

while (iter < maxIter) && (norm(step) > tol) 
   step = df\f;           % calculate step size
   v = v0-step;           % step update
   iter = iter+1;         % iter update
   f = g(v);              % f update

   % Damped Newton
   dstep = step;          % initialize damped step size
   diter = 0;             % iteration counter for Damped Newton
   % when damping is needed
   while ((norm(f) > norm(f0)) && (diter < maxDIter))|| (isnan(norm(f))) || (isinf(norm(f)))
       dstep = dstep/2;   % half step size
       v = v0-dstep;      % step update
       f = g(v);          % f update
       diter = diter+1;   % counter + 1
   end

   f0 = f;                % update previous f
   v0 = v;                % update previous v
   df = J(v);             % update J

end
end

