%-----------------------------------------------------
% PBMMI Assignment 4 BTB1 
% ----- Reset parameters, frequency-dependent loss, variation in T60 and guitar simulation
%
% (took me around 10s to run, but with reasons)
%
% Author: B200247 
% Date: 18/3/23
%
% This code made the following changes to the original code in function 
% starting from line 155:
% 
% A) Reset input parameters in terms of fundamental frequency and
%    inharmonicity, hence more perceptually reasonable
%
% B) Introduced frequency dependent loss with changes in:
%    i) adjusted matrices B and C
%    ii) adjusted numerical stability condition and recalculated 'hmin'
%
% C) Specified a variation in T60 - with the 'loss' vector
% (and implfied density to linear density in function to reduce 1
% parameter...)
%
% And made a guitar simulation with steel light strings:
% - calculated inharmonicity from other physical parameters
% - tuned tension based on fundamental frequencies
% - added some randoness to input force and location 
% - outputted plucks of 6 strings (itype1) in panning channels, and a down
%   strum (itype2) with a slight delay in each string - hence more like real situation 
% - plotted transforms of each string's output
% and as in KS, added a guitar body IR filter.
%
%
% Explanation for possible concerns:
% > for efficiency, as also explained in main script, the 220200x6 y matrix
%   could actually be omitted, while here I included it in just for plot. 
% > for simulation time, Tf had better be larger than T60 of low E to
%   have full output, while this results in very long running time. I still
%   prioritize more accurate plotting here, otherwise Tf could be made shorter.
% > for consolidating all 6 strings, it doesn't look possible in the 
%   current setting as B and C also differ with strings, and as for now 6
%   strings are independent from each other and for writing music, one's
%   enough to deal with all different frequencies, timbres and play chords (by summing) etc..   
%
% reference (for inharmonicity formula and new hmin with frequency-dependent loss respectively:)
% 1 Inharmonicity in plucked guitar strings, Chris J. Murray and Scott B. Whitfield (2022)
% 2 Numerical Modeling of String/Barrier Collisions: The Fretboard. Bilbao, Stefan & Torin, Alberto. (2014). 
%-----------------------------------------------------

clc
clear all
close all

%% set parameters
% steel light string, diameters from https://www.elixirstrings.com/support/string-tension-for-tuning-guitar
% High E (1st) - 0.012"
% B (2nd) - 0.016"
% G (3rd) - 0.024"
% D (4th) - 0.032"
% A (5th) - 0.042" 
% Low E (6th) - 0.053"
% -----------------------------------------------------------------

% flags -----------------------------------------------------
plot_on = 0;                % in-loop plotting on (1) or off (0)
itype = 1;                  % pluck
itype2 = 2;                 % strike

% I/O -------------------------------------------------------
SR = 44100;                 % sample rate (Hz)
Tf = 8.5;                   % duration of simulation (s)
dur = 0.001;                % duration of excitation (s)
exc_st = 0.001;             % start time of excitation (s)
xo = 0.5;                   % coordinate of output (normalised, 0-1)

% coordinate of excitation (normalised, 0-1)
xi = [0.55+0.02*rand(1) 0.50+0.02*rand(1) 0.45+0.02*rand(1) 0.40+0.02*rand(1) 0.35+0.02*rand(1) 0.30+0.02*rand(1)];               

% physical parameters ---------------------------------------
L = 0.65;                          % length (m)
E = 2e11;                          % Young's modulus (Pa) for steel
rho0 = 7850;                       % density (kg/m^3)
loss = [100, 8; 1000, 6];          % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
% this shall differ with strings while I am being lazy here...

% desired string radius (m)
r = [6.731e-4 5.334e-4 4.064e-4 3.048e-4 2.032e-4 1.524e-4];

% desired fundamental frequencies (Hz)
% - here I set frequency first to calculate the needed tension instead of
% the opposite order, for that in real guitar tension is something to be
% adjusted (tuned) to meet frequency requirement too
f0 = [82.41 110.00 146.83 196.00 246.94 329.63]; 

% peak amplitude of excitation (N)
famp = [4+0.5*rand(1) 3.5+0.5*rand(1) 3+0.5*rand(1) 2.5+0.5*rand(1) 2+0.5*rand(1) 1.5+0.5*rand(1)];

% specifical parameters for the current output ---------------
t_total = 8;                       % total output time for out1 (s)
gap = 0.5;                         % gap between notes (s)


%% derived parameters and vectors

rho = 7850*pi*r.^2;                % linear density (kg/m)
T = 7850*(2*L*f0).^2*pi.*r.^2;     % tension (N)
B = pi^3*E*r.^4./(4*T*L^2);        % inharmonicity parameter

Nf = floor(SR*Tf);                 % number of time steps


%% initialise outputs

y = zeros(Nf,6);                   % vector to store each note's output 
% - for plot purpose, otherwise directly storing them into desired 
% locations of out1 might be more efficient

out1 = zeros(floor(SR*t_total),2); % output for separate plucks
out2 = zeros(Nf,2);                % output for a strum

%% main

for i = 1:6

    % calculate each pluck and write output with a panning channel choice
    y(:,i) = makestring(plot_on,itype,rho(i),loss,f0(i),B(i),SR,Tf,xi,famp(i),dur,exc_st,xo);
    out1((i-1)*floor(gap*SR)+1:(i-1)*floor(gap*SR)+Nf, mod(i+1,2)+1) = y(:,i);

    % calculate a strum with slight delay (exc_st+0.025*i) in each string
    out2 = out2+makestring(plot_on,itype2,rho(i),loss,f0(i),B(i),SR,Tf,xi,famp(i),dur,exc_st+0.025*i,xo);
end

% combine output
out = out1;
out(7*floor(gap*SR)+1:7*floor(gap*SR)+length(out2),:) = out2;

% add IR
out = myIR(out);


%% make sound
out = out./max(abs(out));
sound(out,SR);
audiowrite('FDTDstring.wav',out,SR)
%% plot
yfft = 10*log10(abs(fft(y)));
plot([0:Nf-1]'*(SR/Nf), yfft)
xlim([0 SR/2])
legend('E2','A2','D3','G3','B3','E4')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Transform of output')


%% function for stiff string with frequency independent and dependent loss
function y = makestring(plot_on,itype,rho,loss,f0,B,SR,Tf,xi,famp,dur,exc_st,xo)
%
% this function inputs the following parameters:
%
% plot_on: in-loop plotting on (1) or off (0)
% itype: type of input: 1: pluck, 2: strike
% rho: linear density (kg/m)
% loss: loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
% f0: fundamental frequencies(Hz)
% B: inharmonicity parameter
% SR: sampling rate
% Tf: duration of simulation (s)
% xi: coordinate of excitation (normalised, 0-1)
% famp: peak amplitude of excitation (N)
% dur = 0.001;                
% exc_st = 0.001;            
% xo = 0.5;
%
% and output a plucked/stroke string simulation


% derived parameters ------------------------------------------------------

gamma = 2*f0;               % wave speed
K = sqrt(B)*(gamma/pi);     % stiffness constant 
    
%loss parameters
zeta1 = (-gamma^2+sqrt(gamma^4+4*K^2*(2*pi*loss(1,1))^2))/(2*K^2);
zeta2 = (-gamma^2+sqrt(gamma^4+4*K^2*(2*pi*loss(2,1))^2))/(2*K^2);
sig = 6*log(10)*(-zeta2/loss(1,2)+zeta1/loss(2,2))/(zeta1-zeta2);   % frequency independent loss parameter
sig1 = 6*log(10)*(1/loss(1,2)-1/loss(2,2))/(zeta1-zeta2);           % frequency dependent loss parameter

k = 1/SR;                   % time step
Nf = floor(SR*Tf);          % number of time steps

% grid spacing ------------------------------------------------------------
hmin = sqrt(0.5*(gamma^2*k^2+4*sig1*k+sqrt((gamma^2*k^2+4*sig1*k)^2+16*K^2*k^2)));  % minimal grid spacing
N = floor(1/hmin);          % number of grid spacings
h = 1/N;                    % reset grid spacing

% input -------------------------------------------------------------------
t = [0:Nf-1]'/SR;           % time vector
ind = sign(max(-(t-exc_st).*(t-exc_st-dur),0));     % select time period
f = 0.5*famp*ind.*(1-cos(itype*pi*(t-exc_st)/dur)); % calculate force

% matrices for update -----------------------------------------------------
% u = B*u1-C*u2+J*f(n)

% input/output position in grid point index
xi_m = floor(N*xi);
xo_m = floor(N*xo);     

% matrices of spatial derivatives
e = ones(N-1,1);
Dxx = spdiags([e -2*e e], -1:1, N-1,N-1)/h^2;
Dxxxx = spdiags([e -4*e 6*e -4*e e], -2:2, N-1, N-1)/h^4;
Dxxxx(1,1) = 5/h^4; Dxxxx(N-1, N-1)=5/h^4;

% matrices for state update
B = 1/(1+sig*k)*(2*speye(N-1)+gamma^2*k^2*Dxx-K^2*k^2*Dxxxx+2*sig1*k*Dxx);
C = 1/(1+sig*k)*((-1+sig*k)*speye(N-1)-2*sig1*k*Dxx);

% input force vector
J = sparse(N-1,1);
J(xi_m) = 1/(1+sig*k)*k^2/(rho*h);

% output vector
c = sparse(1,N-1);              
c(xo_m) = 1;

% initialise scheme variables ---------------------------------------------

u2 = zeros(N-1,1);              % state
u1 = u2;                        % state
u = u2;                         % state

y = zeros(Nf,1);                % output
xax = (1:N-1)'*h;               % x-axis for plotting

% main loop ---------------------------------------------------------------

for n=1:Nf
    
    % update state, and insert current value of f.
    u = B*u1 + C*u2 + J*f(n);
    
    % read output
    y(n) = c*u;
    
    % plot
    if(plot_on==1)
        % draw current state
        figure(1)
        plot(xax, u, 'k');
        xlabel('x (m)')
        ylabel('u (m)')
        axis([0 1 -0.005 0.005])
        drawnow
    end
    
    % shift state
    u2 = u1;
    u1 = u;  
end

end

%% function for guitar body impulse response filter -----------------------
% This function adds a guitar body impluse response to the input and
% outputs the IR filtered signal 

function y = myIR(x)
    [ir,fs]=audioread('piezo.wav');      % read IR file
    ir = 0.5*sum(ir,2);                  % put to mono in case of stereo
    MFFT = max(length(x),length(ir));    % decide fft length
    y = ifft(fft(x,MFFT).*fft(ir,MFFT)); % convolve for output
    % sound source : https://ccrma.stanford.edu/~jiffer8/420/source/piezo.wav
end