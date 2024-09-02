% ------------------------------------------------
% PBMMI Assignment 5 BTB1 Prune, Generalization and Customization
% 
% Author: B200247 
% Date: 1/4/23
%
% This script may be viewed as a backend of a GUI though GUI is not
% written. It allows users to:
%
% a) Customize audio input
% b) Pick custom materials for the plate 
% c) Customize T60 for anchor frequencies (both T60 and anchor frequencies
%    can be changed), correponding loss parameters are calculated 
%    - so there is EQ
% d) Customize input and output positions with no limit in numbers or
%    output channels
% e) you may try to decide how much to prune the modal frequencies by 
%    yourself too by changeing the 'cents' in first section 
%
% One important feature is the MODAL PRUNING that accelerates the speed,
% and this generally includes pruning according to:
%
% (starting from line 171, 202, 225
% for impulse force, mono input and mono output, the final modes after 
% prunning in current 'cents' setting is around 3800, and around 0.5s in
% the loop
%
% 1) output (Phi_o): prune modes that have least contributions to outputs
%
% 2) input (vector/matrix D): prune the modes whose sum of modal shape
%    functions in all input positions can still be counted as having 
%    negligible contribution.  
%
% 3) frequency (omega): grouping similar frequencies and pruning them to
%    take only one place. It did took me much time to think of and find the
%    accumarray function but the result is gladly compact. 
%
% Details are explained in the main script. 
% 
% ALL THE ABOVE ARE GENERALIZED TO SUPPORT USER DEFINED NUMBER OF CHANNELS
% 
% Plots spectrogram and the loss parameter according to modal frequencies.
% Note that the x axis frequencies are scaled in different ways for the
% plots - I used just default setting for spectrogram.
%
% most brief instruction: default setting is audio input, 
%                         you may change audio_input to 0 for impulse.
% Have fun!
% 
% ------------------------------------------------
clear all; clc; close all;  

%% customized parameters
% input: 1: read in audio, 0/else: impulse force
audio_input = 0;      

% audio input - tranfer all input into mono audio, as I assume this reverb
% is applied to only one channel one track at a time
if audio_input == 1
    [in,SR] = audioread('Guitar_dry.wav');      % read and store audio file, Fs: sample rate [Hz] x: sample   
    [Nf,nc] = size(in);                         % get length and number of channels of input    
    if nc ~=1
        in = sum(in,2)/nc;                      % Average channels to mono when the input audio is not mono
    end
end

% sample rate if not taken from audio
if audio_input ~= 1
    SR = 44.1e3;
end

% material: 1: steel, 2: glass, 3: stone
material = 3;

% T60 time of the corresponding modal frequency - please don't change the
% first and last frequency value.
T60 = [3 4.5 6 1.5 1];
freq = [1, 110, 440, 1760, 2*SR];

% input and output positions [x1,y1;...;xn,yn]'
pin = [1.0 0.5; 0.1 0.5]';      % input positions 
pout = [0.2 0.2; 0.5 0.9]';     % output positions 

cents = 0.15;                   % reference cents for pruning in modal frequencies
% which differs from person to person and this does not follow the same
% convention as perceptual difference in cents for pure sounds. When it is
% 0.1 it sounds good to me, 0.15 sounds and looks same, when it's 0.2 I
% could hear the difference. And though there is of course trade off issue
% - the higher the value is, the faster this code runs, but could be less
% authentic.


%% input parameters

% physical parameters
Lx = 2;                           % x dimensions [m]
Ly = 1;                           % y dimensions [m]
H = 5e-4;                         % thickness [m]
T = 700;                          % tension [N/m]

switch material
    case 1 %steel
        rho = 7.87e3;             % density [kg/m^3]
        E = 200e9;                % youngs modulus [N/m^2]
        v = 0.29;                 % poisons ratio
    case 2 % glass
        rho = 2.50e3;             % density [kg/m^3]
        E = 80e9;                 % youngs modulus [N/m^2]
        v = 0.25;                 % poisons ratio
    case 3 % stone
        rho = 2900;               % density [kg/m^3]
        E = 70e9;                 % youngs modulus [N/m^2]
        v = 0.25;                 % poisons ratio
end


%% derived parameters

% derived paramaters
k = 1/SR;                         % sample period
K = sqrt(E*H^2/(12*rho*(1-v^2))); % kappa
c = sqrt(T/(rho*H));              % wavespeed

if audio_input ~= 1
    Nf = ceil(SR*max(T60));        % number of samples
end

% maximums for stability --------------------------------------------------

% maximum modal frequency
w_max = 2/k;   

% maximum modal wavenumber
beta_max = sqrt((-c^2+sqrt(c^4+4*w_max^2*K^2)/(2*K^2)));  
% maximum mode shape index
Mx = floor(sqrt(beta_max^2-pi^2/Ly^2)*Lx/pi);
My = floor(sqrt(beta_max^2-pi^2/Lx^2)*Ly/pi);

% mode shapes, wavenumbers, frequencies -----------------------------------

[mx,my]=meshgrid(1:Mx,1:My);                    % all possible index
beta = sqrt(pi^2/Lx^2*mx.^2+pi^2/Ly^2*my.^2);   % all modal wavenumbers

% find index pairs and wavenumbers that meets stability constraints
beta(beta>=beta_max) = 0; 
[mx,my,beta] = find(beta);

% modal frequencies
omega = sqrt(K^2*beta.^4+c^2*beta.^2);

% sort for interpolation and keep the index
[omega, idx0] = sort(omega);

% modal damping coefficients ----------------------------------------------
T60_in = interp1(log2(freq), T60, log2(omega/(2*pi)),'pchip');
sigma = 6*log(10)./T60_in;

%sigma = sigma0+sigma1*beta.^2;

% make sure damping coefficients are larger than 0
assert(min(sigma)>0,'error');

% modal shape functions  --------------------------------------------------
Phi_i = 2/sqrt(Lx*Ly)*sin(mx*pi*pin(1,:)/Lx).*sin(my*pi*pin(2,:)/Ly);
Phi_o = 2/sqrt(Lx*Ly)*sin(mx*pi*pout(1,:)/Lx).*sin(my*pi*pout(2,:)/Ly);

% sort
Phi_i = Phi_i(idx0,:);
Phi_o = Phi_o(idx0,:);

%% mode prune 1 - prune output

% prune modes that have least contributions to outputs - either approximate
% to 0 or too small to be percepted. As to the ratio, the number of Phi_os
% that have values less than max*exp(-4) (or even exp(-2)) turns out to be
% the same with those small than max*exp(-12), and I can't hear/see the difference.
% Hence it shall be safe to prune up to this.

% considered if the output is multichannel, in this case a mode is pruned
% only when it can be counted as having negligible contribution in all channels

% find index of modes that can be pruned
idx = find(all(abs(Phi_o)<max(abs(Phi_o))*exp(-4),2));

% apply pruning to all necessaray and relevant parameters
sigma(idx)=[];
omega(idx)=[];
Phi_i(idx,:)=[];
Phi_o(idx,:)=[];

%% parameters for state update

B = (2-k^2*omega.^2)./(1+k*sigma);
C = -(1-k*sigma)./(1+k*sigma);

% when multipul input positions are selected, the effect that is exerted to
% the plate is their sum. Under the assumption that the force (input audio)
% is the same, it can be equivalent to summation in D before the loop
D = sum(Phi_i./(1+k*sigma),2);


%% mode prune 2 - prune input

% prune modes that have least contributions to inputs, also considering
% the possibility of having multiple inputs, we prune only the modes whose 
% sum of modal shape functions in all input positions can still be counted
% as having negligible contribution. This can be done by checking sum of
% Phi_i, while checking D is equivalent.

% same as output, a threshold of max*exp(-4) prunes the same modes as
% max*exp(-12), and it shall be safe to do so.

% find index of modes that can be pruned
idx2 = find(abs(D)<max(abs(D))*exp(-4));

% apply pruning to all necessaray and relevant parameters
omega(idx2)=[];
sigma(idx2)=[];
Phi_o(idx2,:)=[];
B(idx2)=[];
C(idx2)=[];
D(idx2)=[];


%% mode prune 3 - similar frequency
% here I select model frequency as reference - which should be the same
% with refering to beta or sigma - as can be seen from the formula, same 
% omega should mean same beta when beta >= 0, and also same sigma.

% the idea is to take modes that have similar frequencies as the same one,
% pruning them to only one place (with their mean) in B, C, D and
% compensating for the weights in output shape function by summing their
% Phi_os up.

% the approach is to 
% i) first label the omegas into different groups by 1) taking
% log 2) with the division and floor function, putting frequencies that have
% a difference within the the specified 'cents' into the same group with
% 'ref' as its label.
% ii) prune modes that have omegas in the same group by keeping only their
% mean in B, C, D and taking their sum in output modal shape function. As
% far as we have the correct label, this can be done directly with
% accumarray - though I did take much time to find this function.

% One issue here is similar omega doesn't necessarily mean similar Phi, but
% in real practice I found that this difference could be perceptually small 
% enough to be neglected too - or as far as the 'cents' chosen are small
% enough and in that case we are likely pruning similar omegas with also
% similar mxs and mys.

% Create mask according to similarity in omega
ref = floor(log2(omega/(2*pi))*1200/cents);

% make its index starting from 1 to have miminum vector size in process
ref = ref-min(ref)+1;

% omega has already been sorted to calculate sigma before (so as other
% relevant parameters), and this makes things easier to explain. The label
% should follow an order and be something of the form [1;2;5;7;7;9...] etc,
% with accumarray(ref, a, [], @mean), the 4th and 5th value of a is taken
% the mean and stored in the 7th place of new a... 
% However, even if sigma is calculated in another way (original way for
% example), and parameters are not sorted, this pruning still works.

% Do the accumarray - this prunes similar modes
% a is also a reference mask for later pruning - as accumarrary function
% leaves 0s for index that does not have a corresponding value and we want
% to remove these 0s later. This a makes sure only this kind of 0s are
% removed, while in B and D there might be values that are themselves 0.
a = ones(length(B),1);
a = accumarray(ref, a, [], @mean);
B = accumarray(ref, B, [], @mean);
C = accumarray(ref, C, [], @mean);
D = accumarray(ref, D, [], @mean);

% operation in output - accumarrary does not deal with data more than 1 
% dimension, so loop in columns is necessary when user chooses multiple 
% output points.
% the loop may be replaced by cell but may not be more efficient
% when it is mono, it is exactly the same with simply
% phi_o = accumarray(ref, phi_o, [], @sum);
Phi_oo = zeros(max(ref),width(pout));
for i = 1:width(pout)
    phi_o = Phi_o(:,i);
    phi_o = accumarray(ref, phi_o, [], @sum);
    Phi_oo(:,i) = phi_o;
end

% This prunes redundant 0s
idx3 = find(a==0);
B(idx3)=[];
C(idx3)=[];
D(idx3)=[];
Phi_oo(idx3,:)=[];


%% initialization

if audio_input ~= 1
    in = [1; zeros(Nf-1,1)];     % input vector (impulse)
end
out = zeros(Nf,width(pout));     % output vector
 
% states
p0 = zeros(length(B),1);         
p1 = p0;
p2 = p1;

%% main loop
tic
for n = 1:Nf

    % update states
    p0 = B.*p1+C.*p2+D.*in(n);
    p2 = p1;
    p1 = p0;

    % sum and add to output
    % out(n,:) = sum(Phi_oo.*p0);
    % for stereo and multichannel now it is the following that runs faster
    out(n,:) = p0'*Phi_oo;
end
toc

% normalize and differentiate output
v = diff(normalize(out));

%% make sound
sound(v,SR);
v = v/max(abs(v));
audiowrite('plate.wav',v,SR)
%% draw spectrogram
figure
spectrogram(sum(v,2))
title('Spectrogram for output')
figure
plot(log(omega/(2*pi)),sigma);
title('Loss parameter according to modal frequencies')
xlabel('modal frequencies(in log)')
ylabel('sigma')