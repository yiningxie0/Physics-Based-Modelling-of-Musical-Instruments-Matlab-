%-------------------------------------------------
% PBMMI Assignment 1 BTB1 Extensions, Guitar body IR and convolution, Musical sequences
% 
% Author: B200247 
% Date: 5/2/23
%
%
% This scrip involves a musical sequence (main script), and 3 functions:
%
% A) Extended KS Algorithm 
% - involves user-selected t60, and automated calculation of loss parameter
% - decay shortening, stretching and tuning, moving pick, pick direction 
%   etc, effects like slurs and glissandis are implemented in the main 
%   function for reasons explained then  
%
% ! Change parameters for y{1} and try in line 114 to see the difference 
%   they make, listening final output may be less ideal as they're afffected 
%   by other operations later
%
% B) function for guitar body impulse response filter with fast convolution
% written by myself 
%
% C) a midi to frequency tranfer function
%
%
%
% The sequence involves for example ('note' here refers to the ones in
% final output):
% 
% i) note 1-3 can be written in loop ,down picksin different strings
% ii) note 4, up pick in the same string with 3
% iii) note 5, a chord
% iv) note 6, mute before next note, common practice for electric guitar
% v) note 8, muting/replacing previous note due to either in purpose or
% pluck in the same string 
% vi) note 11&12 slur (line 151 for clearer sound)
% vii) note 16 glissandi (line 174 for clearer sound)
% Besides, more flexible rhythm can be created by changing position index k
% Rest can be made by changing a 'p' to 1 (easy and smooth) or mannually
% leaving blank (costs less computation)
%
%
% User-friendly issue:
% The displaying is designed with complexity in purpose to just show what
% the function can do, and is prevented from being more clean by Matlab and
% file number limits. If this can be put into other language like C++, with
% parent and child class, a possible UI, and a DAW, it won't require more 
% operations than those necessary for making music with other synthesisers. 
% 
%
% Musically, take it as a guitar beginner just having fun!
%
% references: 1 D. A. Jaffe and J. O. Smith, "Extensions of the Karplus-Strong plucked string algorithm", Com- puter Music Journal 7(2), pp. 56-69, 1983.
%             2 https://ccrma.stanford.edu/~jos/pasp/Extended_Karplus_Strong_Algorithm.html
% sound source of the guitar body response: https://ccrma.stanford.edu/~jiffer8/420/source/piezo.wav
% notes taken from Canon in D (Pachelbel)
%-------------------------------------------------
clear all; clc; close all;  

%//////////////////////////////////////////////////////////////////////////
%% main script - musical sequences from extended KS algorithm
%//////////////////////////////////////////////////////////////////////////

%% basic parameters -------------------------------------------------------
% shared parameters
Fs = 44.1e3;              % sampling rate [Hz]
R = 0.95;                 % dynamics parameter
S = 0.5;                  % stretch factor [0,1]
k = floor(Fs/1.5);        % position index for notes       

% notes in midi - same number in the next line means chords, emited 12
% forms a slur with 11 and will be shown later
%        1  2  3  4  5  5  6  7  8  9  9  9  10 11 12(62) 13 13 13 14 15 16
notes = [50 57 62 66 45 52 49 52 57 47 59 62 54 59        41 57 60 54 57 61]';

% set t60 for each note (though some would be muted in implementation)
% t60 within member of chords should be the same
t = [2 2 2 2 3 3 3 2 2 4 4 4 4 4 5 5 5 4 4 3]';

% pluck direction - 0 are down and non-zero are up
p = [0 0 0 0.4 0 0 0.5 0 0.7 0 0 0 0.4 0 0.8 0 0 0 0 0]';

% pluck position in fraction of the string
mu =[0.5 0.2 0.5 0.4 0.5 0.5 0.5 0.5 0.4 0.4 0.4 0.4 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]';



%% pre-allocate vectors ---------------------------------------------------

y = cell(length(notes),1);     % cell to store outputs of KS function in different lengths

yfinal = zeros(ceil(14*Fs),1); % extimated maximum length for the overall output



%% calculate output for each note -----------------------------------------

% transfer midi to frequency
f = midi2freq(notes);

% /////////////// generate output for each note ///////////////////////////
for i = 1:length(f)
    [z,y0] = KS(Fs,f(i),R,t(i),S,mu(i),p(i));
    
    % y{14} will be further modified for slur, y{20} for glissandi
    if i~=14 && i~=20 
        y{i} = y0;
    else
        y{i} = z;
    end
end

% try here if you want to see the difference of parameters
% soundsc(y{1},Fs) 

% ///////////////////// slurs and glissandis //////////////////////////////
% these 2 are recalculated independently as they require partly different
% parameters, here I change the value of N for slur and make N an gradually
% changing vector for glissandis. 
%
% Actually all effects including unchanged ones can be generalised to a
% vector N, while as most notes doesn't need it to change, calculating
% independently is more efficient
%
% in other software like C++ with parent and child class this can be done
% more efficiently both in running time and in coding

% --------------------- recalculation for y{14} ---------------------------
% recap of parameters, see function session for detailed notes
N = floor(Fs/f(14)-S); 
N1 = floor(Fs/midi2freq(62)-S);
C = (1-N+floor(N))/((1+N-floor(N)));
M = Fs*t(14);
yp1 = 0;
rho = exp(-log(1000)/(f(14)*t(14)))/(((1-S)^2+2*S*(1-S)*cos(2*pi*f(14)/Fs)+S^2)^0.5);

% 1st part in basic frequency
for n = N+1:floor((M-1-N)/6.5)
    yp0 = y{14}(n-N)+C*y{14}(n-N+1)-C*yp1;
    y{14}(n+1) = rho*(S*yp1+(1-S)*yp0);
    yp1 = yp0;
end

% 2nd part in modified frequency
for n = floor((M-1-N)/6.5)+1:M-1
    yp0 = y{14}(n-N1)+C*y{14}(n-N1+1)-C*yp1;
    y{14}(n+1) = rho*(S*yp1+(1-S)*yp0);
    yp1 = yp0;
end

% soundsc(y{14},Fs)


% --------------------- recalculation for y{20} ---------------------------
% recap of parameters
N = floor(Fs/f(20)-S); 
N1 = floor(Fs/midi2freq(73)-S);
M = Fs*t(20);
yp1 = 0;
rho = exp(-log(1000)/(f(20)*t(20)))/(((1-S)^2+2*S*(1-S)*cos(2*pi*f(20)/Fs)+S^2)^0.5);

% make delayline length and tuning coefficient gradually changing vectors
Nv = linspace(N,N1,floor((M-N-1)*8/10));
Nv = [N*ones(1,ceil((M-N-1)*1/10)) Nv N1*ones(1,ceil((M-N-1)*1/10))];
C = (1-Nv+floor(Nv))./((1+Nv-floor(Nv)));

% generate output
for n = N+1:M-1
    yp0 = y{20}(n-floor(Nv(n-N)))+C(n-N)*y{20}(n-floor(Nv(n-N))+1)-C(n-N)*yp1;
    y{20}(n+1) = rho*(S*yp1+(1-S)*yp0);
    yp1 = yp0;
end

%soundsc(y{20},Fs)


%% make chords and sequences ----------------------------------------------

% firt 3 notes in sequence without muting the previous note
for j = 1:3
    yfinal((j-1)*k+1:(j-1)*k+length(y{j})) = yfinal((j-1)*k+1:(j-1)*k+length(y{j}))+y{j};
end

% 4th note, up pick
yfinal(3*k+1:3*k+length(y{4})) = y{4};

% a chord in the 5th position
yfinal(4*k+1:4*k+length(y{5})) = (y{5}+y{6})/2;

% 6th note
yfinal(5*k+1:5*k+length(y{7})) = yfinal(5*k+1:5*k+length(y{7}))+y{7};
% 'mute' part of the 6th note
yfinal(ceil(5.6*k+1):6*k) = 0.15 * yfinal(ceil(5.6*k+1):6*k);

% 7th note
yfinal(6*k+1:6*k+length(y{8})) = yfinal(6*k+1:6*k+length(y{8}))+y{8};

% 8th note muting previous note, that is, e.g when playing on the same 
% string, uppick, or just in purpose, a p>0 compensates for the original vibration
yfinal(7*k+1:7*k+length(y{9})) = y{9};

% 9th note - chord again
yfinal(8*k+1:8*k+length(y{10})) = (y{10}+y{11}+y{12})/3;

% 10th, 'mute' part
yfinal(9*k+1:9*k+length(y{13})) = yfinal(9*k+1:9*k+length(y{13}))+y{13};
yfinal(ceil((9+0.6)*k+1):9*k+length(y{13})) = 0.15 * yfinal(ceil((9+0.6)*k+1):9*k+length(y{13}));

% 11th - 12th the slur
yfinal(10*k+1:10*k+length(y{14})) = y{14};

% 13th, chord
yfinal(12*k+1:12*k+length(y{15})) = (y{15}+y{16}+y{17})/3;

% 14 - 16th with 16th the glissandi
for j = 13:15
    yfinal(j*k+1:j*k+length(y{j+5})) = yfinal(j*k+1:j*k+length(y{j+5}))+y{j+5};
end


%% implement the guitar body impulse response -----------------------------
yfinal = myIR(yfinal);


%% make sound -------------------------------------------------------------
%yfinal = normalize(yfinal);
yfinal = yfinal./max(abs(yfinal));
sound(yfinal, Fs);
audiowrite('KS.wav',yfinal,Fs)

%% functions --------------------------------------------------------------
%//////////////////////////////////////////////////////////////////////////
%% Extended KS algorithm - user-selected t60, decay alternation, moving pick, etc.
%//////////////////////////////////////////////////////////////////////////

function [z,y] = KS(Fs,f0,R,t60,S,mu,p)
% this function works as an extention of the basic KS algorithm, involves
% user-selected t60, decay shortening, stretching and tuning, moving pick,
% pick direction etc, effects like slurs and glissandis are implemented in
% the main function for reasons explained above.
%
    % definition of input parameters
    %
    % Fs - sampling rate [Hz]
    % f0 - desired fundamental frequency
    % R  - dynamics parameter
    % t60 - duration in seconds at f0
    % S - stretch factor [0,1]
    % mu - moving pluck (0,1)
    % p - pick direction [0,1)  0 is down
    %
    % definition of output parameters
    % z - output (the algorithm's input) before feedback loop
    % y - final output after whole loop and normalization
    
    % derived parameters///////////////////////////////////////////////////////
    
    M = t60*Fs;          % duration in samples
    Nexact = Fs/f0-S;    % ideal number of samples in the delay line
    N = floor(Fs/f0-S);  % interger part of the delay at initial frequency
    P = Nexact-N;        % the fractional delay
    C = (1-P)/(1+P);     % allpass filter coefficient
    NE = round(mu*N);    % delay length in samples for moving pick
    
    rho = exp(-log(1000)/(f0*t60))/(((1-S)^2+2*S*(1-S)*cos(2*pi*f0/Fs)+S^2)^0.5); % loss parameter
    
    
    % initialisation///////////////////////////////////////////////////////////
    
    x1 = 0;                    % x[n-1]
    y0 = 0;                    % y'[n] for pick direction filter
    yp1 = 0;                   % y''[n-1] for tuning filter
    
    v = -1+2.*rand(N+1,1);     % white noise input
    y = zeros(M,1);            % output
    z = zeros(M,1);            % final average output
    dlinebufE = zeros(NE+1,1); % delay buffer for moving pick, include the situation NE=0
    
    
   
    % filters before delay circle//////////////////////////////////////////////
    for n = 0:N 
       
        % --------------------- dynamics filter -------------------------------
        % lowpass filter when R belongs to (0,1)
        % R = 0 means no filter and 1 means no sound 
        % x[n]=(1-R)*v[n]+R*x[n-1]
    
        x0 = (1-R)*v(n+1)+R*x1; 
        y(n+1) = x0;
        x1 = x0; 
    
        % --------------------- pick-position filter --------------------------
        % comb filter
        % y[n]=x[n]-x[n-mu*N]
    
        dlinebufE(NE+1) = x0;
        dlinebufE = circshift(dlinebufE,1);
        xe = x0-dlinebufE(NE+1);
    
        % ---------------------- pick-direction filter ------------------------
        % lowpass filter
        % p = 0 for one direction, between (1,0) for the other
        % y[n] = (1-p)*x[n]+p*y[n-1]
    
        y(n+1) = (1-p)*xe+p*y0;
        y0 = y(n+1);
    end
    
    z = y;   % interrupt here for more modifications when needed
    
    
    % main Karplus-Strong algorithm/////////////////////////////////////////////
    %
    % with decay-time alternation
    % includes an allpass tuning filter, a 2-point average filter and a delayline
    % for tuning filter y[n] = C*x[n]+x[n-1]-C*y[n-1]
    % for KS, n>N, y[n]=rho*(S*y[n-N-1]+(1-S)*y[n-N])
    
    for n = N+1:M-1
    
        yp0 = y(n-N)+C*y(n-N+1)-C*yp1;
        y(n+1) = rho*(S*yp1+(1-S)*yp0);
        yp1 = yp0;
        
    end
     

    % plot
    % t = (0:(M-1))/M*t60;
    % 
    % subplot(2,1,1)
    % plot(t,y);
    % xlabel('Time(s)');                         
    % ylabel('Amplitude');             
    % title('Output signal');   
    % 
    % 
    % subplot(2,1,2)
    % 
    % NFFT= 2^nextpow2(length(y));  
    % Y = fft(y,NFFT);
    % f = (0:NFFT/2-1)/NFFT*Fs;
    % Ydb = 20*log10(abs(Y(1:NFFT/2,:)));
    % 
    % plot(f,Ydb);
    % hold on
    % xline(f0,'--');
    % legend('Basic Karplus-Strong Output','f=f0')
    % xlabel('Frequency(Hz)');                         
    % ylabel('dB');             
    % title('Output signal');   
end

%//////////////////////////////////////////////////////////////////////////
%% function for guitar body impulse response filter -----------------------
%//////////////////////////////////////////////////////////////////////////
% This function adds a guitar body impluse response to the input and
% outputs the IR filtered signal 
%
function y = myIR(x)
    [ir,fs]=audioread('piezo.wav');      % read IR file
    ir = 0.5*sum(ir,2);                  % put to mono in case of stereo
    MFFT = max(length(x),length(ir));    % decide fft length
    y = ifft(fft(x,MFFT).*fft(ir,MFFT)); % convolve for output
    % sound source : https://ccrma.stanford.edu/~jiffer8/420/source/piezo.wav
end


%//////////////////////////////////////////////////////////////////////////
%% function to transfer midi to frequency ---------------------------------
%//////////////////////////////////////////////////////////////////////////
% This function transfers midi notes to corresponding frequency
%
function y = midi2freq(x)
    y = 440*2.^((x-69)/12);
end

