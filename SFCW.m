
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Entry Here
fc = 24e9;      % 24 Ghz Wave
c = 3e8;        %Speed of light
lambda = c/fc;  %wavelength
Nsweep = 1;     %Perform N sweeps of the radar and overlap the plots

BW = 2e9        %2Ghz BW
Fc = 4e6        %Minimum Freq ex 2Mhz so we have 1000 steps

Phase_NoiseAndOffset    = [-80,100e3] %Noise and Offset 
SystemWhite_Noise       = [-60]       %Iq Noise floor (not yet used)
Circulator_Issolation   = [-20];      %-20dB issolation 

car_dist = 2;                         %Distance to crude_oil                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FreqSteps = BW/Fc 
fs  = BW*2;
smshop = 100%samples per hop 
L = length(linspace(0,smshop/Fc,10001));
wave = zeros(L,FreqSteps);
figure(1)

pnoise = comm.PhaseNoise('Level',[-80],'FrequencyOffset',[100e3], ...
    'SampleRate',2*100e3);

for steps = 1:FreqSteps
   t = [0:1:10000];
   I  = cos(((2*pi*(Fc*steps/(2*BW)*t))));
   Q  = sin(((2*pi*(Fc*steps/(2*BW)*t))));
   Z  = I - 1i*Q;
   %wave = vertcat(wave,Z');
   Z_Phase = phase_noise(Z',Phase_NoiseAndOffset(1),Phase_NoiseAndOffset(2));
   wave(:,steps) = Z_Phase;
   %wave(:,steps) = Z';
   %plot(I)
   %Frq = Fc*steps;
   steps
   %plot(real(wave(:,steps)));
end



%% Target Model
% The target of an ACC radar is usually a car in front of it. This example
% assumes the target car is moving 50 m ahead of the car with the
% radar, at a speed of 96 km/h along the x-axis. 
%
% The radar cross section of a car, according to [1], can be computed based
% on the distance between the radar and the target car.

car_rcs = db2pow(min(10*log10(car_dist)+5,20));

cartarget = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',car_rcs,'PropagationSpeed',c,...
    'OperatingFrequency',fc);

%cartarget = phased.RadarTarget('Model','Swerling2','MeanRCS',car_rcs,'PropagationSpeed',c,...
%   'OperatingFrequency',fc);

carmotion = phased.Platform('InitialPosition',[car_dist;0;0.0]);

%%
% The propagation model is assumed to be free space.

channel = phased.WidebandFreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

%% Radar System Setup
% The rest of the radar system includes the transmitter, the receiver, and
% the antenna. This example uses the parameters presented in [1]. Note that
% this example models only main components and omits the effect from other
% components, such as coupler and mixer. In addition, for the sake of
% simplicity, the antenna is assumed to be isotropic and the gain of the
% antenna is included in the transmitter and the receiver.

ant_aperture = 6.06e-4;                         % in square meter
ant_gain = aperture2gain(ant_aperture,lambda);  % in dB

tx_ppower = db2pow(1)*1e-3;                     % in watts
tx_gain = 9+ant_gain;                           % in dB

rx_gain = 15+ant_gain;                          % in dB
rx_nf = 4.5;                                    % in dB



transmitter = phased.Transmitter('PeakPower',tx_ppower,'Gain',tx_gain);
receiver = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,...
    'SampleRate',fs);
radarmotion = phased.Platform('InitialPosition',[0;0;0]);


%% Actual Sweep

sig = combineSteps(wave,FreqSteps);
plotSweepSpectrum(sig,fs);

tgt_pos = [1;0;0];% x y z position

for m = 1:Nsweep
    disp('Sweeping')
    
    %% Setup the TX signal
    txsig =  transmitter(sig);
    
    %% Propagate the signal and reflect off the target

    [radar_pos,radar_vel] = radarmotion(1);
    txsig = channel(txsig,radar_pos,tgt_pos,radar_vel,radar_vel);
 
    txsig = cartarget(txsig); %step(H,X,UPDATE)
    
    %% Dechirp the received radar return
    txsig = receiver(txsig);    
    
    txsig = circulator(Circulator_Issolation,sig,txsig);
    dechirpsig       = dechirp(txsig,sig);
    dechirpsig       = IQ_filter(dechirpsig); %Fillter the data through a LPF
 
    FFT_range(c,fs,dechirpsig,0)
   

end
    
 clear;

 
    %% FilterDesigner LPF for filter IQ jumps

 function [filtered_data] = IQ_filter(IQ_data)
% All frequency values are in Hz.
Fs = 4000000000;  % Sampling Frequency

N  = 25;      % Order
Fc = 350000;  % Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass('N,F3dB', N, Fc, Fs);
Hd = design(h, 'butter');
output = filter(Hd,IQ_data);
filtered_data = output;
figure(4)
hold on
%plot(real(IQ_data))
%plot(real(output))
%plot(imag(output))
end


 
 function FFT_range (c,Fs,IQ_data,sweep_slope)
    %FFT
    IQ_data = decimate(IQ_data,1000);
    Fs = Fs/1000;
    T = 1/Fs;             % Sampling period
    L = length(IQ_data);  % Length of signal
    t = (0:L-1)*T;        % Time vector
    window = hann(L);
    Y = fft(IQ_data.*window);
    P2 = abs(Y/L);
    P1 = P2(1:floor((L/2+1)));
    P1(2:end-1) = 2*P1(2:end-1);

    L = length(IQ_data);  % Length of signal
    f = Fs*(0:(L/2))/L;
    deltaR = c ./ (2*f*1000);
    deltaR = linspace(0,99,L/2+1);  %  37.5m
    deltaR = deltaR ./2;

    figure(3)
    hold on
    axis([0 10 -120 0])
    %axis([0 6e5 -160 10])
    plot(deltaR,mag2db(P1))
    %title('Range Power Plot')
    xlabel('Range (m)')
    ylabel('|P1 db(m)|')
 
    %Est Range
    [y,x] = max(mag2db(P1)); % find peak FFT point
    rng(x); % map to range array
 end


 function [combined] = combineSteps(wave,steps)
 disp('Combining Waveforms(this may take a bit)')
 wholesig = [] ;  
     for count = 1:steps
     sig =  wave(:,count);
     wholesig = vertcat(wholesig,sig);
     end
combined = wholesig;
 end
 
 
 function plotSweepSpectrum(data,fs)
  spectrogram(data,32,16,32,fs,'yaxis');
 end
  
 function [res] = range_resolution(c,steps,BW,fs)
 L = length(dechirpsig);  % Length of signal
 f = fs*(0:(L/2))/L;
 deltaR = c ./ (2*f);
 end

 function [IQ_Data_noise] = phase_noise(IQ_Data,PhaseNoise,Offset)
 pnoise = comm.PhaseNoise('Level',PhaseNoise,'FrequencyOffset',Offset, ...
     'SampleRate',2*Offset);
 IQ_Data_noise = pnoise(IQ_Data);
 WhiteNoise    = awgn(IQ_Data_noise,1,.01);
 IQ_Data_noise = WhiteNoise;
 end

 function [txsig_out] = circulator(coupling_factor, initial, target)
    coupling_factor = 10^(coupling_factor/10);
    txsig_out = target + coupling_factor * initial;
 end



