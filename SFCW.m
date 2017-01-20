

fc = 24e9; % 24 Ghz Wave
c = 3e8;   %Speed of light
lambda = c/fc;  %wavelength

%%
% The sweep time can be computed based on the time needed for the signal to
% travel the unambiguous maximum range. In general, for an FMCW radar
% system, the sweep time should be at least 5 to 6 times the round trip
% time. This example uses a factor of 5.5.

range_max = 4;
%tm = 5*range2time(range_max,c);
tm = 0.00001; %Use sweep of 1ms
%%
% The sweep bandwidth can be determined according to the range resolution
% and the sweep slope is calculated using both sweep bandwidth and sweep
% time.


smshop = 10%samples per hop 
BW = 2e9
Fc = 1e7
fs  = BW*2;

FreqSteps = BW/Fc
L = length(linspace(0,smshop/Fc,1001));
wave = zeros(L,FreqSteps);
figure(1)
for steps = 1:FreqSteps
   t = [0:1:1000];
   I  = cos(((2*pi*(Fc*steps/(2*BW)*t))));
   Q  = sin(((2*pi*(Fc*steps/(2*BW)*t))));
   Z  = I - 1i*Q;
   %wave = vertcat(wave,Z');
   wave(:,steps) = Z';
   %plot(I)
   %Frq = Fc*steps;
   steps
   plot(real(wave(:,steps)));
end
t = [0:1:1000];
%I  = cos(((2*pi*(Fc*steps/(2*BW)*t))));
%plot(t*1/BW*1000,I);

%wind  = hann(length(wave))';
%wave = wave.*wind';
%spectrogram(wave,32,16,32,2e3,'yaxis');
%plot(real(wave))

%plot(t,I)
%hold on
%plot(t,Q)



sig = wave;


%% Target Model
% The target of an ACC radar is usually a car in front of it. This example
% assumes the target car is moving 50 m ahead of the car with the
% radar, at a speed of 96 km/h along the x-axis. 
%
% The radar cross section of a car, according to [1], can be computed based
% on the distance between the radar and the target car.

car_dist = 2;
car_rcs = db2pow(min(10*log10(car_dist)+5,20));

cartarget = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',car_rcs,'PropagationSpeed',c,...
    'OperatingFrequency',fc);

%cartarget = phased.RadarTarget('Model','Swerling2','MeanRCS',car_rcs,'PropagationSpeed',c,...
%   'OperatingFrequency',fc);

carmotion = phased.Platform('InitialPosition',[car_dist;0;0.0]);

%%
% The propagation model is assumed to be free space.

channel = phased.FreeSpace('PropagationSpeed',c,...
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


%rng(2012);
Nsweep = 4;
%xr = complex(zeros(waveform.SampleRate*waveform.SweepTime,Nsweep));

for m = 1:Nsweep

    dechirp_Data = [];
    dechirp_Data2 = [];

    for step_hop = 1:steps
        
    sig =  wave(:,step_hop);
    txsig = transmitter(sig);
    
    % Propagate the signal and reflect off the target
    tgt_pos = [20;0;0];% x y z position
    [radar_pos,radar_vel] = radarmotion(1);
    txsig = channel(txsig,radar_pos,tgt_pos,radar_vel,radar_vel);
 
    txsig = cartarget(txsig); %step(H,X,UPDATE)
    
    %Dechirp the received radar return
    txsig = receiver(txsig);    
    

    dechirpsig       = dechirp(txsig,sig);
   % plot(real(dechirpsig))
    dechirp_Data     = vertcat(dechirp_Data,dechirpsig(600:1000)); %BRKPNT
    dechirp_Data2(:,step_hop) = dechirpsig;
    
    figure(12)
    hold on
    plot(real(dechirp_Data))
    
    figure(2)
    hold on
    plot(real(txsig));
    plot(real(sig));
    figure(1)
    hold on
    plot(real(dechirpsig));
    I_avg(step_hop)= mean(real(dechirpsig(600:1000)));
    Q_avg(step_hop)= mean(imag(dechirpsig(600:1000)));
    step_hop
    end
    figure(13)
    plot(I_avg)
    hold on
    plot(Q_avg)
    
    %dechirpsig = decimate(dechirpsig,10000);
    %FFT
    dechirp_Data = I_avg + 1i*Q_avg ;
    plot(real(dechirp_Data))
    Fs = fs;              % Sampling frequency
    T = 1/Fs;             % Sampling period
    L = length(dechirp_Data);             % Length of signal
    t = (0:L-1)*T;        % Time vector
    window = hann(L)';
    Y = ifft(dechirp_Data.*window);
    f = Fs*(0:(L/2))/L;
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    
    %rng = c*f/sweep_slope/2; %  RANGE PLOT IN M 
    figure(3)
    hold on
    %axis([0 6 -120 10])
    plot(f,mag2db(P1))
    title('Range Power Plot')
    xlabel('Range (m)')
    ylabel('|P1 db(m)|')
    rng = c*f/sweep_slope/2; %  beat2range(f,sweep_slope,c)
     
    % Visualize the spectrum
    %specanalyzer([txsig dechirpsig]);
    
    %xr(:,m) = dechirpsig;
    
%Est Range
[y,x] = max(mag2db(P1)) % find peak FFT point
rng(x) % map to range array

figure(5)
hold on
[Pyy,F] = periodogram(dechirpsig,[],2048,Fs,'centered');
%plot(F/1000,10*log10(Pyy));
xlabel('Frequency (kHz)');
ylabel('Power/Frequency (dB/Hz)');
title('Periodogram Power Spectral Density Estimate After Dechirping');
end
    figure(20)
    hold on
    sinewave = []
    for step_hop = 1:steps
    sinewave = vertcat(sinewave,wave(:,step_hop));
    end
    spectrogram(sinewave,32,16,32,fs,'yaxis');
    
    
%% END EDS CODE %%

  sig =  wave(:,1);
  Fs = fs;              % Sampling frequency
     T = 1/Fs;             % Sampling period
    L = length(sig);             % Length of signal
    t = (0:L-1)*T;        % Time vector
    window = hann(L);
    Y = ifft(sig);
    f = Fs*(0:(L/2))/L;
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    
    %rng = c*f/sweep_slope/2; %  RANGE PLOT IN M 
    figure(3)
    hold on
    %axis([0 6 -120 10])
    plot(f,mag2db(P1))
    title('Range Power Plot')
    xlabel('Range (m)')
    ylabel('|P1 db(m)|')


