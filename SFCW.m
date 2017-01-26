
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Entry Here


fc = 24e9;      %24 Ghz Wave
c = 3e8;        %Speed of light
lambda = c/fc;  %wavelength
Nsweep = 1;     %Perform N sweeps of the radar and overlap the plots

BW = 2e9        %2Ghz BW
Fc = 4e6        %Minimum Freq ex 2Mhz so we hvae 1000 steps

Phase_NoiseAndOffset    = [-80,100e3] %Noise and Offset 
SystemWhite_Noise       = [-60]       %Iq Noise floor
Circulator_Issolation   = [-20];      %-20dB issolation

distance_comm   = 4;    % (m) distance between the radar and commodity surface
tot_sweep_time  = 1e-3  % (s) long sweep times create large signal arrays (slow)
comm_perm       = 2.3;  % (e) Commodity permitivity

%  End User Entry                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%% Begin Test Code
FreqSteps = BW/Fc 
fs  = BW*2;
tot_points = fs*tot_sweep_time
points_per_step = tot_points/FreqSteps
L = points_per_step;
wave = zeros(L,FreqSteps);
figure(1)

for steps = 1:FreqSteps
   t = [0:1:points_per_step-1];
   I  = cos(((2*pi*(Fc*steps/(2*BW)*t))));
   Q  = sin(((2*pi*(Fc*steps/(2*BW)*t))));
   Z  = I - 1i*Q;
   %wave = vertcat(wave,Z');
   %Z_Phase = phase_noise(Z',Phase_NoiseAndOffset(1),Phase_NoiseAndOffset(2));
   %wave(:,steps) = Z_Phase;
   wave(:,steps) = Z';
   %plot(I)
   %Frq = Fc*steps;
   %steps
   %plot(real(wave(:,steps)));
end
steps


%% Target Model
rcs_comm  = db2pow(min(10*log10(distance_comm)+5,20)); %RCS
c1 = 1/ sqrt((4*pi*10e-7)*(8.854*10e-12)*(comm_perm)); %Propagation speed calculation
target_comm = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',rcs_comm,'PropagationSpeed',c1,...
    'OperatingFrequency',fc);

%target_comm = phased.RadarTarget('Model','Swerling2','MeanRCS',rcs_comm,'PropagationSpeed',c1,...
%   'OperatingFrequency',fc);


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

tx_power = db2pow(ant_gain)*db2pow(1)*1e-3;     % in watts
tx_gain  = 9+ant_gain;                          % in dB

rx_gain = 15+ant_gain;                          % RX LNA gain in dB
rx_nf   = 3;                                    % Noise Figure in dB

tgt_pos = [distance_comm;0;0];                  % x y z position of target object.


transmitter = phased.Transmitter('PeakPower',tx_power,'Gain',tx_gain);
receiver = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,...
    'SampleRate',fs);
radarmotion = phased.Platform('InitialPosition',[0;0;0]);


%% Actual Sweep

sig_combined = combineSteps(wave,FreqSteps); %Combine all steps into one wavefform


for m = 1:Nsweep
    %% Add any phase noise
    sig = phase_noise(sig_combined,Phase_NoiseAndOffset(1),Phase_NoiseAndOffset(2));
    plotSweepSpectrum(sig,fs); %Plot the Spectrogram
    disp('Sweeping')
    Nsweep
    
    %% Setup the TX signal
    txsig =  transmitter(sig);
    
    %% Propagate the signal and reflect off the target
    [radar_pos,radar_vel] = radarmotion(1);
    txsig = channel(txsig,radar_pos,tgt_pos,radar_vel,radar_vel);
    txsig = target_comm(txsig); %step(H,X,UPDATE)
    
    %% Dechirp the received radar return
    txsig = receiver(txsig);    
    
    %% Add Coupling dechirp and LPF
    txsig = circulator(Circulator_Issolation,sig,txsig);
    dechirpsig       = dechirp(txsig,sig);
    dechirpsig       = IQ_filter(dechirpsig); %Fillter the data through a LPF

    %% Plot FFT
    FFT_range(c,fs,dechirpsig,FreqSteps,BW)
   

end
    
 clear all;

 
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
%figure(4)
%hold on
%plot(real(IQ_data))
%plot(real(output))
%plot(imag(output))
end


 %% FFT Plotting and range
 function FFT_range (speedOfLight,Fs,IQ_data,steps,BW)
    %FFT
    %IQ_data = decimate(IQ_data,1000);
    %Fs = Fs/1000;
    T = 1/Fs;             % Sampling period
    L = length(IQ_data);  % Length of signal
    t = (0:L-1)*T;        % Time vector
    window = hann(L);
    Y = fft(IQ_data.*window);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    L = length(IQ_data);  % Length of signal
    f = Fs*(0:(L/2))/L;
   
    dR = speedOfLight/(2*steps*(BW/steps));
    dR = dR/L;
    deltaR = linspace(dR,dR*steps*L,L/2+1); 
    deltaR = (deltaR.*8000);
    
    figure(3)
    hold on
    axis([0 10 -100 -10])
    plot(deltaR,mag2db(P1))
    title('Range Power Plot')
    xlabel('Range (m)')
    ylabel('|P1 db(m)|')
    title('FFT Object Range and Magnitude');
    %Est Range
    [y,x] = max(mag2db(P1)); % find peak FFT point
    deltaR(x)
 end

 %% Combining the steps in the waveform
 function [combined] = combineSteps(wave,steps)
 disp('Combining Waveforms(this may take a bit)')
 wholesig = [] ;  
     for count = 1:steps
     sig =  wave(:,count);
     wholesig = vertcat(wholesig,sig);
     end
combined = wholesig;
 end
 
 %% Plotting Spectrogram
 function plotSweepSpectrum(data,fs)
 figure(1)
 spectrogram(data,32,16,32,fs,'yaxis');
 title('SFCW Signal Spectrogram/Sweep-time');
 end
  
 %% Estimate Range
 function [res] = range_resolution(c,steps,BW,fs)
 L = length(dechirpsig);  % Length of signal
 f = fs*(0:(L/2))/L;
 deltaR = c ./ (2*f);
 end

 %% Adding IQ phasenoise
 function [IQ_Data_noise] = phase_noise(IQ_Data,PhaseNoise,Offset)
 pnoise = comm.PhaseNoise('Level',PhaseNoise,'FrequencyOffset',Offset, ...
     'SampleRate',2*Offset);
 IQ_Data_noise = pnoise(IQ_Data);
 WhiteNoise    = awgn(IQ_Data_noise,1,.01);
 IQ_Data_noise = WhiteNoise;
 end

 %% Adding Circulator Coupling
 function [txsig_out] = circulator(coupling_factor, initial, target)
    coupling_factor = 10^(coupling_factor/10);	
    txsig_out = target + coupling_factor * initial;
 end



