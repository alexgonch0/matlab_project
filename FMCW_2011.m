function FMCW_2011
%% FMCW Radar Signal Simulation
% # The waveform generator generates the FMCW signal.
% # The transmitter and the antenna amplify the signal and radiate the
%   signal into space.
% # The signal propagates to the target, gets reflected by the target, and
%   travels back to the radar.
% # The receiving antenna collects the signal.
% # The received signal is dechirped and saved in a buffer.
% # FFT of beat freq and range plotted 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% User Entry Section
fc = 24e9;                    % 24 Ghz is the system operating frequency
c = 3e8;                      % Speed of light 
Nsweep = 1;

range_max_meters   = 4;       % Bottom of tank in rail car
tot_sweep_time     = 1e-3;    % Use sweep of 1ms (long sweeps create large arrays at high range resolution
range_res_meters   = 0.05;    % 5 cm resolution

PhaseNoiseLevel         = [-93, -120];  % Phase Noise levels in dBc/Hz
                                        %%% Try yourself [-91, -120] and [-95, -120]
PhaseNoiseOffset        = [2e9,  8e9];  % Phase Noise frequency offsets in Hz
SystemWhite_Noise       = -60;          % IQ Noise floor NOT USED IN THIS VERSION
Circulator_Isolation    = -20;          % Isolation in TX RX circulator coupling

distance_comm   = 1.5;        % (m) distance between the radar and commodity surface
comm_perm       = 2.3;        % (e) Commodity permitivity
%  End of User Entry




%% Start Sweep Code here:

lambda = c/fc;                       % wavelength
bw = range2bw(range_res_meters,c);   % bandwidth 
sweep_slope = bw/tot_sweep_time;     

% Maximum beat frequency
fr_max = range2beat(range_max_meters,sweep_slope,c);
fb_max = fr_max;

% We use sample rate of the larger of twice the maximum beat frequency and 
% the bandwidth
fs = max(2*fb_max,bw);


%% FMCW wave
waveform = phased.FMCWWaveform('SweepTime',tot_sweep_time,'SweepBandwidth',bw,...
    'SampleRate',fs);
sig = step(waveform);


%% Target Model
c1 = 1/ sqrt((4*pi*10e-7)*(8.854*10e-12)*(comm_perm));  % propagation speed calculation
rcs_comm = db2pow(min(10*log10(distance_comm)+5,20));   % cross-section of the commodity under the radar
target_comm = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',rcs_comm,...
    'PropagationSpeed',c1,'OperatingFrequency',fc);

%% Channel
% The propagation model is assumed to be free space.

channel = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

%% Radar System Setup
% The rest of the radar system includes the transmitter, the receiver, and
% the antenna. Currently the antenna is assumed to be isotropic and the
% gain of the antenna is included in the transmitter and the receiver.

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

%% Sweep:

        for m = 1:Nsweep
            %% Propagate the signal and reflect off the target
            [radar_pos,radar_vel] = step(radarmotion,(1));
            
            %% Add any phase noise
            sig = phase_noise(sig,PhaseNoiseLevel,PhaseNoiseOffset);
            plotSweepSpectrum(sig,fs); %Plot the Spectrogram
            disp('Sweeping')
            Nsweep

            
            %% Setup the TX signal
            txsig = step(transmitter,sig);


            %% Propagate the signal and reflect off the target
            txsig = step(channel,txsig,radar_pos,tgt_pos,radar_vel,radar_vel);
            txsig = step(target_comm,txsig); 


            %% Add Coupling, dechirp and LPF
            txsig = circulator(Circulator_Isolation,sig,txsig);
        

            %% Received radar return
            txsig = step(receiver,txsig);


            %% Dechirp 
            dechirpsig       = dechirp(txsig,sig);

            %% FFT Plot
            FFT_range(c,fs,dechirpsig,sweep_slope)


        end
end

%% FFT Plotting and range
function FFT_range (c,Fs,IQ_data,sweep_slope)
    % FFT
    T = 1/Fs;             % Sampling period
    L = length(IQ_data);  % Length of signal
    t = (0:L-1)*T;        % Time vector
    window = hann(L);
    Y = fft(IQ_data.*window);
    f = Fs*(0:(L/2))/L;
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    rng = c*f/sweep_slope/2; % RANGE PLOT IN M 
    figure(3)
    hold on
    axis([0 20 -100 10])
    plot(rng,mag2db(P1))
    title('Range Power Plot')
    xlabel('Range (m)')
    ylabel('|P1 db(m)|')
    title('FMCW FFT Object Range and Magnitude');
    
    % Range estimation
    [y,x] = max(mag2db(P1)); % find peak FFT point
    disp('Distance of object based on FFT (m):')
    rng(x) % map to range array
end

 %% Adding Circulator Coupling
 function [txsig_out] = circulator(isolation, initial, target)
    isolation = 10^(isolation/10);	
    txsig_out = target + isolation * initial;
 end


 %% Adding IQ phasenoise
 function [IQ_Data_noise] = phase_noise(IQ_Data,Level,Offset)
 pnoise = comm.PhaseNoise('Level',Level,'FrequencyOffset',Offset,...
          'SampleRate', 4*Offset(2));
 IQ_Data_noise = step(pnoise,IQ_Data);
 WhiteNoise    = awgn(IQ_Data_noise,1,.01);
 IQ_Data_noise = WhiteNoise;
 end
 
 %% Plotting Spectrogram
 function plotSweepSpectrum(data,fs)
 figure(1)
 spectrogram(data,32,16,32,fs,'yaxis');
 title('FMCW Signal Spectrogram/Sweep-time');
 end
 
 %% Not Yet Used
function dechirped_output = mixer(Transmit_Waveform, Received_Waveform)
% This is a stand-alone mixer function. It's not yet complete, I'm still
% working on it, but if you see anything that can be fixed or if you have
% any suggestions to improve it, just let me know.

% Transmit_Waveform and Received_Waveform must be FMCW Waveforms.

%Fs = 24e9; Tm = 0.0001
%hwav = phased.FMCWWaveform('SampleRate',Fs,'SweepTime',Tm);

%transmit_waveform = phased.FMCWWaveform('SampleRate',Fs,'SweepTime',Tm,'SweepBandwidth',100.0e3,'OutputFormat','Sweeps','NumSweeps',2);
figure()
plot(Transmit_Waveform)

xref = step(Transmit_Waveform);

%received_waveform = phased.FMCWWaveform('SampleRate',Fs,'SweepTime',Tm, 'SweepBandwidth',150.0e3,'OutputFormat','Sweeps','NumSweeps',2);
figure()
plot(Received_Waveform)

x = step(Received_Waveform);

dechirped_output = dechirp(x, xref);

phase_noise(dechirped_output);

figure()
plot(real(dechirped_output))

[Pxx,F] = periodogram(dechirped_output,[],1024,Fs,'centered');
figure()
plot(F/1000,10*log10(Pxx)); grid;
xlabel('Frequency (kHz)');
ylabel('Power/Frequency (dB/Hz)');
title('Periodogram Power Spectral Density Estimate After Dechirping');
end