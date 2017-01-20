function FMCW

% The following table summarizes the radar parameters.
% 
%  System parameters            Value
%  ----------------------------------
%  Operating frequency (GHz)    24
%  Maximum target range (m)     4
%  Range resolution (m)         0.01
%  Maximum target speed (km/h)  0
%  Sweep time (seconds)         0.0001
         
fc = 24e9;  %24 Ghz Wave
c  = 3e8;   %Speed of light

range_max_meters   = 4;        %Bottom of tank in rail car
sweep_time_seconds = 0.0001;   %Use sweep of 0.1ms
range_res_meters   = 0.01;     %1 cm resolution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine Sweep bandwidth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = c/fc;  %wavelength
bw = range2bw(range_res_meters,c); %required bandwidth 
sweep_slope = bw/sweep_time_seconds; 

%maximum beat frequency can be calcualted
fr_max = range2beat(range_max_meters,sweep_slope,c);
fb_max = fr_max;


%%
% This example adopts a sample rate of the larger of twice the maximum beat
% frequency and the bandwidth.

fs = max(2*fb_max,bw);


%%
%Setup the radar FMCW wave
waveform = phased.FMCWWaveform('SweepTime',sweep_time_seconds,'SweepBandwidth',bw,...
    'SampleRate',fs);

sig = waveform();
figure(1)
spectrogram(sig,32,16,32,fs,'yaxis');
title('FMCW signal spectrogram');

%% Target Model
car_bottom = 2;
car_rcs = db2pow(min(10*log10(car_bottom)+5,20));

cartarget = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',car_rcs,'PropagationSpeed',c,...
    'OperatingFrequency',fc);

%cartarget = phased.RadarTarget('Model','Swerling2','MeanRCS',car_rcs,'PropagationSpeed',c,...
%   'OperatingFrequency',fc);

carmotion = phased.Platform();

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

%%
% Automotive radars are generally mounted on vehicles, so they are often in
% motion. This example assumes the radar is traveling at a speed of 100
% km/h along x-axis. So the target car is approaching the radar at a
% relative speed of 4 km/h.

radarmotion = phased.Platform();

%% Radar Signal Simulation
% # The waveform generator generates the FMCW signal.
% # The transmitter and the antenna amplify the signal and radiate the
% signal into space.
% # The signal propagates to the target, gets reflected by the target, and
% travels back to the radar.
% # The receiving antenna collects the signal.
% # The received signal is dechirped and saved in a buffer.
% # Once a certain number of sweeps fill the buffer, the Fourier transform
% is performed in both range and Doppler to extract the beat frequency as
% well as the Doppler shift. One can then estimate the range and speed of
% the target using these results. Range and Doppler can also be shown as an
% image and give an intuitive indication of where the target is in the
% range and speed domain.



        Nsweep = 4;
        % Alex comment
        % #justcapstonethings

        for m = 1:Nsweep
            % Update radar and target positions
            [radar_pos,radar_vel] = radarmotion(waveform.SweepTime);
            [tgt_pos,tgt_vel] = carmotion(waveform.SweepTime);

            % Transmit FMCW waveform
            sig = waveform();
            
            sig = phase_noise(sig);
            
            txsig_init = transmitter(sig);

            % Propagate the signal and reflect off the target
            tgt_pos(1) = 2;
            txsig = channel(txsig_init,radar_pos,tgt_pos,radar_vel,tgt_vel);
            UPDATE = false; %only used for Swirl1 stuff

            % txsig = cartarget(txsig,UPDATE); %step(H,X,UPDATE)
            txsig = cartarget(txsig); %step(H,X,UPDATE)

            % Dechirp the received radar return
            txsig = receiver(txsig);
            
            %   Mix original signal with received one:
            txsig = circulator(5, txsig_init, txsig);

            dechirpsig = dechirp(txsig,sig);

            %FFT and range of object
            FFT_range(c,fs,dechirpsig,sweep_slope)
            
            figure(2)
            spectrogram(sig,32,16,32,fs,'yaxis');
        end

end


function FFT_range (c,Fs,IQ_data,sweep_slope)
    %FFT
    T = 1/Fs;             % Sampling period
    L = length(IQ_data);  % Length of signal
    t = (0:L-1)*T;        % Time vector
    window = hann(L);
    Y = fft(IQ_data.*window);
    f = Fs*(0:(L/2))/L;
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    rng = c*f/sweep_slope/2; %  RANGE PLOT IN M 
    figure(3)
    hold on
    axis([0 6 -170 10])
    plot(rng,mag2db(P1))
    title('Range Power Plot')
    xlabel('Range (m)')
    ylabel('|P1 db(m)|')

    
    %Est Range
    [y,x] = max(mag2db(P1)); % find peak FFT point
    rng(x) % map to range array
end

function [txsig_out] = circulator(coupling_factor, initial, target)
    txsig_out = target + coupling_factor * initial;
end

function [IQ_Data] = phase_noise(IQ_Data)
    pnoise = comm.PhaseNoise('Level',-40,'FrequencyOffset',10);
    IQ_Data = pnoise(IQ_Data);
    IQ_Data = awgn(IQ_Data,40,'measured','db');
end

function dechirped_output = Mixer(Transmit_Waveform, Received_Waveform)
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

dechirped_output = dechirp(xref, x);

figure()
plot(real(dechirped_output))


[Pxx,F] = periodogram(dechirped_output,[],1024,Fs,'centered');
figure()
plot(F/1000,10*log10(Pxx)); grid;
xlabel('Frequency (kHz)');
ylabel('Power/Frequency (dB/Hz)');
title('Periodogram Power Spectral Density Estimate After Dechirping');
end
