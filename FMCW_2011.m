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


%% User Entry Here
fc = 24e9;       %24 Ghz is the system operating frequency
c  = 3e8;         %Speed of light 
Nsweep = 5;

range_max_meters   = 4;       % Bottom of tank in rail car
tot_sweep_time     = 1e-3;    % Use sweep of 1ms (long sweeps create large arrays at high range resolution
range_res_meters   = 0.001;   % 1 cm resolution

Phase_NoiseAndOffset    = [-65,100e3]; %Noise and Offset 
SystemWhite_Noise       = -60;       %Iq Noise floor NOT USED IN THIS VERSION
Circulator_Issolation   = -30;       %Issolation in TX RX circulator coupling

slant_length    = 0.0115787; % (m) slant lenght of antenna
beam_shape      = 'Circular'; %Circular or Elliptic choice
aspect_angle    = 0; %Define angle of incident beam
phys_ant_d       = 0.0381;

distance_comm   = 2.5;    % (m) distance between the radar and commodity surface
comm_perm       = 2.3;    % (e) Commodity permitivity

%Decimation will allow us to cut the data into 'decimation" number of
%points to replicate SFCW behavior and use SFCW type of FFT
BOOL_DECIMATE = true;     % true enables decimation, false disables it
PointsAfterdecimation    = 2000;      % make it such that there is a 'decimation' number of points as if SFCW


%  End User Entry                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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
sig_combined = step(waveform);


%% Target Model
c1 = 1/ sqrt((4*pi*10e-7)*(8.854*10e-12)*(comm_perm));  %Propagation speed calculation
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

ant_diameter = sqrt(3*lambda*slant_length);
effective_d = ant_diameter/phys_ant_d;
ant_gain = ((pi*ant_diameter)/lambda)^2 * effective_d;

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
            sig = phase_noise(sig_combined,Phase_NoiseAndOffset(1),Phase_NoiseAndOffset(2));
            plotSweepSpectrum(sig,fs); %Plot the Spectrogram
            disp('Sweeping')
            Nsweep;

            
            %% Setup the TX signal
            txsig = step(transmitter,sig);


            %% Propagate the signal and reflect off the target
            txsig = step(channel,txsig,radar_pos,tgt_pos,radar_vel,radar_vel);
            txsig = step(target_comm,txsig); 


            %% Add Coupling, dechirp and LPF
            txsig = circulator(Circulator_Issolation,sig,txsig);
        

            %% Received radar return
            txsig = step(receiver,txsig);


            %% Dechirp 
            dechirpsig       = dechirp_mixer(txsig, sig, m, Nsweep);

            %% FFT Plot
            FFT_range(c,fs,dechirpsig,sweep_slope,bw,BOOL_DECIMATE,PointsAfterdecimation)


        end
end

%% FFT Plotting and range
function FFT_range (c,Fs,IQ_data,sweep_slope,bw,BOOL_DECIMATE,decimation)

    if BOOL_DECIMATE == true
    speedOfLight = c;
    Delta_F_Hz = bw/decimation;
    
    decimationFactor = length(IQ_data)/decimation;
    IQ_data = decimate(IQ_data,decimationFactor); %
    B = IQ_data;
    L = length(B);        % Length of signal
    window = hann(L);
    B =  B.*window;  
    L = round(speedOfLight/(Delta_F_Hz * 0.001 )+1); %res = 0.001 m
    B = [B; complex(zeros(L-length(IQ_data), 1))];
    
    Xaxis = 1:1:L;
    Xaxis = (((Xaxis*speedOfLight)/(Delta_F_Hz*(L - 1))));
    Xaxis = Xaxis./2;
    Xaxis = Xaxis - Xaxis(1);
    
    Y = fft(B); 
    
    P2 = abs(Y/L);
    figure(3)
    hold on
    axis([0 5 -120 -20])
    plot(Xaxis,mag2db(P2)) 
    title('Range Power Plot')
    xlabel('Range (m)')
    ylabel('|P1 db(m)|')
    title('FMCW FFT Decimated Object Range and Magnitude');
    
    else
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
    
    rng = c*f/sweep_slope/2; % rng is RANGE PLOT IN M 
    figure(3)
    hold on
    axis([0 5 -100 10])
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
end

 %% Adding Circulator Coupling
 function [txsig_out] = circulator(coupling_factor, initial, target)
    coupling_factor = 10^(coupling_factor/10);	
    txsig_out = target + coupling_factor * initial;
 end


 %% Adding IQ phasenoise
 function [IQ_Data_noise] = phase_noise(IQ_Data,PhaseNoise,Offset)
 pnoise = comm.PhaseNoise('Level',PhaseNoise,'FrequencyOffset',Offset, ...
     'SampleRate',2*Offset);
 IQ_Data_noise = step(pnoise,IQ_Data);
 %WhiteNoise    = awgn(IQ_Data_noise,1,.01);
 %IQ_Data_noise = WhiteNoise;
 end
 
 %% Plotting Spectrogram
 function plotSweepSpectrum(data,fs)
 figure(1)
 spectrogram(data,32,16,32,fs,'yaxis');
 title('FMCW Signal Spectrogram/Sweep-time');
 end
 
 %% Not Yet Used
function [dechirp_out] = dechirp_mixer(received, transmitted, m, Nsweep)
% Input variables "received" and "transmitted" must be Nx1 complex doubles of the
% same length.
% Input variables 'm' and 'Nsweep' are used to create an output periodogram after
% the final sweep.
% Input variables 'noise_lv' and 'noise_off' are used to represent the
% phase noise level and phase noise offset, respectively. [-180, -190], [2e9,  8e9]

Fs = 24e9; %Sample Frequency set to 24 GHz.

clean_out = dechirp(received, transmitted); % Mixes received and transmitted signals


pnoise = comm.PhaseNoise('Level', [-180, -190], 'FrequencyOffset', [2e9,  8e9],'SampleRate', 4*(8e9));

conv_loss = 0.93; % Conversion Loss from mixer

dechirp_out = conv_loss*(step(pnoise, clean_out)); % Adds phase noise, saves to output


% Prints output periodogram on final sweep
if m == Nsweep

    [Pyy,F] = periodogram(dechirp_out,[],1024,Fs,'centered');
    figure();
    plot(F/1000,1*log10(Pyy));
    grid;
    xlabel('Frequency (kHz)');
    ylabel('Power/Frequency (dB/Hz)');
    title('Power Spectral Density Estimate After Mixer Stage'); 
    
end


end
