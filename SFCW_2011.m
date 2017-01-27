function SFCW_2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Entry Here


fc = 24e9;       %24 Ghz is the system operating frequency
c = 3e8;         %Speed of light 
Nsweep = 1;      %Number of sweep for the radar to perform with the radar (overlap the plots)

BW = 2e9;        %System Bandwidth (higer bandwidth provide better resolution for target range
Fc = 4e6;        %Minimum frequency which dictates the number of steps ex 2Mhz so we hvae 1000 steps
                    %which also increases the range resolution 
                    
tot_sweep_time  = 1e-3;  % (s) long sweep times create large signal arrays (slow) 

Phase_NoiseAndOffset    = [-80,100e3]; %Noise and Offset 
SystemWhite_Noise       = -60;       %Iq Noise floor NOT USED IN THIS VERSION
Circulator_Issolation   = -30;       %Issolation in TX RX circulator coupling

distance_comm   = 1.5;    % (m) distance between the radar and commodity surface
comm_perm       = 2.3;    % (e) Commodity permitivity

%  End User Entry                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%% Begin Test Code
lambda = c/fc;          %wavelength
FreqSteps = BW/Fc;       %calculate number of steps
fs  = BW*2;             %Samplign frequency at 4Ghz
tot_points = fs*tot_sweep_time; %points for given sweep time
points_per_step = tot_points/FreqSteps;
L = points_per_step;
wave = zeros(L,FreqSteps);
figure(1)



%create a sine wave for every step
for steps = 1:FreqSteps
   t = [0:1:points_per_step-1];
   I  = cos(((2*pi*(Fc*steps/(2*BW)*t))));
   Q  = sin(((2*pi*(Fc*steps/(2*BW)*t))));
   Z  = I - 1i*Q;
   wave(:,steps) = Z';
end
steps


%% Target Model
rcs_comm  = db2pow(min(10*log10(distance_comm)+5,20)); %RCS
c1 = 1/ sqrt((4*pi*10e-7)*(8.854*10e-12)*(comm_perm)); %Propagation speed calculation
target_comm = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',rcs_comm,'PropagationSpeed',c1,...
    'OperatingFrequency',fc);


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
         txsig = step(transmitter,sig);

        %% Propagate the signal and reflect off the target
        [radar_pos,radar_vel] = step(radarmotion,(1));

        txsig = step(channel,txsig,radar_pos,tgt_pos,radar_vel,radar_vel);

        txsig = step(target_comm,txsig); %step(H,X,UPDATE)

        %% Dechirp the received radar return
        txsig = step(receiver,txsig);


        %% Add Coupling, dechirp and LPF
        txsig = circulator(Circulator_Issolation,sig,txsig);
        dechirpsig       = dechirp(txsig,sig);
        dechirpsig       = IQ_filter(dechirpsig); %Fillter the data through a LPF

        %% Plot FFT
        FFT_range(c,fs,dechirpsig,FreqSteps,BW,tot_sweep_time)
        end
end   


 
%% FilterDesigner LPF for filtering IQ jumps
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
 function FFT_range (speedOfLight,Fs,IQ_data,steps,BW,sweeptime)
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
    
    dR = sweeptime*(speedOfLight/2);
    dR = dR/L;
    deltaR = linspace(dR,dR*L,L/2+1); %frequency to range conversion
    deltaR = deltaR - (deltaR(1));    %fixes distance offset - 0Hz = 0m
    
    figure(3)
    hold on
    axis([0 10 -100 10])
    plot(deltaR,mag2db(P1))
    title('Range Power Plot')
    xlabel('Range (m)')
    ylabel('|P1 db(m)|')
    title('FFT Object Range and Magnitude');
    %Est Range
    [y,x] = max(mag2db(P1)); % find peak FFT point
    disp('Distance of object based on FFT (m):')
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
  
 %% Adding IQ phasenoise
 function [IQ_Data_noise] = phase_noise(IQ_Data,PhaseNoise,Offset)
 pnoise = comm.PhaseNoise('Level',PhaseNoise,'FrequencyOffset',Offset, ...
     'SampleRate',2*Offset);
 IQ_Data_noise = step(pnoise,IQ_Data);
 WhiteNoise    = awgn(IQ_Data_noise,1,.01);
 IQ_Data_noise = WhiteNoise;
 end

 %% Adding Circulator Coupling
 function [txsig_out] = circulator(coupling_factor, initial, target)
    coupling_factor = 10^(coupling_factor/10);	
    txsig_out = target + coupling_factor * initial;
 end



