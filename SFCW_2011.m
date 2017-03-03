function SFCW_2011
%% SFCW Radar Signal Simulation
% # The waveform generator combines the SFCW steps
% # The transmitter and the antenna amplify the signal and radiate the
%   signal into space.
% # The signal propagates to the target, gets reflected by the target, and
%   travels back to the radar.
% # The receiving antenna collects the signal.
% # The received signal is dechirped and saved in a buffer.
% # The signal is LPF filltered to remove the IQ DC jumps in the FFT.
% # FFT and range plotted 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% User Entry Here
fc = 24e9;       %24 Ghz is the system operating frequency
c  = 3e8;        %Speed of light 
Nsweep = 5;      %Number of sweep for the radar to perform with the radar (overlap the plots)

BW = 2e9;        %2Ghz System Bandwidth (higer bandwidth provide better resolution for target range
Fc = 4e6;        %Minimum frequency which dictates the number of steps ex 2Mhz so we hvae 1000 steps
                 %which also increases the range resolution 
                    
tot_sweep_time  = 1e-3;  % (s) long sweep times create large signal arrays (slow) 

Phase_NoiseAndOffset    = [-90,100e3]; %Noise and Offset 
SystemWhite_Noise       = -58;       %Iq Noise floor NOT USED IN THIS VERSION
Circulator_Isolation    = -30;       %Issolation in TX RX circulator coupling

distance_comm   = 2.5;    % (m) distance between the radar and commodity surface
comm_perm       = 2.3;    % (e) Commodity permitivity

%  End User Entry                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%% Start Sweep Code here:

lambda = c/fc;          %wavelength
FreqSteps = BW/Fc;       %calculate number of steps
fs  = BW*2;             %Samplign frequency at 4Ghz
tot_points = fs*tot_sweep_time; %points for given sweep time
points_per_step = tot_points/FreqSteps;
L = points_per_step;
wave = zeros(L,FreqSteps);
figure(1)



%create a sine wave for every step with the given number of points
for steps = 1:FreqSteps
   t = [0:1:points_per_step-1];
   I  = cos(((2*pi*(Fc*steps/(2*BW)*t))));
   Q  = sin(((2*pi*(Fc*steps/(2*BW)*t))));
   Z  = I + 1i*Q; %combine into a I Q type waveform
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



sig_combined = combineSteps(wave,FreqSteps); %Combine all steps into one wavefform

%% Sweep:

        for m = 1:Nsweep
            
        %% Add any phase noise
        sig = phase_noise(sig_combined,Phase_NoiseAndOffset(1),Phase_NoiseAndOffset(2));
        plotSweepSpectrum(sig,fs); %Plot the Spectrogram
        disp('Sweeping')
        m

        %% Setup the TX signal
        txsig = step(transmitter,sig);

        %% Propagate the signal and reflect off the target
        [radar_pos,radar_vel] = step(radarmotion,(1));
        txsig = step(channel,txsig,radar_pos,tgt_pos,radar_vel,radar_vel);
        txsig = step(target_comm,txsig); 
        
        
        %% Add Coupling
        txsig = circulator(Circulator_Isolation,sig,txsig);
        

        %% Received radar return
        txsig = step(receiver,txsig);
        txsig = phase_noise(txsig,Phase_NoiseAndOffset(1),Phase_NoiseAndOffset(2));

        %% Dechirp and LPF
        dechirpsig       = dechirp(txsig,sig);
        dechirpsig       = IQ_filter(dechirpsig); %Fillter the data through a LPF

        %% Plot FFT
        FFT_range(c,fs,dechirpsig,FreqSteps,BW,tot_sweep_time,Fc)
        end
end   


 
%% FilterDesigner LPF for filtering IQ jumps
%Filter created in MATLAB filterdesigner to filter square jumps 
%IQ_data: data to be filtered
%Returns: IQ_data passed through a LPF
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


 %% FFT Plotting and range (decimate,window,IFFT and plot the data)
 %speedOfLight: factor that was predefined
 %Fs: sample frequency
 %IQ_data: recived data
 %steps: number of steps the user entered
 %sweeptime: the sweeptime the user entered in seconds
 %stepSizeHz: step size in Hz of each step
 %Returns: combined signal at all steps 
 function FFT_range (speedOfLight,Fs,IQ_data,steps,BW,sweeptime,stepSizeHz)
    
    
    %TRIG CODE
    decimationFactor = length(IQ_data)/steps;
    IQ_data = decimate(IQ_data,decimationFactor); %Apply decimation
    
    Delta_F_Hz = stepSizeHz;
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

    Y = ifft((B*L*decimationFactor)); %undo the IFFT division by N and decimation division
    P2 = abs(Y/L);

    % TRIG VERSION
    figure(3)    
    plot(Xaxis,mag2db(P2)) 
    hold on
    axis([0 10 -100 10])
    title('Range Power Plot')
    xlabel('Range (m)')
    ylabel('|P1 db(m)|')
    title('SFCW IFFT Object Range and Magnitude');
    
    %Est Range
    [y,x] = max(mag2db(P2(1000:4000))); % find peak FFT point
    disp('Distance of object based on FFT (m):')
    Xaxis(x+1000)
 end

 %% Combining the steps in the waveform by combinbing each step
 %steps: number of steps in the waveform
 %wave: the IQ waveform data
 %Returns: combined signal at all steps 
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
 %fs: is sampling frequency
 %data: is the IQ data to be ploted
 %Returns: nothing
 function plotSweepSpectrum(data,fs)
 figure(2)
 data = complex(imag(data),real(data)); %Swap needed becuse spectrogram does FFT not IFFT
 spectrogram(data,32,16,32,fs,'yaxis');
 title('SFCW Signal Spectrogram/Sweep-time');
 end
  
 %% Adding IQ phasenoise
 %IQ_Data: is the original data to apply phase noise to
 %PhaseNoise: the ammount of phase noise to add in db
 %Offset: frequency offsets to apply phase noise to (Hz)
 %Returns: a phase noise mixed version of the IQ data
 function [IQ_Data_noise] = phase_noise(IQ_Data,PhaseNoise,Offset)
 pnoise = comm.PhaseNoise('Level',PhaseNoise,'FrequencyOffset',Offset, ...
     'SampleRate',2*Offset);
 IQ_Data_noise = step(pnoise,IQ_Data);
 %WhiteNoise    = awgn(IQ_Data_noise,1,.01);
 %IQ_Data_noise = WhiteNoise;
 end

 %% Adding Circulator Coupling (RX TX coupling)
 %isolation: dB of issolation between the circulator
 %initial: the TX signal
 %target: the RX signal
 %Returns: the initial recived signal with a portion of the TX signal 
 
 function [txsig_out] = circulator(isolation, initial, target)
    isolation = 10^(isolation/10); % convert from db to linear
    txsig_out = target + isolation * initial;
 end



