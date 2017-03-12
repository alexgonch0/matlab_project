
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
fc = 24e9;       % 24 Ghz is the system operating frequency
c  = 3e8;        % Speed of light 
Nsweep = 1;      % Number of sweep for the radar to perform with the radar (overlap the plots)

BW = 2e9;        % 2GHz System Bandwidth (higer bandwidth provides better resolution for target range)
Fc = 1e6;        % Minimum frequency which dictates the number of steps (e.g. 2MHz) so we have
                 % 1000 steps which also increases the range resolution 
                    
tot_sweep_time = 1e-3;              % (s) long sweep times create large signal arrays (slow) 

Phase_NoiseAndOffset = [-80,100e3]; % Noise and Offset 
SystemWhite_Noise = -58;            % Iq Noise floor [NOT USED IN THIS VERSION]
Circulator_Isolation = -20;         % Issolation in TX RX circulator coupling

slant_length = 0.0115787;           % (m) slant lenght of antenna
slant_angle = 22;                   % (degress) 
phys_ant_d = 0.0381;

dist_comm = 1;                      % (m) distance between the radar and commodity surface
tank_h = 4.71;
comm_perm = 2.3;                    % (e) Commodity permitivity
air_perm = 1;
metal_perm = 999;
CALERROR = true;
%  End User Entry                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%% Start Sweep Code here:

lambda = c/fc;                  % wavelength
FreqSteps = BW/Fc;              % calculate number of steps
fs  = BW*2;                     % sampling frequency at 4 GHz
tot_points = fs*tot_sweep_time; % points for given sweep time
points_per_step = tot_points/FreqSteps;
L = points_per_step;
wave = zeros(L,FreqSteps);
% figure(1) <---- what is this? (Alex)



% Create a sine wave for every step with the given number of points:
for steps = 1:FreqSteps
   t = [0:1:points_per_step-1];
       if CALERROR % simulate small calibration errors if bool is true
           randomCallError = -1e4 + (1e4).*rand(1,1);
       else
           randomCallError = 0;
       end
   I  = cos(((2*pi*(((Fc+randomCallError)*steps)/(2*BW)*t))));
   Q  = sin(((2*pi*(((Fc+randomCallError)*steps)/(2*BW)*t))));
   Z  = I + 1i*Q; % combine into an IQ type waveform
   wave(:,steps) = Z';
end
steps


%% Target Model
rcs_comm  = db2pow(min(10*log10(dist_comm)+5,20));     % RCS
c1 = 1/ sqrt((4*pi*10e-7)*(8.854*10e-12)*(comm_perm)); % Propagation speed calculation
target_comm = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',rcs_comm,'PropagationSpeed',c1,...
    'OperatingFrequency',fc);


%% Channel
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

k = 2*pi/lambda;
r = dist_comm*tan(slant_angle*pi/180);

ant_diameter = sqrt(3*lambda*slant_length);
effective_d = ant_diameter/phys_ant_d;
ant_gain = ((pi*ant_diameter)/lambda)^2 * effective_d;

tx_power = db2pow(ant_gain)*db2pow(1)*1e-3;     % in watts
tx_gain  = 15+ant_gain;                         % in dB

rx_gain = 22+ant_gain;                          % RX LNA gain in dB
rx_nf   = 3;                                    % Noise Figure in dB

tgt_pos = [dist_comm;0;0];                  % x y z position of target object.


transmitter = phased.Transmitter('PeakPower',tx_power,'Gain',tx_gain);
receiver = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,...
    'SampleRate',fs);
radarmotion = phased.Platform('InitialPosition',[0;0;0]);


sig_combined = combineSteps(wave,FreqSteps); % ?ombine all steps into one wavefform
aspect_angle = 0;
i = 0;

%% Sweep:

        for m = 1:Nsweep
        
        motion = 0.5;
        stop_condition = 5;
        if m ~= 1 && i == 0
        aspect_angle = aspect_angle + motion;
            if aspect_angle == stop_condition
                i = stop_condition/motion;
            end
        elseif m ~= 1 
        aspect_angle = aspect_angle - motion;
        i = i-1;
        end
            
        if aspect_angle ~= 0
        rcs_comm = pi*(k^2)*(r^4)*((2*besselj(1,2*k*r*sin(aspect_angle*pi/180))/...
            (2*k*r*sin(aspect_angle*pi/180)))^2)*cos(aspect_angle*pi/180)^2;
        else
        rcs_comm = (4*pi^3*r^4)/ (lambda^2);
        end
            
        %% Add any phase noise
        sig = phase_noise(sig_combined,Phase_NoiseAndOffset(1),Phase_NoiseAndOffset(2));
        plotSweepSpectrum(sig,fs); % plot the Spectrogram
        disp('Sweeping')

        %% Setup the TX signal
        txsig = step(transmitter,sig);
        
        
        % the delay ? is given by R/c, where R is the propagation distance
        % and c is the propagation speed. The free-space path loss is given by
        
        %Lfsp = 20*log10(distance_comm*2)+20*log10(fc)+20*log10((4*pi)/c); % tx/rx gain already in signal
        Lfsp_oneway = 20*log10(dist_comm)+20*log10(fc)+20*log10((4*pi)/c)
        
        powInterface = (10^((-Lfsp_oneway)/20)); %path loss
        txInterface = txsig * powInterface; %power at interface using through path loss
        refCoeff  = (air_perm - comm_perm)/(air_perm + comm_perm) %reflect signal
        returnsig =  txInterface*abs(refCoeff); %return signal 
        
        trmCoeff  = (1 - abs(refCoeff)) %transmition coeff
        txInterfaceOil = trmCoeff * txInterface; %transmition wave
        refCoeff  = (comm_perm - metal_perm)/(comm_perm + metal_perm) %reflect with bottom
        txInterfaceOil = txInterfaceOil*abs(refCoeff);
        delayOil  =  ((tank_h-dist_comm)*2)/(c/1.46);
        delayOil  =  delayOil + ((dist_comm*2)/c);
        delayOilFs = delayOil*fs;
        NOil = round(delayOilFs); %ammount to delay
        rxsigOil = [zeros(NOil,1);txInterfaceOil];
        rxsigOil = rxsigOil(1:end-NOil);
        
        
        lfactor = 10^((-Lfsp_oneway)/20); %next trip back free space path loss
        returnsig = returnsig*lfactor;
        delay =  (dist_comm*2)/c;
        delayFs = delay*fs;
        
        N = round(delayFs); %ammount to delay
        rxsig = [zeros(N,1);returnsig];
        rxsig = rxsig(1:end-N);
        
        %rxsig = rxsig.*lfactor; %apply loss
        %apply backscatter
        %Y = sqrt(G)*X
        rxsig = sqrt(rcs_comm)*rxsig;
        
        rxsig = rxsig + rxsigOil;
        
        %% Received radar return with gain
        rxsig = step(receiver,rxsig);
        
        %% Add Coupling
        rxsig = circulator(Circulator_Isolation,txsig,rxsig);

        dechirpsig = dechirp(rxsig,txsig);
        dechirpsig = IQ_filter(dechirpsig); %Fillter the data through a LPF

        %% Plot FFT
        FFT_range(c,fs,dechirpsig,FreqSteps,BW,tot_sweep_time,Fc)
        end
   


 
%% FilterDesigner LPF for filtering IQ jumps
%Filter created in MATLAB filterdesigner to filter square jumps 
%IQ_data: data to be filtered
%Returns: IQ_data passed through a LPF
function [filtered_data] = IQ_filter(IQ_data)
% All frequency values are in Hz.
Fs = 4e9;     % Sampling Frequency

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
 % speedOfLight: factor that was predefined
 % Fs: sample frequency
 % IQ_data: recived data
 % steps: number of steps the user entered
 % sweeptime: the sweeptime the user entered in seconds
 % stepSizeHz: step size in Hz of each step
 % Returns: combined signal at all steps 
 function FFT_range(speedOfLight,Fs,IQ_data,steps,BW,sweeptime,stepSizeHz)
    % TRIG CODE
    decimationFactor = length(IQ_data)/steps;
    IQ_data = decimate(IQ_data,decimationFactor); % Apply decimation
    
    Delta_F_Hz = stepSizeHz;
    B = IQ_data;
    L = length(B); % Length of signal
    window = hann(L);
    B = B.*window;  
    
    L = round(speedOfLight/(Delta_F_Hz * 0.001) + 1); % res = 0.001 m

    B = [B; complex(zeros(L-length(IQ_data), 1))];
    
    Xaxis = 1:1:L;
    Xaxis = (((Xaxis*speedOfLight)/(Delta_F_Hz*(L - 1))));
    Xaxis = Xaxis./2;
    Xaxis = Xaxis - Xaxis(1);

    Y = ifft((B*L*decimationFactor)); % undo the IFFT division by N and decimation division
    P2 = abs(Y/L);
    
    % TRIG VERSION
    figure(2)    
    plot(Xaxis,mag2db(P2)) 
    hold on
    axis([0 10 -40 40])
    title('Range Power Plot')
    xlabel('Range (m)')
    ylabel('|P1 db(m)|')
    title('SFCW IFFT Object Range and Magnitude');
    
    % Estimate Range
    [y,x] = max(mag2db(P2(round(1/Xaxis(2)):round(5/Xaxis(2))))); % find peak FFT point 1m to 5m
    disp('Distance of object based on FFT (m):')
    Xaxis(x+round(1/Xaxis(2)))
    
    % Calculate SNR (you can read the details in the function scope below)
    peak = Xaxis(x+round(1/Xaxis(2)));
    window = 1.0;
    snr = calculateSNR(P2, Xaxis, peak, window);
    snr_disp = ['SNR: ', num2str(round(snr*100)/100), ' dB']; % round to 2 digits after decimal point
    disp(snr_disp)
 end

 %% Combining the steps in the waveform by combinbing each step
 % steps: number of steps in the waveform
 % wave: the IQ waveform data
 % Returns: combined signal at all steps 
 function [combined] = combineSteps(wave,steps)
 disp('Combining Waveforms... (this may take a bit)')
 wholesig = [] ;  
     for count = 1:steps
     sig = wave(:,count);
     wholesig = vertcat(wholesig,sig);
     end
 combined = wholesig;
 end
 
 %% Plotting Spectrogram
 % fs: is sampling frequency
 % data: is the IQ data to be ploted
 % Returns: nothing
 function plotSweepSpectrum(data,fs)
 figure(1)
 data = complex(imag(data),real(data)); %Swap needed becuse spectrogram does FFT not IFFT
 spectrogram(data,32,16,32,fs,'yaxis');
 title('SFCW Signal Spectrogram/Sweep-time');
 end
  
 %% Adding IQ phasenoise
 % IQ_Data: is the original data to apply phase noise to
 % PhaseNoise: the ammount of phase noise to add in db
 % Offset: frequency offsets to apply phase noise to (Hz)
 % Returns: a phase noise mixed version of the IQ data
 function [IQ_Data_noise] = phase_noise(IQ_Data,PhaseNoise,Offset)
 pnoise = comm.PhaseNoise('Level',PhaseNoise,'FrequencyOffset',Offset, ...
     'SampleRate',2*Offset);
 IQ_Data_noise = step(pnoise,IQ_Data);
 %WhiteNoise    = awgn(IQ_Data_noise,1,.01);
 %IQ_Data_noise = WhiteNoise;
 end

 %% Adding Circulator Coupling (RX TX coupling)
 % isolation: dB of issolation between the circulator
 % initial: the TX signal
 % target: the RX signal
 % Returns: the initial recived signal with a portion of the TX signal 
 
 function [txsig_out] = circulator(isolation, initial, target)
    isolation = 10^(isolation/10); % convert from db to linear
    txsig_out = target + isolation * initial;
 end

 %% SNR smart calculation
 % y_data (array, ?) : power from FFT
 % x_data (array, m) : values of x required to apply the window
 % peak_position (m) : actual distance to commodity acquired from FFT
 % window (m)        : the range away from the peak to be used in avg noise power calculation
 % Returns: SNR in dB.
 function [snr_out] = calculateSNR(y_data, x_data, peak_position, window)
    figure(3)
    plot(x_data, y_data);
    axis([0 (peak_position + window) -1 9])
    
    %% Define window in terms of y_data indices:
    length_of_y        = size(y_data, 1);    
    peak_index         = round(peak_position/x_data(end)*length_of_y);
    window_index_left  = round(peak_index - (window/2)/x_data(end)*length_of_y);
    window_index_right = round(peak_index + (window/2)/x_data(end)*length_of_y);
    
    debug_display = ['index_left: ', num2str(window_index_left), ' peak_index: ', num2str(peak_index), ' index_right: ', num2str(window_index_right)];
    disp(debug_display)
    
    %% Collect noise floor data:
    marker_left  = peak_index;   % left boundary of the peak
    marker_right = peak_index;   % right boundary of the peak
    
    % Acquire data from the left of the peak:
    for i = (peak_index - 1) : -1 : window_index_left
        if y_data(i) < y_data(i-1)
            marker_left = i;
            break
        end
    end
    
    % Acquire data from the right of the peak:
    for i = peak_index : 1 : window_index_right
        if y_data(i) < y_data(i+1)
            marker_right = i;
            break
        end
    end
    
    % Calculate average noise power:
    noise_data_left  = y_data(window_index_left : marker_left);
    noise_data_right = y_data(marker_right : window_index_right);
    noise_average    = (mean(noise_data_left) + mean(noise_data_right))/2;
    
    % Output SNR in dB:
    snr_out = mag2db(y_data(peak_index)) - mag2db(noise_average);
 end
