
%% SFCW Radar Signal Simulation
% # The waveform generator combines the SFCW steps
% # The transmitter and the antenna amplify the signal and radiate the
%   signal into space.
% # The signal propagates to the target, gets reflected by the target, and
%   travels back to the radar.
% # The signal propagates through the medium, gets reflected by the bottom, and
%   travels back to the radar.
% # The signal delays are accounted for and shift the IQ signal accrodingly
% # The receiving antenna collects the signal.
% # The received signal is dechirped and saved in a buffer.
% # The signal is LPF filltered to remove the IQ DC jumps in the FFT.
% # FFT and range plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% User Entry Here
fc = 24e9;       % 24 GHz is the system operating frequency
c = 1/sqrt((4*pi*10^-7)*(8.854187*10^-12)); % propagation speed calculation = 3e8; % speed of light 
Nsweep = 3;      % Number of sweep for the radar to perform with the radar (overlap the plots)

BW = 2e9;        % 2 GHz System Bandwidth (higer bandwidth provides better resolution for target range)
Fc = 2e6;        % Frequency step size
                    
tot_sweep_time  = 3e-4;  % (s) long sweep times create large signal arrays (slow) 

Phase_NoiseAndOffset = [-80,100e3]; % Noise and Offset taken form data sheet
Circulator_Isolation = -20;         % Isolation in TX RX circulator coupling

slant_length    = 0.0115787; % (m) slant lenght of antenna
slant_angle     = 22;        % (degrees) slant angle
phys_ant_d      = 0.0381;

dist_comm       = 0.80;      % (m) distance between the radar and commodity surface
tank_h          = 3.20;      % (m) height of the tank
comm_perm       = 2.30;      % (e) Commodity permitivity
air_perm        = 1.00;
metal_perm      = 999;
sweepType       = 'up';      % quad_up, up, down
CALERROR        = true;      % non-linear calibration (deviations in calibration)
call_dev        = 3.0e4;     % (Hz) Calibration deviation form ideal (random)
% End User Entry                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Start Sweep Code
lambda = c/fc;          % wavelength
freqSteps = BW/Fc;      % calculate number of steps
fs = BW*2;              % sampling frequency at 2Fs = 4Ghz
tot_points = fs*tot_sweep_time; % points for given sweep time
points_per_step = tot_points/freqSteps;
L = points_per_step;
%wave = zeros(L,FreqSteps);

%% Acquire a sine way for the sweep
[wave, frequencyForCal] = generateSweepWaveform(Fc, BW, freqSteps, points_per_step, call_dev, CALERROR, sweepType, tot_sweep_time);

%% Target Model
rcs_comm = db2pow(min(10*log10(dist_comm)+5,20)); %RCS
c_comm = 1/ sqrt((4*pi*10^-7)*(8.854*10^-12)*(comm_perm)); %Propagation speed calculation for comm
target_comm = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',rcs_comm,'PropagationSpeed',c_comm,...
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
r = dist_comm*tan(slant_angle*pi/180); % for sloshing

ant_diameter = sqrt(3*lambda*slant_length);
effective_d = ant_diameter/phys_ant_d;
ant_gain = ((pi*ant_diameter)/lambda)^2 * effective_d;

tx_power = db2pow(ant_gain)*db2pow(1)*1e-3;     % in watts
tx_gain  = 15+ant_gain;                         % in dB

rx_gain  = 25+ant_gain;                         % RX LNA gain in dB
rx_nf    = 3;                                   % Noise Figure in dB

transmitter = phased.Transmitter('PeakPower',tx_power,'Gain',tx_gain);
receiver = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,'SampleRate',fs);
sig_combined = wave;% combineSteps(wave,FreqSteps); % combine all steps into one wavefform


%% Sweep:
for stepNumber = 1:Nsweep
    
    %% Add cal drift
    sig_combined = DriftCalibraton(wave,2.0e6,frequencyForCal,points_per_step,BW,CALERROR);

    %% Add phase noise
    sig = phase_noise(sig_combined,Phase_NoiseAndOffset(1),Phase_NoiseAndOffset(2));
    plotSweepSpectrum(sig,fs); % plot the Spectrogram
    disp('Sweeping')
    disp(num2str(stepNumber));
    

    
    %% Setup the TX signal
    txsig = step(transmitter,sig);

    %% Calcualate and apply pathloss in air to the transmitted signal from ant to commdity
    LfspOneWay  = pathLoss(0,dist_comm,fc,c); 
    txInterface = txsig * LfspOneWay; 

    %% If there is slowhing, get an RCS value
    rcs_comm = 4*pi^3*.45^4/lambda^2; % rcsSlosh(lambda,stepNumber,r,k)

    %% Calculate and apply the return signal from commdity back to ant and the delay
    returnsig  = txInterface*reflectionCoeff(air_perm,comm_perm); % return signal 
    LfspOneWay = pathLoss(0,dist_comm,fc,c); 
    returnsig  = returnsig*LfspOneWay;
    returnsig  = delaySignal(returnsig,dist_comm*2,fs,c);
    returnsig  = returnsig*rcs_comm; % rcs here



    %% Calculate and apply the signal from commdity to bottom with delays
    txInterfaceOil = txInterface*transmissionCoeff(air_perm,comm_perm); % transmition wave
    LfspTwoWay     = pathLoss(dist_comm,tank_h*2,fc,c_comm); % path loss in medium 2 way
    txInterfaceOil = txInterfaceOil*LfspTwoWay;
    txInterfaceOil = txInterfaceOil*reflectionCoeff(comm_perm,metal_perm); % reflect from bottom (mostly phase change of 180)
    txInterfaceOil = txInterfaceOil*transmissionCoeff(comm_perm,air_perm); % reflect oil to air boundry going back 
    txInterfaceOil = delaySignal(txInterfaceOil,(tank_h-dist_comm)*2,fs,c_comm); % delay in medium 2-way
    txInterfaceOil = delaySignal(txInterfaceOil,(dist_comm)*2,fs,c); % return delay air
    rcs = 4*pi^3*.45^4/lambda^2
    txInterfaceOil = txInterfaceOil*rcs; % rcs would go here


    %% Combined recived signal from commodity and bottom
    rxsig = txInterfaceOil + returnsig;

    %% Received radar return with gain
    rxsig = step(receiver,rxsig);

    %% Add coupling
    rxsig = circulator(Circulator_Isolation,txsig,rxsig);

    dechirpsig = dechirp(rxsig,txsig);
    %dechirpsig = IQ_filter(dechirpsig); % filter the data through a LPF
    %FIX FOR ANY SWEEP LENGTH

    %% Plot FFT
    FFT_range(c,fs,dechirpsig,freqSteps,BW,tot_sweep_time,Fc,tank_h,stepNumber,Nsweep)
end


 %% FilterDesigner LPF for filtering IQ jumps
 % Filter created in MATLAB filterDesigner to filter square jumps 
 % IQ_data: data to be filtered
 % Returns: IQ_data passed through a LPF
 function [filtered_data] = IQ_filter(IQ_data)
 % All frequency values are in Hz.
 Fs = 4e9;    % Sampling Frequency

 N  = 25;     % Order
 Fc = 350000; % Cutoff Frequency

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

 
 %% Creating a sine wave for every step with the given number of points
 % Fc (Hz): frequency step size
 % BW (Hz): system Bandwidth
 % freqSteps: number of frequency steps in one sweep
 % points_per_step: number of points in one step
 % call_dev (Hz): Calibration deviation form ideal (random) 
 % CALERROR (bool): calibration ON/OFF flag
 % sweepType: waveform of the sweep
 % Returns (V): wave
 function [wave, frequencyForCal] = generateSweepWaveform(Fc, BW, freqSteps, points_per_step, call_dev, CALERROR, sweepType, tot_sweep_time)
    %frequencyForCal = [];   
    wave = [];
  
    switch sweepType
        case 'up'
            for steps = 1:freqSteps
               t = 0:1:(points_per_step - 1);
                   if CALERROR % simulate small calibration errors if bool is true
                       randomCallError = -call_dev + (call_dev).*rand(1,1);
                   else
                       randomCallError = 0;
                   end
               I = cos(((2*pi*(((Fc+randomCallError)*steps)/(2*BW)*t))));
               Q = sin(((2*pi*(((Fc+randomCallError)*steps)/(2*BW)*t))));
               Z = I + 1i*Q; % combine into an IQ type waveform
               frequencyForCal(steps) = (Fc+randomCallError)*steps; % keep freqs
               wave = vertcat(wave,Z');
            end
        case 'down'
            for steps = 1:freqSteps
               t = 0:1:(points_per_step - 1);
                   if CALERROR % simulate small calibration errors if bool is true
                       randomCallError = -call_dev + (call_dev).*rand(1,1);
                   else
                       randomCallError = 0;
                   end
               I = cos(((2*pi*((((BW - ((steps-1)*Fc))*randomCallError))/(2*BW)*t))));
               Q = sin(((2*pi*((((BW - ((steps-1)*Fc))*randomCallError))/(2*BW)*t))));
               Z = I + 1i*Q; % combine into an IQ type waveform
               frequencyForCal(steps) = BW - ((steps-1)*Fc);
               wave = vertcat(wave,Z');
            end
        case 'quad_up'
           for steps = 1:freqSteps
               t = 0:1:(points_per_step - 1);
                   if CALERROR % simulate small calibration errors if bool is true
                       randomCallError = -call_dev + (call_dev).*rand(1,1);
                   else
                       randomCallError = 0;
                   end
               Fc = ((BW -0)/(tot_sweep_time^2))*(steps*tot_sweep_time/freqSteps)^2;
               I = cos(((2*pi*(((Fc+randomCallError))/(2*BW)*t))));
               Q = sin(((2*pi*(((Fc+randomCallError))/(2*BW)*t))));
               Z = I + 1i*Q; % combine into an IQ type waveform
               frequencyForCal(steps) = (Fc); % keep freqs
               wave = vertcat(wave,Z');
            end
        otherwise
            error('Non-existent sweep type has been defined. Terminating...');
    end     
 end

 %% FFT Plotting and range (decimate,window,IFFT and plot the data)
 % speedOfLight: factor that was predefined
 % Fs: sample frequency
 % IQ_data: recived data
 % steps: number of steps the user entered
 % sweeptime: the sweeptime the user entered in seconds
 % tank_h: height of the tank
 % stepSizeHz: step size in Hz of each step
 % stepNumber: current step
 % Nsweep: number of sweeps
 % Returns: combined signal at all steps 
 function FFT_range (speedOfLight,Fs,IQ_data,steps,BW,sweeptime,stepSizeHz,tank_h,stepNumber,Nsweep)
    %% Statistical data (Alex):
    % This stuff is populated with every sweep and processed later to obtain
    % standard deviations and average values.
    persistent peak_positions peak_magnitudes SNRs
    if stepNumber == 1
        peak_positions  = zeros(Nsweep, 1);
        peak_magnitudes = zeros(Nsweep, 1);
        SNRs            = zeros(Nsweep, 1);
    end
    
    %% FFT:
    decimationFactor = length(IQ_data)/steps;
    IQ_data = decimate(IQ_data,decimationFactor); % apply decimation
    
    Delta_F_Hz = stepSizeHz;
    B = IQ_data;
    L = length(B); % length of signal
    window = hann(L);
    B =  B.*window;
    
    L = round(speedOfLight/(Delta_F_Hz * 0.001) + 1); % res = 0.001 m
    
    B = [B; complex(zeros(L-length(IQ_data), 1))];
    
    Xaxis = 1:1:L;
    Xaxis = (((Xaxis*speedOfLight)/(Delta_F_Hz*(L - 1))));
    Xaxis = Xaxis./2;
    Xaxis = Xaxis - Xaxis(1);

    Y = ifft((B*L*decimationFactor)); % undo the IFFT division by N and decimation division
    P2 = abs(Y/L);
    
    figure(2)
    plot(Xaxis,mag2db(P2))
    hold on
    axis([0 10 -40 40])
    title('Range Power Plot')
    xlabel('Range (m)')
    ylabel('|P1 db(m)|')
    title('SFCW IFFT Object Range and Magnitude');
    
    %% Estimate Range:
    [y,x] = max(mag2db(P2(round(1/Xaxis(2)):round(5/Xaxis(2))))); % find peak FFT point between 1m to 5m
    peak_x = Xaxis(x+round(1/Xaxis(2)));
    disp('Distance of object based on FFT (m):')
    disp(num2str(peak_x))
    
    %% Calculate SNR:
    snr = calculateSNR(P2, Xaxis, peak_x, tank_h);
    snr_disp = ['SNR: ', num2str(round(snr*100)/100), ' dB']; % round to 2 digits after decimal point
    disp(snr_disp)
    
    %% Store statistical data:
    peak_positions(stepNumber)  = peak_x;
    peak_magnitudes(stepNumber) = y;
    SNRs(stepNumber)            = snr;
    
    %% Output statistical data:
    if stepNumber == Nsweep
        peak_positions_std_dev = std(peak_positions)*1000;
        peak_magnitudes_std_dev = std(peak_magnitudes);
        SNRs_mean = round(mean(SNRs)*100)/100;
        stats_disp = ['Peak positions std dev: ', num2str(peak_positions_std_dev), ' mm', ...
                     ' Peak magnitudes std dev: ', num2str(peak_magnitudes_std_dev), ' dB', ...
                     ' Average SNR: ', num2str(SNRs_mean), ' dB'];
        disp(stats_disp)
    end
 end

 %% Combining the steps in the waveform by combinbing each step
 % steps: number of steps in the waveform
 % wave: the IQ waveform data
 % Returns: combined signal at all steps 
 function [combined] = combineSteps(wave,steps)
 disp('Combining Waveforms... (this may take a bit)')
 wholesig = [] ;  
     for count = 1:steps
     sig =  wave(:,count);
     wholesig = vertcat(wholesig,sig);
     end
 combined = wholesig;
 end
 
 %% Plotting Spectrogram
 % fs: sampling frequency
 % data: the IQ data to be ploted
 % Returns: nothing
 function plotSweepSpectrum(data,fs)
 figure(1)
 data = complex(imag(data),real(data)); % swap needed because spectrogram does FFT not IFFT
 spectrogram(data,32,16,32,fs,'yaxis');
 title('SFCW Signal Spectrogram/Sweep-time');
 end
  
 %% Adding IQ phasenoise
 % IQ_Data: is the original data to apply phase noise to
 % PhaseNoise: the ammount of phase noise to add in db
 % Offset: frequency offsets to apply phase noise to (Hz)
 % Returns: a phase noise mixed version of the IQ data
 function [IQ_Data_noise] = phase_noise(IQ_Data,PhaseNoise,Offset)
 pnoise = comm.PhaseNoise('Level',PhaseNoise,'FrequencyOffset',Offset, 'SampleRate',2*Offset);
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
 
 %% Calculate Free Space Loss
 % startingDistance: dB of issolation between the circulator
 % endingDistance: the TX signal
 % frequency: the RX signal
 % c: speed of light in medium
 % Returns: the initial recived signal with a portion of the TX signal 
 function [pathloss_linear,pathloss_db] = pathLoss(startingDistance, endingDistance, frequency,c)
 pathloss_db  = 20*log10(endingDistance-startingDistance)+20*log10(frequency)+20*log10((4*pi)/c);
 pathloss_linear = (10^((-pathloss_db)/20));
 end
 
 %% Calculate reflection coeffeciant
 % n1: dielectric on side(1) entering
 % n2: dielectric on side(2) exiting
 % Returns: the reflection coeffeciant
 function [refCoeff] = reflectionCoeff(n1_enter, n2_exit)
 refCoeff = (sqrt(n1_enter) - sqrt(n2_exit))/(sqrt(n1_enter) + sqrt(n2_exit)); 
 end
 
 
 %% Calculate transmision coeffeciant
 % n1: dielectric on side(1) entering
 % n2: dielectric on side(2) exiting
 % Returns: the transmition coeffeciant
 function [trnCoeff] = transmissionCoeff(n1_enter, n2_exit)
 trnCoeff = (2*sqrt(n1_enter)/(sqrt(n1_enter) + sqrt(n2_exit))); 
 end
 
  
 %% Calculate signal delay and delay and array
 % distance: distance to apply delay
 % SignalIQ: Signal to apply fractional delay
 % fs: system sample rate
 % c: speed of light in medium
 % Returns: a delayed version of the signal
 function [delayedIQ] = delaySignal(SignalIQ,distance,fs,c)
 delay = distance/c;
 delayFs = delay*fs;
 fracDelay = dsp.VariableFractionalDelay;
 delayedIQ = fracDelay(SignalIQ,delayFs);
 end
 
 
 %% Sloshing
 % Dani code... need to comment
 function [rcs] = rcsSlosh(lambda,m,r,k)
        aspect_angle = 0;
        i = 0;
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
        rcs = pi*(k^2)*(r^4)*((2*besselj(1,2*k*r*sin(aspect_angle*pi/180))/...
            (2*k*r*sin(aspect_angle*pi/180)))^2)*cos(aspect_angle*pi/180)^2;
        else
        rcs = (4*pi^3*r^4)/ (lambda^2);
        end
 
 end
 
 %% SNR smart calculation (Alex)
 % y_data (array, ?) : ? from FFT
 % x_data (array, m) : values of x required to apply the window
 % peak_position (m) : actual distance to commodity acquired from FFT
 % tank_h (m)        : height of the tank
 % Returns: SNR in dB.
 function [snr_out] = calculateSNR(y_data, x_data, peak_position, tank_h)
    window = 0.5;
    
    %% Define window in terms of y_data indices:
    length_of_y = size(y_data, 1);    
    peak_index  = round(peak_position/x_data(end)*length_of_y);
    
    %% Collect noise floor data:
    noise_data = [];    
    
    % Window on the right from the peak:
    if peak_position < (0.5 + window) % first 0.5 m are assumed to be distorted by coupler isolation non-idealities
        marker_left = peak_index; % left boundary of the peak
        % Acquire data from the right of the peak:
        for i = peak_index : 1 : round(peak_index + 3*window/x_data(end)*length_of_y)
            if y_data(i) < y_data(i+1)
                marker_left = i;
                break
            end
        end
        window_index_right = round(marker_left + window/x_data(end)*length_of_y);
        noise_data = y_data(marker_left : window_index_right);
        disp(num2str(marker_left));
    % Window on the left from the peak (default):
    else
        marker_right = peak_index; % right boundary of the peak
        % Acquire data from the right of the peak:
        for i = peak_index : -1 : round(peak_index - 3*window/x_data(end)*length_of_y)
            if y_data(i) < y_data(i+1)
                marker_right = i;
                break
            end
        end
        window_index_left = round(marker_right - window/x_data(end)*length_of_y);
        noise_data = y_data(window_index_left : marker_right);
    end  
    
    %% Output SNR in dB:
    noise_average = mean(noise_data);
    snr_out = mag2db(y_data(peak_index)) - mag2db(noise_average);
 end
 
  %% Drift
 % drift the pre calibrated wave by a very small random frequency
 % divieation
 function [IQwave] = DriftCalibraton(wave,driftError,frequencyForCal,points_per_step,BW,CALERROR)
%% Create a sine wave for every step with the given number of points
    if CALERROR
    wave = [];
    L = length(frequencyForCal);
    driftError = 0.50;
    for steps = 1:L
       t = [0:1:points_per_step-1];
       randomCallError = (-(rand(1,1)+rand(1,1)))/250;

       I = cos(2*pi*(((frequencyForCal(steps)+randomCallError*frequencyForCal(steps))/(2*BW)*t)));
       Q = sin(2*pi*(((frequencyForCal(steps)+randomCallError*frequencyForCal(steps))/(2*BW)*t)));
       Z = I + 1i*Q; % combine into an IQ type waveform
       wave = vertcat(wave,Z');
       frequencyCal(steps) =  frequencyForCal(steps)+randomCallError*frequencyForCal(steps);
    end
    figure(11);
    hold on
    plot(frequencyCal);
    IQwave = wave;
    else
        IQwave = wave; %return old one
    end
    
 end
 
 
 
 

 