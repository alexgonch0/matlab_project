
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
Nsweep = 10;     % Number of sweep for the radar to perform with the radar (overlap the plots)

BW = 2e9;        % 2 GHz System Bandwidth (higer bandwidth provides better resolution for target range)
freqStepSize = 1e6;        % Frequency step size
                    
tot_sweep_time  = 1e-4;  % (s) long sweep times create large signal arrays (slow) 

Phase_NoiseAndOffset = [-150,100e3]; % Noise and Offset taken form data sheet
Circulator_Isolation = -20;         % Issolation in TX RX circulator coupling

slant_length    = 0.0115787; % (m) slant lenght of antenna
slant_angle     = 22;        % (degrees) slant angle
phys_ant_d      = 0.0381;

dist_comm       = 2.00;      % (m) distance between the radar and commodity surface
tank_h          = 3.20;
comm_perm       = 2.30;      % (e) Commodity permitivity
air_perm        = 1.00;
metal_perm      = 999;
sweepType       = 'up'; %quad_up , up ,down
CALERROR        = true;      % non-linear calibration (deviations in calibration)
call_dev        = 1.0e4;     % (Hz) Calibration deviation form ideal (random)
% End User Entry                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Start Sweep Code
lambda = c/fc;          % wavelength
freqSteps = BW/freqStepSize;      % calculate number of steps
fs =BW*2;               % sampling frequency at 2Fs = 4Ghz
tot_points = fs*tot_sweep_time; % points for given sweep time
points_per_step = tot_points/freqSteps;
L = points_per_step;
%wave = zeros(L,FreqSteps);



%% Target Model
rcs_comm = db2pow(min(10*log10(dist_comm)+5,20)); %RCS
c_comm = 1/ sqrt((4*pi*10^-7)*(8.854187*10^-12)*(comm_perm)); %Propagation speed calculation for comm
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

%% Sweep:
%http://micsymposium.org/mics_2009_proceedings/mics2009_submission_64.pdf
%% Acquire a sine way for the sweep
[frequencyForCal] = generateSweepWaveform(BW,fc,freqSteps,call_dev, CALERROR, sweepType,tot_sweep_time);

error = 0.0001;
driftedCalFreq = [];
Idata = [];
Qdata = [];

for sweepNumber = 1:Nsweep
driftedCalFreq = DriftCalibraton(error,frequencyForCal,CALERROR);
    for stepNumber = 1:freqSteps
    %frequencyForCal
    propFreq = driftedCalFreq(stepNumber)+fc;
    txsig = step(transmitter,1); %1 volt is base power
    

    %% Calcualate and apply pathloss in air to the transmitted signal from ant to commdity
    LfspOneWay  = pathLoss(error,dist_comm,propFreq,c); 
    txInterface = txsig * LfspOneWay; 

    %% If there is slowhing, get an RCS value
    rcs_comm = 4*pi^3*(.1*dist_comm)^4/lambda^2;

    %% Calculate and apply the return signal from commdity back to ant and the delay
    returnsig  = txInterface*reflectionCoeff(air_perm,comm_perm); % return signal 
    LfspOneWay = pathLoss(0,dist_comm,propFreq,c); 
    returnsig  = returnsig*LfspOneWay;
    returnsig  = returnsig*rcs_comm; % rcs here

    
    %% Calculate and apply the signal from commdity to bottom with delays
    txInterfaceOil = txInterface*transmissionCoeff(air_perm,comm_perm); % transmition wave
    LfspOneWay     = pathLoss(dist_comm,tank_h,fc,c_comm); % path loss in medium 1 way
    txInterfaceOil = txInterfaceOil*LfspOneWay;
    txInterfaceOil = txInterfaceOil*reflectionCoeff(comm_perm,metal_perm); % reflect from bottom (mostly phase change of 180)
    txInterfaceOil = txInterfaceOil*LfspOneWay;
    txInterfaceOil = txInterfaceOil*transmissionCoeff(comm_perm,air_perm); % reflect oil to air boundry going back 
    rcs = 4*pi^3*(.2*tank_h)^4/lambda^2;
    txInterfaceOil = txInterfaceOil*rcs; % rcs would go here
    
    %% Combined recived signal from commodity and bottom
    Idata(stepNumber) = returnsig*(cos(-2*pi*propFreq*dist_comm*2/c));
    Qdata(stepNumber) = returnsig*(sin(-2*pi*propFreq*dist_comm*2/c));
    
    Idata(stepNumber) =  Idata(stepNumber) + txInterfaceOil*(cos(-2*pi*propFreq*(((tank_h-dist_comm)*2/c_comm)+(dist_comm*2/c))));
    Qdata(stepNumber) =  Qdata(stepNumber) + txInterfaceOil*(sin(-2*pi*propFreq*(((tank_h-dist_comm)*2/c_comm)+(dist_comm*2/c))));
    


    end
%% Combine into one vector
IQ_data = Idata - 1i*Qdata;    
%% Add circulator coupling as DC offset
%rxsig = circulator(Circulator_Isolation,txsig,rxsig);
%% Received radar return with gain
rxsig = step(receiver,IQ_data);
%% Plot FFT and Range data
FFT_range(c,rxsig,freqStepSize,sweepNumber,Nsweep);
end   
    
    

 %% FilterDesigner LPF for filtering IQ jumps
 % Filter created in MATLAB filterDesigner to filter square jumps 
 % IQ_data: data to be filtered
 % Returns: IQ_data passed through a LPF
 function [filtered_data] = IQ_filter(IQ_data)
 % All frequency values are in Hz.
 Fs = 4e9;    % Sampling Frequency

 N  = 25;     % Order
 Fc = 450000; % Cutoff Frequency

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
 function [frequencyForCal] = generateSweepWaveform(BW,fc,freqSteps,call_dev,CALERROR,sweepType,tot_sweep_time)
    stepSize = BW/freqSteps;
    switch sweepType
        case 'up'
            for steps = 1:freqSteps
                   if CALERROR % simulate small calibration errors if bool is true
                       randomCallError = -call_dev + (call_dev).*rand(1,1);
                   else
                       randomCallError = 0;
                   end
               frequencyForCal(steps) = fc+(stepSize+randomCallError)*steps; % keep freqs
            end
        case 'down'
            for steps = 1:freqSteps
                   if CALERROR % simulate small calibration errors if bool is true
                       randomCallError = -call_dev + (call_dev).*rand(1,1);
                   else
                       randomCallError = 0;
                   end
             frequencyForCal(steps) = fc+(stepSize+randomCallError)*steps; % keep freqs
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
 % stepSizeHz: step size in Hz of each step
 % stepNumber: current step
 % Nsweep: number of sweeps
 % Returns: combined signal at all steps 
 function FFT_range (speedOfLight,IQ_data,Delta_F_Hz,stepNumber,Nsweep)
    %% Statistical data (Alex):
    % This stuff is populated with every sweep and processed later to obtain
    % standard deviations and average values.
    persistent peak_positions peak_magnitudes SNRs
    if stepNumber == 1
        peak_positions  = zeros(Nsweep, 1);
        peak_magnitudes = zeros(Nsweep, 1);
        SNRs            = zeros(Nsweep, 1);
    end

    B = IQ_data';
    L = length(B); % length of signal
    window = hann(L);
    B =  B.*window;
    L = round(speedOfLight/(Delta_F_Hz * 0.001) + 1); % res = 0.001 m
    B = [B; complex(zeros(L-length(IQ_data), 1))];
    Xaxis = 1:1:L;
    Xaxis = (((Xaxis*speedOfLight)/(Delta_F_Hz*(L - 1))));
    Xaxis = Xaxis./2;

    Y = ifft((B*L)); % undo the IFFT division by N and decimation division
    P2 = abs(Y);
    
    figure(2)
    hold on
    plot(Xaxis,mag2db(P2))
    axis([0 10 -40 30])
    title('Range Power Plot')
    xlabel('Range (m)')
    ylabel('|P1 db(m)|')
    title('SFCW IFFT Object Range and Magnitude');
   

    
    %% Estimate Range:
    [y,x] = max(mag2db(P2(round(1/Xaxis(1)):round(5/Xaxis(1))))); % find peak FFT point between 1m to 5m
    peak_x = Xaxis(x-1+round(1/Xaxis(1)));
    disp('Distance of object based on FFT (m):')
    disp(num2str(peak_x))
    
    %% Calculate SNR:
    window = 0.8;
    snr = calculateSNR(P2, Xaxis, peak_x, window);
    snr_disp = ['SNR: ', num2str(round(snr*100)/100), ' dB']; % round to 2 digits after decimal point
    disp(snr_disp)
    
    %% Store statistical data:
    peak_positions(stepNumber)  = peak_x;
    peak_magnitudes(stepNumber) = y;
    SNRs(stepNumber)            = snr;
    
    %% Output statistical data:
    if stepNumber == Nsweep
        % Everything is rounded to 2 digits after decimal point.
        peak_positions_std_dev = std(peak_positions)*1000;
        peak_magnitudes_std_dev = std(peak_magnitudes);
        SNRs_mean = round(mean(SNRs)*100)/100;
        stats_disp = ['Peak positions std dev: ', num2str(peak_positions_std_dev), ' mm', ...
                     ' Peak magnitudes std dev: ', num2str(peak_magnitudes_std_dev), ' dB', ...
                     ' Average SNR: ', num2str(SNRs_mean), ' dB'];
        disp(stats_disp)
    end
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
 % window (m)        : the range away from the peak to be used in avg noise power calculation
 % Returns: SNR in dB.
 function [snr_out] = calculateSNR(y_data, x_data, peak_position, window)
    % Temporary plot for debugging purposes:
    %figure(10)
    %plot(x_data, y_data)
    %hold on
    %axis([0 (peak_position + window) -1 19])
    
    %% Define window in terms of y_data indices:
    length_of_y        = size(y_data, 1);    
    peak_index         = round(peak_position/x_data(end)*length_of_y);
    window_index_left  = round(peak_index - (window/2)/x_data(end)*length_of_y);
    window_index_right = round(peak_index + (window/2)/x_data(end)*length_of_y);
    
    % Temporary console output for debugging purposes:
    debug_display = ['index_left: ', num2str(window_index_left), ' peak_index: ', num2str(peak_index), ' index_right: ', num2str(window_index_right)];
    %disp(debug_display)
    
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
    
    %% Output SNR in dB:
    snr_out = mag2db(y_data(peak_index)) - mag2db(noise_average);
 end
 
  %% Drift
 % drift the pre calibrated wave by a very small random frequency
 % divieation
 function [driftedCalFreq] = DriftCalibraton(driftError,frequencyForCal,CALERROR)
%% Create a sine wave for every step with the given number of points
    if CALERROR
    L = length(frequencyForCal);
    for steps = 1:L
       randomCallError = (-(rand(1,1)+rand(1,1)))/(1/driftError);
       frequencyReCal(steps) =  frequencyForCal(steps)+randomCallError*frequencyForCal(steps);
    end
    %figure(11);
    %hold on
    %plot(frequencyReCal);
    %plot(frequencyForCal);
    driftedCalFreq = frequencyReCal;
    else
    driftedCalFreq = frequencyForCal;
    end
    
 end
 
 
 
 

 