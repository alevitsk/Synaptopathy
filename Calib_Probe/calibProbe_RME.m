%% Sound Source, Microphone Probe Thevenin Calibration
% Note: Calibrations need to be run separately for each sound source

% Initialize ER-10X  (Also needed for ER-10C for calibrator)
p10x = genpath('C:\Users\Lab User\Desktop\Experiments\Calibration\ER10X\Matlab');
addpath(p10x);
loaded = ER10XConnectLoadLib('C:\Users\Lab User\Desktop\Experiments\Calibration\ER10X\ER10X API\');
[err, ER10XHandle] = er10x_open();
fprintf(1, 'Result of ER10_X_OPEN: %s\n', err);
if strcmp(err, 'ER10X_ERR_OK')
    fprintf('Continuing...\n');
else
    error('Something wrong! Calibration aborted!');
end
err = er10x_connect(ER10XHandle);
fprintf(1, 'Result of ER10_X_CONNECT: %s\n', err);
if strcmp(err, 'ER10X_ERR_OK')
    fprintf('Continuing...\n');
else
    error('Something wrong! Calibration aborted!');
end

% Initializing Calibration
calib = calibSetDefaults_RME();
Fs = calib.SamplingRate * 1000;

%% Added to account for the installation of the new PTB version.
% addPTB = 1; % if 1, then assume PTB isn't in the path and add it
% ptbPath = 'C:\oldtoolbox\Psychtoolbox';
%
% if addPTB
%     if ~exist(ptbPath,'dir')
%         fprintf(['\n\nCouldn''t find the specified PsychToolBox directory: %s.\n'...
%             'Please adjust the ptbPath variable in the %s function to point to your PTB directory.\n'], ptbPath, mfilename);
%         return;
%     end
%
%     addpath(genpath(ptbPath));
%     rmpath('C:\toolbox\Psychtoolbox');
%     PsychtoolboxPostInstallRoutine(1); % run a slightly edited version of the PTB post-install routine
% end

% Initializing RME
InitializePsychSound(1);
% Open device in duplex mode (3), with full device control (2)
% Use 2 output channels and 1 input channel at sampling rate of Fs.
devices = PsychPortAudio('GetDevices');
deviceindex = [];
for k = 1:numel(devices)
    if strcmp(devices(k).DeviceName, 'ASIO MADIface USB')
        deviceindex = devices(k).DeviceIndex;
    end
end

% Channel index in PyschPortAudio starts at zero whereas the RME box
% numbering starts at 1.
% So 2, 3 here would be the "balanced line outs" #3 and #4 respectively
% Zero on the second row corresponds to "balanced line in" channel #1.
delay = 4174;
pahandle = PsychPortAudio('Open', deviceindex, 3, 2, Fs, [4, 4], calib.BufferSize, []);

% Allocate input Buffer Size
PsychPortAudio('GetAudioData', pahandle, (calib.BufferSize+delay) * 1/Fs * (calib.Averages + calib.ThrowAway));

deviceflag = 1;
while deviceflag == 1
    device = input('Please enter X or C for ER-10X/ER-10C respectively:', 's');
    switch device
        case {'X', 'x'}
            device = 'ER-10X';
            deviceflag = 0;
        case {'C', 'c'}
            device = 'ER-10C';
            deviceflag = 0;
            % ER-10C has more distortion, hence attenuate by another 15 dB
            calib.Attenuation = calib.Attenuation + 15;
        otherwise
            fprintf(2, 'Unrecognized device! Try again!');
    end
end

driverflag = 1;
while driverflag == 1
    driver = input('Please enter whether you want driver 1, 2, or 3 (Aux on ER-10X):');
    switch driver
        case {1, 2}
            drivername = strcat('Ph',num2str(driver));
            driverflag = 0;
        case 3
            if strcmp(device, 'ER-10X')
                drivername = 'PhAux';
                driverflag = 0;
            else
                fprintf(2, 'Unrecognized driver! Try again!');
            end
        otherwise
            fprintf(2, 'Unrecognized driver! Try again!');
    end
end
calib.device = device;
calib.drivername = drivername;
calib.driver = driver;

% Make linear chirp upto nyquist
vo = chirpStimulus(calib.BufferSize, 0.90);
buffdata = zeros(4, numel(vo)+delay);
buffdata(driver, 1:length(vo)) = vo; % The other source plays nothing

% Check for clipping and load to buffer
if(any(abs(buffdata(1, :)) > 1))
    error('What did you do!? Sound is clipping!! Cannot Continue!!\n');
end

% Fill the audio playback buffer with the audio data 'wavedata':
PsychPortAudio('FillBuffer', pahandle, buffdata);

playrecTrigger = 1;
resplength = numel(vo)+delay; % How many samples to read from OAE buffer


%% Set attenuation and play
calib.vo = vo;
calib.vins = zeros(calib.CavNumb, calib.Averages, calib.BufferSize);
% drop = db2mag(-1 * calib.Attenuation);
drop = db2mag(-1 * 39);
PsychPortAudio('Volume', pahandle, drop);

err = er10x_move_to_position_and_wait(ER10XHandle, 0, 20000);
fprintf(1, 'Result of moving to position 1: %s\n', err);
if strcmp(err, 'ER10X_ERR_OK')
    fprintf('Continuing...\n');
else
    error('Something wrong! Calibration aborted!');
end

for m = 1:calib.CavNumb
    %     %Play chirp and measure response
    %     for n = 1: (calib.Averages + calib.ThrowAway)
    %         %Start playing from the buffer:
    %         % 1 Repetition (repeat handled in this script rather than
    %         % PsychPortAudio)
    %         startTime = PsychPortAudio('Start', pahandle, 1);
    %         WaitSecs(calib.BufferSize * 1/Fs);
    %         % Stop playback:
    %         PsychPortAudio('Stop', pahandle);
    %         vin = PsychPortAudio('GetAudioData', pahandle);
    %         %Accumluate the time waveform - no artifact rejection
    %         if (n > calib.ThrowAway)
    %             calib.vins(m, n-calib.ThrowAway, :) = vin;
    %         end
    %     end
    
    
    %Start playing from the buffer:
    startTime = PsychPortAudio('Start', pahandle, calib.Averages + calib.ThrowAway);
%     WaitSecs(calib.BufferSize * 1/Fs * (calib.Averages + calib.ThrowAway));
    WaitSecs(length(buffdata)* 1/Fs * (calib.Averages + calib.ThrowAway));
    % Stop playback:
    PsychPortAudio('Stop', pahandle, 1);
    vin = PsychPortAudio('GetAudioData', pahandle);
    vin = vin(3,:); %Tim broke this
    
    %Compute the average
    vin_reshaped = reshape(vin,...
        (calib.BufferSize+delay), (calib.Averages + calib.ThrowAway))';
    calib.vins(m, :, :) = vin_reshaped( (calib.ThrowAway + 1):end, delay+1:end);
    energy = squeeze(sum(calib.vins(m, :, :).^2, 3));
    good = energy < median(energy) + 2*mad(energy);
    vavg = squeeze(mean(calib.vins(m, good, :), 2));
    Vavg = rfft(vavg);% ./ rfft(vo);
    
    % Apply calibrations to convert voltage to pressure
    % For ER-10X, this is approximate
    mic_sens = 50e-3; % mV/Pa. TO DO: change after calibration
%         mic_gain = db2mag(ER10X_info.gain+6);
    mic_gain = db2mag(40);
    P_ref = 20e-6;
    DR_onesided = 1;
    mic_output_V = Vavg / (DR_onesided * mic_gain);
    output_Pa = mic_output_V/mic_sens;
    output_Pa_20uPa_per_Vpp = output_Pa / P_ref; % unit: 20 uPa / Vpeak-peak
    
    freq = 1000*linspace(0,calib.SamplingRate/2,length(Vavg))';
    calib.freq = freq;
    
    % Note no attenuation gives 4.75 V peak for the chirp
    Vo = rfft(calib.vo) * 4.75 * db2mag(-1 * calib.Attenuation);
    
    calib.CavRespH(:,m) =  output_Pa_20uPa_per_Vpp ./ Vo; % unit = 20 uPa/Vpp
    calib.unit = '20 uPa / 1 Vpeak';
    
    if m < calib.CavNumb
        err = er10x_move_to_position_and_wait(ER10XHandle, m, 20000);
        fprintf(1, 'Result of moving to position %d: %s\n', m+1, err);
        if strcmp(err, 'ER10X_ERR_OK')
            fprintf('Continuing...\n');
        else
            error('Something wrong! Calibration aborted!');
        end
        
    else
        if (calib.doInfResp == 1)
            out2 = input(['Done with ER-10X cavities.. Move to infinite tube!\n',...
                'Continue? Press n to stop or any other key to go on:'], 's');
        end
    end
end

if(calib.doInfResp == 1)
    % One more measurement in infinite pipe
    vins_inf = zeros(calib.Averages, calib.BufferSize);
    for n = 1: (calib.Averages + calib.ThrowAway)
        %Start playing from the buffer:
        invoke(RZ, 'SoftTrg', playrecTrigger);
        currindex = invoke(RZ, 'GetTagVal', 'indexin');
        while(currindex < resplength)
            currindex=invoke(RZ, 'GetTagVal', 'indexin');
        end
        
        vin = invoke(RZ, 'ReadTagVex', 'dataout', 0, resplength,...
            'F32','F64',1);
        %Accumluate the time waveform - no artifact rejection
        if (n > calib.ThrowAway)
            vins_inf(n, :) = vin;
        end
        % Get ready for next trial
        invoke(RZ, 'SoftTrg', 8); % Stop and clear "OAE" buffer
        %Reset the play index to zero:
        invoke(RZ, 'SoftTrg', 5); %Reset Trigger
        
    end
    energy = squeeze(sum(vins_inf.^2, 2));
    good = energy < median(energy) + 2*mad(energy);
    vavg = squeeze(mean(vins_inf(good, :), 1));
    Vavg = rfft(vavg')';
    
    
    % Apply calibrations to convert voltage to pressure
    % For ER-10X, this is approximate
    mic_sens = 50e-3; % mV/Pa. TO DO: change after calibration
    mic_gain = db2mag(40);
    P_ref = 20e-6;
    DR_onesided = 1;
    mic_output_V = Vavg / (DR_onesided * mic_gain);
    output_Pa = mic_output_V/mic_sens;
    output_Pa_20uPa_per_Vpp = output_Pa / P_ref; % unit: 20 uPa
    
    freq = 1000*linspace(0,calib.SamplingRate/2,length(Vavg))';
    calib.freq = freq;
    
    % Note no attenuation gives 4.75 V peak for the chirp
    Vo = rfft(calib.vo) * 4.75 * db2mag(-1 * calib.Attenuation);
    
    calib.InfRespH = output_Pa_20uPa_per_Vpp ./ Vo; % unit = 20 uPa/Vpp
end


%% Plot data
figure(1);
ax(1) = subplot(2, 1, 1);
semilogx(calib.freq, db(abs(calib.CavRespH)), 'linew', 2);
ylabel('Response (dB re: 20 \mu Pa / V_{peak})', 'FontSize', 16);
ax(2) = subplot(2, 1, 2);
semilogx(calib.freq, unwrap(angle(calib.CavRespH), [], 1), 'linew', 2);
xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('Phase (rad)', 'FontSize', 16);
linkaxes(ax, 'x');
legend('show');
xlim([100, 24e3]);
%% Compute Thevenin Equivalent Pressure and Impedance

%set up some variables
irr = 1; %ideal cavity reflection

%  calc the cavity length
calib.CavLength = cavlen(calib.SamplingRate,calib.CavRespH, calib.CavTemp);
if (irr)
    la = [calib.CavLength 1]; %the one is reflection fo perfect cavit
else
    la = calib.CavLength; %#ok<UNRCH>
end

df=freq(2)-freq(1);
jef1=1+round(calib.f_err(1)*1000/df);
jef2=1+round(calib.f_err(2)*1000/df);
ej=jef1:jef2; %limit freq range for error calc

calib.Zc = cavimp(freq, la, irr, calib.CavDiam, calib.CavTemp); %calc cavity impedances

%% Plot impedances
% It's best to have the set of half-wave resonant peaks (combined across
% all cavities and including all harmonics) distributed as uniformly as
% possible across the frequency range of interest.
figure(2)
plot(freq/1000,dB(calib.Zc)); hold on
xlabel('Frequency kHz')
ylabel('Impedance dB')
%
pcav = calib.CavRespH;
options = optimset('TolFun', 1e-12, 'MaxIter', 1e5, 'MaxFunEvals', 1e5);
la=fminsearch(@ (la) thverr(la,ej, freq, pcav, irr, calib.CavDiam, calib.CavTemp),la, options);
calib.Error = thverr(la, ej, freq, pcav, irr, calib.CavDiam, calib.CavTemp);

calib.Zc=cavimp(freq,la, irr, calib.CavDiam, calib.CavTemp);  % calculate cavity impedances
[calib.Zs,calib.Ps]=thvsrc(calib.Zc,pcav); % estimate zs & ps

plot(freq/1000,dB(calib.Zc),'--'); %plot estimated Zc

calib.CavLength = la;

if ~(calib.Error >= 0 && calib.Error <=1)
    h = warndlg ('Calibration error out of range!');
    waitfor(h);
end

%% Save calib.Zs and Ps - you can measure them weekly/daily and load
datetag = datestr(clock);
calib.date = datetag;
datetag(strfind(datetag,' ')) = '_';
datetag(strfind(datetag,':')) = '_';
fname = strcat('Calib_',drivername,device,datetag, '_RME.mat');
save(fname,'calib');
%% Close ER-10X, psychPortAudio connections etc.
PsychPortAudio('Close', pahandle);
err = er10x_disconnect(ER10XHandle);
if strcmp(err, 'ER10X_ERR_OK')
    fprintf('Continuing...\n');
else
    error('Something wrong! Could not close ER10X!');
end
[err, ER10XHandle] = er10x_close(ER10XHandle);
if strcmp(err, 'ER10X_ERR_OK')
    fprintf('Continuing...\n');
else
    error('Something wrong! Could not close ER10X!');
end
ER10XConnectUnloadLib()
rmpath(p10x);
% just before the subject arrives