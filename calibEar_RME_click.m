%% Ear calibration following probe calibration
[FileName,PathName,FilterIndex] = uigetfile('Calib*RME*.mat',...
    'Please pick PROBE CALIBRATION file to use');
probefile = fullfile(PathName, FileName);
load(probefile);

% Initializing RME
InitializePsychSound(2);
% Open device in duplex mode (3), with full device control (2)
% Use 2 output channels and 1 input channel at sampling rate of Fs.
devices = PsychPortAudio('GetDevices');
deviceindex = [];
for k = 1:numel(devices)
    if strcmp(devices(k).DeviceName, 'ASIO MADIface USB')
        deviceindex = devices(k).DeviceIndex;
    end
end

Fs = calib.SamplingRate * 1000;
driver = calib.driver;
delay = 4174;
pahandle = PsychPortAudio('Open', deviceindex, 3, 2, Fs, [4, 4], calib.BufferSize, []);
% pahandle = PsychPortAudio('Open', deviceindex, 3, 2, Fs, [4, 4], calib.BufferSize,devices(deviceindex).LowInputLatency, []);

% Allocate input Buffer Size
PsychPortAudio('GetAudioData', pahandle, (calib.BufferSize + delay) * 1/Fs);

%% Get subject and ear info
subj = input('Please subject ID:', 's');
earflag = 1;
while earflag == 1
    ear = input('Please enter which ear (L or R):', 's');
    switch ear
        case {'L', 'R', 'l', 'r', 'Left', 'Right', 'left', 'right', 'LEFT', 'RIGHT'}
            earname = strcat('Ear-', upper(ear(1)));
            earflag = 0;
        otherwise
            fprintf(2, 'Unrecognized ear type! Try again!');
    end
end

calib.subj = subj;
calib.ear = ear;

% Make click
vo = clickStimulus(calib.BufferSize);
buffdata = zeros(4, calib.BufferSize+delay);
buffdata(driver, 1:length(vo)) = vo; % The other source plays nothing

% Check for clipping and load to buffer
if(any(abs(buffdata(driver, :)) > 1))
    error('What did you do!? Sound is clipping!! Cannot Continue!!\n');
end

% Fill the audio playback buffer with the audio data 'wavedata':
PsychPortAudio('FillBuffer', pahandle, buffdata);

playrecTrigger = 1;
resplength = numel(vo)+delay; % How many samples to read from OAE buffer

%% Set attenuation and play
drop = db2mag(-1 * calib.Attenuation);
PsychPortAudio('Volume', pahandle, drop);

vins_ear = zeros(calib.Averages, calib.BufferSize);
for n = 1: (calib.Averages + calib.ThrowAway)
    %Start playing from the buffer:
    % 1 Repetition (repeat handled in this script rather than
    % PsychPortAudio)
    startTime = PsychPortAudio('Start', pahandle, 1);
    % Stop playback:Audio('Start', pahandle, 1);
    WaitSecs(calib.BufferSize * 1/Fs);
    PsychPortAudio('Stop', pahandle);
    vin = PsychPortAudio('GetAudioData', pahandle);
    
    %vin = vin(3,:); %Tim broke this
    %Accumulate the time waveform - no artifact rejection
    if (n > calib.ThrowAway) && ~isempty(vin)
        vins_ear(n, :) = vin(3,:);
    end
end

energy = squeeze(sum(vins_ear.^2, 2));
good = (energy < median(energy) + 2*mad(energy)) & (energy > median(energy) - 2*mad(energy)) ;
vavg = squeeze(mean(vins_ear(good, :), 1));
Vavg = rfft(vavg');
calib.vavg_ear = vavg;

% Apply calibrations to convert voltage to pressure
% For ER-10X, this is approximate
mic_sens = 50e-3; % mV/Pa. TO DO: change after calibration
mic_gain = db2mag(40);
P_ref = 20e-6;
DR_onesided = 1;
mic_output_V = Vavg / (DR_onesided * mic_gain);
output_Pa = mic_output_V/mic_sens;
output_Pa_20uPa_per_Vpp = output_Pa / P_ref; % unit: 20 uPa / Vpeak-peak

freq = 1000*linspace(0,calib.SamplingRate/2,length(Vavg))';
calib.vins_ear = vins_ear;

% Note no attenuation gives 4.75 V peak for the chirp
Vo = rfft(calib.vo) * 4.75 * db2mag(-1 * calib.Attenuation);

calib.EarRespH =  output_Pa_20uPa_per_Vpp ./ Vo; %save for later


PsychPortAudio('Close', pahandle);

%% Plot data
figure(1);
ax(1) = subplot(2, 1, 1);
semilogx(calib.freq, db(abs(calib.EarRespH)), 'linew', 2);
ylabel('Response (dB re: 20 \mu Pa / V_{peak})', 'FontSize', 16);
ax(2) = subplot(2, 1, 2);
semilogx(calib.freq, unwrap(angle(calib.EarRespH), [], 1), 'linew', 2);
xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('Phase (rad)', 'FontSize', 16);
linkaxes(ax, 'x');
legend('show');
xlim([100, 24e3]);

%% Calculate Ear properties
calib = findHalfWaveRes(calib);
calib.Zec_raw = ldimp(calib.Zs, calib.Ps, calib.EarRespH);
% calib.Zec = zsmo(calib.Zec, z_tube(calib.CavTemp, calib.CavDiam),...
%     calib.SamplingRate * 1000);
calib.Zec = calib.Zec_raw;

% decompose pressures
calib.fwb = 0.55;% bandwidth/Nyquist freq of freq.domain window

%%

% *ec: Ear canal
% *s: Source
% R*: Reflectance
% Z*: Impedance
% Pfor: Forward pressure
% Prev: Reverse pressure
% Pinc: Incident pressure

[calib.Rec, calib.Rs, calib.Rx, calib.Pfor, calib.Prev, calib.Pinc, ...
    calib.Px, calib.Z0, calib.Zi, calib.Zx] = decompose(calib.Zec, ...
    calib.Zs, calib.EarRespH, calib.Ps, calib.fwb, ...
    calib.CavTemp, calib.CavDiam);

% Check for leaks as in Groon et al
ok = find (calib.freq >= 200 & calib.freq <= 500);
calib.A_lf =  mean(1-(abs(calib.Rec(ok))).^2);
fprintf(1, 'Low-frequency absorbance: %2.3f\n', calib.A_lf);
calib.Yphase_lf = mean(cycs(1./calib.Zec(ok)))*360;
fprintf(1, 'Low-frequency admittance phase: %2.3f%c\n',...
    calib.Yphase_lf, char(176));

if (calib.A_lf > 0.29)
    h = warndlg ('Sound-leak alert! Low-frequency absorbance > 0.29');
    waitfor(h);
end

if (calib.Yphase_lf < 44)
    h = warndlg ('Sound-leak alert! Low-frequency admittance phase < 44 degrees');
    waitfor(h);
end

%% Plot Ear Absorbance
figure(2);
semilogx(calib.freq * 1e-3, 100*(1 - abs(calib.Rec).^2), 'linew', 2);
xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('Absorbance (%)', 'FontSize', 16);
xlim([0.2, 8]); ylim([0, 100]);
set(gca, 'FontSize', 16, 'XTick',[0.25, 0.5, 1, 2, 4, 8]);

%% Save Ear Calculations
datetag = datestr(clock);
calib.date = datetag;
datetag(strfind(datetag,' ')) = '_';
datetag(strfind(datetag,':')) = '_';
fname = strcat('Calib_',calib.drivername,calib.device,'_',subj,earname,'_',date, '_RME.mat');
save(fname,'calib');

% just before the subject arrives
