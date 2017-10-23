function varargout = AudioEffectsGenerationFinal(varargin)
% AUDIOEFFECTSGENERATIONFINAL MATLAB code for AudioEffectsGenerationFinal.fig
%      AUDIOEFFECTSGENERATIONFINAL, by itself, creates a new AUDIOEFFECTSGENERATIONFINAL or raises the existing
%      singleton*.
%
%      H = AUDIOEFFECTSGENERATIONFINAL returns the handle to a new AUDIOEFFECTSGENERATIONFINAL or the handle to
%      the existing singleton*.
%
%      AUDIOEFFECTSGENERATIONFINAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUDIOEFFECTSGENERATIONFINAL.M with the given input arguments.
%
%      AUDIOEFFECTSGENERATIONFINAL('Property','Value',...) creates a new AUDIOEFFECTSGENERATIONFINAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AudioEffectsGenerationFinal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AudioEffectsGenerationFinal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AudioEffectsGenerationFinal

% Last Modified by GUIDE v2.5 27-Nov-2013 17:16:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AudioEffectsGenerationFinal_OpeningFcn, ...
                   'gui_OutputFcn',  @AudioEffectsGenerationFinal_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before AudioEffectsGenerationFinal is made visible.
function AudioEffectsGenerationFinal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AudioEffectsGenerationFinal (see VARARGIN)

% Choose default command line output for AudioEffectsGenerationFinal
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AudioEffectsGenerationFinal wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AudioEffectsGenerationFinal_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnOpenWaveFile.
function btnOpenWaveFile_Callback(hObject, eventdata, handles)
% hObject    handle to btnOpenWaveFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname]= uigetfile({'*.wav'},'File Selector');
pathname=strcat(pathname,filename);
set(handles.txtFileName,'String',filename);
[signal, Fs] = wavread(pathname);
set(handles.txtSamplingRate,'String',Fs);
set(handles.txtPath,'String',pathname);


% --- Executes on button press in btnPlayOriginal.
function btnPlayOriginal_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlayOriginal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x=get(handles.txtPath,'String');
[signal, Fs]=wavread(x);
axes(handles.axsOriginalSignalAT);
T = 1/Fs;                     % Sample time
L=size(signal,1);
t = (0:L-1)*T;                % Time vector
plot(t,signal)
axis tight;
xlabel('Time (s)');
ylabel('Amplitude');
axes(handles.axsOriginalSignalFS);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
fftsignal = fft(signal,NFFT)/L;
xpts = Fs/2*linspace(0,1,NFFT/2+1);
plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
axis tight;
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
wavplay(signal,Fs);

% --- Executes on button press in btnPlaySoundEffect.
function btnPlaySoundEffect_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlaySoundEffect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    contents = get(handles.popSoundEffect,'Value');
    switch contents
        case 1
            x = get(handles.txtPath,'String');
            [signal, Fs] = wavread(x);
            y = equalizer(signal,Fs);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);    
        case 2
            x = get(handles.txtPath,'String');
            [signal, Fs] = wavread(x);
            delay = 0.06;  % 60ms dealy
            gain = 0.7;    
            y = reverberation(signal,Fs,gain,delay);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);    
        case 3
            x = get(handles.txtPath,'String');
            [signal, Fs] = wavread(x);
            f = 10; % frequency in Hz
            y=Tremolo(signal,Fs,f);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);        
        case 4
            x = get(handles.txtPath,'String');
            [signal, Fs] = wavread(x);
            Modfreq=10;     % 10 KHz
            Width=0.0008;   % 8 Ms
            y = vibrato(signal,Fs,Modfreq,Width);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);        
        case 5
            x = get(handles.txtPath,'String');
            [signal, Fs] = wavread(x);
            y = chorus(signal,Fs);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);        
        case 6
            x = get(handles.txtPath,'String');
            [signal, Fs] = wavread(x);
            delay = 250;
            n = 10;
            y = tunnelecho(signal,Fs,delay,n);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);        
        case 7
            x = get(handles.txtPath,'String');
            [signal, Fs] = wavread(x);
            delay = 0.003;   % 3ms delay
            y = flanging(signal,Fs,delay);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);        
        case 8
            x=get(handles.txtPath,'String');
            [signal, Fs]=wavread(x);
            delay=100;
            n=2;
            y=slapback(signal,Fs,delay,n);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);
        case 9 
            x=get(handles.txtPath,'String');
            [signal, Fs]=wavread(x);
            y = decimation(signal);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);
        case 10
            x=get(handles.txtPath,'String');
            [signal, Fs]=wavread(x);
            a = 5; % amount of fuzz in signal
            y = fuzz(signal,a);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);
        case 11
            x=get(handles.txtPath,'String');
            [signal, Fs]=wavread(x);
            y = interpolation(signal);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);
        case 12
            x=get(handles.txtPath,'String');
            [signal, Fs]=wavread(x);
            num_sig = 3;
            delay = 500;
            y = openecho(signal, Fs, delay, num_sig);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);
        case 13
            x=get(handles.txtPath,'String');
            [signal, Fs]=wavread(x);
            y = overdrive(signal);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);   
        case 14
            x=get(handles.txtPath,'String');
            [signal, Fs]=wavread(x);
            f = 300; % frequency 300 Hz
            y = ringModulation(signal, Fs, f);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);
        case 15
             x=get(handles.txtPath,'String');
            [signal, Fs]=wavread(x);
            y = reverseecho(signal, Fs);
            axes(handles.axsSoundEffectAT);
            T = 1/Fs;                     % Sample time
            L=size(y,1);
            t = (0:L-1)*T;                % Time vector
            plot(t,y)
            axis tight;
            xlabel('Time (s)');
            ylabel('Amplitude');
            axes(handles.axsSoundEffectFS);
            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
            fftsignal = fft(y,NFFT)/L;
            xpts = Fs/2*linspace(0,1,NFFT/2+1);
            plot(xpts,2*abs(fftsignal(1:NFFT/2+1)));
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('|Y(f)|');
            wavplay(y,Fs);
        otherwise  
    end


%Frequency band filtering/ equalizer
function y = equalizer(signal,Fs)
    filt_gain = [6 3 2 1.7 1.3 1 0.8 0.6];
    
    [b,a] = butter(2,200/(Fs/2),'low');
    y1 = filt_gain(1).*filter(b,a,signal);

    [b,a] = butter(2,[(200/(Fs/2)) (500/(Fs/2))],'bandpass');
    y2 = filt_gain(2).*filter(b,a,signal);

    [b,a] = butter(2,[(500/(Fs/2)) (1000/(Fs/2))],'bandpass');
    y3 = filt_gain(3).*filter(b,a,signal);

    [b,a] = butter(2,[(1000/(Fs/2)) (2000/(Fs/2))],'bandpass');
    y4 = filt_gain(4).*filter(b,a,signal);

    [b,a] = butter(2,[(2000/(Fs/2)) (4000/(Fs/2))],'bandpass');
    y5 = filt_gain(5).*filter(b,a,signal);

    [b,a] = butter(2,[(4000/(Fs/2)) (7000/(Fs/2))],'bandpass');
    y6 = filt_gain(6).*filter(b,a,signal);

    [b,a] = butter(2,[(7000/(Fs/2)) (12000/(Fs/2))],'bandpass');
    y7 = filt_gain(7).*filter(b,a,signal);

    [b,a] = butter(2,[(12000/(Fs/2)) (20000/(Fs/2))],'bandpass');
    y8 = filt_gain(8).*filter(b,a,signal);
    
    y = y1+y2+y3+y4+y5+y6+y7+y8;


%Reverberation
function y = reverberation(signal, Fs, gain, delay_time)

len=length(signal);
 
delay = round(Fs*delay_time);

y = signal;

for i=(delay+1):1:len
    y(i,1) = (gain*signal(i,1)) + signal((i - delay),1) - (gain * y((i - delay),1));
end
    
%Tremolo
function y = Tremolo(signal, Fs, freq)
    len=length(signal);

    t=1:len;
    modulator = sin(2*pi*t*freq/Fs);

    y = signal.*modulator';


%vibrato
function y=vibrato(signal,Fs,freq,Width)
    Delay=Width; % basic delay of input sample in sec
    DELAY=round(Delay*Fs); % basic delay in # samples
    WIDTH=round(Width*Fs); % modulation width in # samples

    f=freq/Fs; % modulation frequency in # samples
    len=length(signal); % # of samples in WAV-file
    l=2+DELAY+WIDTH*2; % length of the entire delay
    Delayline=zeros(l,1); % memory allocation for delay
    y=zeros(size(signal)); % memory allocation for output vectorCM0268

    for n=1:(len-1)
        z=1+DELAY+WIDTH*sin(2*pi*n*f);
        i=floor(z);
        frac=z-i;
        Delayline=[signal(n);Delayline(1:l-1)];
        
        y(n,1)=Delayline(i+1)*frac+Delayline(i)*(1-frac);
    end    
%refernce: http://www.cs.cf.ac.uk/Dave/CM0268/PDF/10_CM0268_Audio_FX.pdf
%Page 376
 
    
%Chorus
function y = chorus(signal, Fs)
    len=length(signal);
    num_sigs = 5;

    rng('shuffle');
    delay_time = randi([10 25],num_sigs,1);
    delay = round(Fs.*(delay_time./1000));
    sum_delays = sum(delay);

    x=zeros(num_sigs,len);

    for i=1:num_sigs
        x(i,:) = signal(:,1)';
    end

    y=zeros(len+(sum_delays*(num_sigs-1)),1);
    amp = 1.0;
    for i=1:num_sigs
        y(1+delay(i)*(i-1):len+delay(i)*(i-1),1) = amp.*(y(1+delay(i)*(i-1):len+delay(i)*(i-1),1)+x(i,:)');
    end
    

%Tunnelecho
function y = tunnelecho(signal,Fs,delay_time,n)
    len=length(signal);

    delay = round(Fs*(delay_time/1000));

    x=zeros(n,len);

    for i=1:n
        x(i,:) = signal(:,1)';
    end

    y=zeros(len+(delay*(n-1)),1);
    amp = 1.0;
    for i=1:n
        if amp<0.1
            amp = 0.01;
        end
        y(1+delay*(i-1):len+delay*(i-1),1) = amp.*(y(1+delay*(i-1):len+delay*(i-1),1)+x(i,:)');
        amp = amp - (1.0/n);
    end
    
    
%Flanging
function y = flanging(signal, Fs, delay_time)

    len=length(signal);

    samp_delay = round(Fs*delay_time);

    t=1:1:len;
    modulator = sin(2*pi*t/Fs);

    y(1:samp_delay,1) = signal(1:samp_delay,1);

    for i=(samp_delay+1):1:len
        delay = abs(round(modulator(i)*samp_delay));
        y(i,1) = signal(i,1) + signal((i - delay),1);
    end
%refernce: http://www.cs.cf.ac.uk/Dave/CM0268/PDF/10_CM0268_Audio_FX.pdf
%Page 382
    
    
%Slapback
function y = slapback(signal, Fs, delay_time, n)

    len=length(signal);
    delay = round(Fs*(delay_time/1000));
    x=zeros(n,len);

    for i = 1:n
        x(i,:) = signal(:,1)';
    end

    y=zeros(len+(delay*(n-1)),1);
    amp = 1.0;
    for i=1:n
        disp(amp);
        y(1+delay*(i-1):len+delay*(i-1),1) = amp.*(y(1+delay*(i-1):len+delay*(i-1),1)+x(i,:)');
    end


%Decimation
function y = decimation(signal)

    len = length(signal);
    y = zeros(len/2,1);
    for i = 1:(len/2)
        k = 2*i;
        value = signal(k,1);
        y(i,1) = value;
    end
    

%Fuzz
function y = fuzz(signal, alpha)

    len=length(signal);
    y = signal;
    for i=1:1:len
        y(i,1) = (signal(i,1)./abs(signal(i,1))).*(1-exp(alpha.*signal(i,1).*signal(i,1)./abs(signal(i,1))));
    end
   
    
%Interpolation
function y = interpolation(signal)

    len = length(signal);
    y = zeros(len*2,1);
    for i =1:len
        k = 2*i;
        value = signal(i,1);
        y(k,1) = value;
    end
    

%Openecho
function y = openecho(signal, Fs, delay_time, n)

    len=length(signal);
    delay = round(Fs*(delay_time/1000));
    x=zeros(n,len);

    for i=1:n
        x(i,:) = signal(:,1)';
    end

    y=zeros(len+(delay*(n-1)),1);
    amp = [3.2 1.7 0.5];
    for i=1:n
        y(1+delay*(i-1):len+delay*(i-1),1) = amp(i).*(y(1+delay*(i-1):len+delay*(i-1),1)+x(i,:)');
    end
    
    
% Overdrive
function y = overdrive(signal)

    len=length(signal);
    y = signal;
    %f(x) = 2x, for 0<=x<1/3
    %f(x) = (3-(2-3x)^2)/3 for 1/3<=x<2/3
    %f(x) = 1 for 2/3<=x<=1
    for i = 1:len
       curr_input = abs(signal(i,1));

       if curr_input >=0 && curr_input<(1/3)
           y(i,1) = 2*curr_input;
       end

       if curr_input >=(1/3) && curr_input<(2/3)
           y(i,1) = (3-(2-(3*curr_input))^2)/3;
       end

       if curr_input >=(2/3) && curr_input<=1
           y(i,1) = 1;
       end
    end
    

%Ring Modulation
function y = ringModulation(signal, Fs, freq)

    len=length(signal);
    t=1:len;
    modulator = sin(2*pi*t*freq/Fs);
    y = signal(:,1).*modulator';
    
%ReverseEcho   
    function [ outputwav ] = reverseecho( inputwav, Fs)
%RECHO is a function that creates a reverse echo effect on an input .wav file
%   inputwav is a filename of a .wav file and must be a character string
%   nechos is the number of reverse echos 
%   delay is the amount of delay in milliseconds and should be at least 35
%   build is the attentuation of the signal at each echo i.e. the first
%       reverse echo would be build^necho*original signal
%The output is of a signal that can be played with Matlab's sound function
%and read into wavwrite with sampling frequency Fs equal to input signal in
%order to create the output .wav file.

build = 0.2;
delay = 80;
necho = 4;

y=inputwav;
length=size(y,1);

reversey=zeros(length,1);
for i=0:length-1
    reversey(i+1)=y(length-i);
end

delaysample = round(Fs*delay/1000);
n=round(length+necho*delaysample);

routputwav=zeros(n,1);
routputwav(1:length)=reversey(:);

for i=1:necho;
    start=i*delaysample;
    stop=start+length-1;
    for j=start:stop
    routputwav(j,1)=routputwav(j,1) + build*reversey(j-start+1,1);
    end
end

outputwav=zeros(n,1);
for i=0:n-1;
    outputwav(i+1)=routputwav(n-i);
end


% --- Executes on selection change in popSoundEffect.
function popSoundEffect_Callback(hObject, eventdata, handles)
% hObject    handle to popSoundEffect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popSoundEffect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popSoundEffect

    

% --- Executes during object creation, after setting all properties.
function popSoundEffect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popSoundEffect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function txtFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in btnPlaySoundEffect.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlaySoundEffect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
