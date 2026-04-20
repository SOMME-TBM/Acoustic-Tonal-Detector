function [raw_tonals,smooth_tonals,energetic_detector,energetic_duration_detector] = func_tonal_detector(sig_acou,fs,dt_aPSD,bandwidth_call_freq,window_call_durations,space_btw_diff_call,thresh_sigma,appli_fit_type,fit_type,thresh_Rsquare)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tonal detector - description
% This tonal detector finds energy peaks over the signal and then assessed  tonal contour.

% The first step is an energy detector, which estimates the average sound
% level exposure (aPSD) as the power spectral density level integrated over a
% frequency bandwidths (bandwidth_call_freq) for a given duration (dt_aPSD)
% The detector considers the distribution of aPSD values over the signal
% with the hypothesis that the distribtuion should be normal, and thus
% based upon the mode, median or mean of aPSD values set a threshold above
% which values are considered as overdisperes (i.e. few occurence of high
% energy):

% First step doesn't discriminate energetic sounds, thus a first
% selection is based upon the tonal durations with a minimum and a maximum
% duration possible for the tonals (window_call_durations)

% The second step is a contour-based detector.
% The contour (i.e. frequency) is recovered from the spectrogram of the
% potential tonals detected based upon the energy criteria, by finding which
% frequency has the maximum energy for every time unit of the spectrogram.
% From tonals definition, the contour should be smooth and thus a fit model
% could be applied (fit_type). Then, a R² threshold is set to considered
% whether the energetic sound follows accurately a descriptive fit model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see Caron Delbosc et al. 2025 for more details and please cite this paper
% if you use this function
% for any question or bug reporting please contact:
% Gaëtan RICHARD richard.somme@orange.fr
% or Naïs CARON DELBOSC : nais.caron.delbosc@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detector output
% 'raw_tonals' is a cell array with all contours detected and extracted:
% every cell (tonal) is a double-precision array [time frequency]
% the time is in second from the begininig of the signal (sig_acou)

% 'smooth_tonals' is the same as raw_tonals but with smoothed contours (using a pchip relationship)

% 'energetic_detector' is a matrix (type double) with two column:
% First col: time in second when the average PSD was found above a threshold (resolution of this column is at 'dt_aPSD' seconds)
% Second col: numer of the event (several line with the same number = same event

% 'energetic_duration_detector' is a matrix (type double) with three column derived from 'energetic_detector' by considering onmy events long enough and not too long (i.e within [window_call_durations]):
% 1st col: begining time in second of the event kept from 'energetic_detector'
% 2nd col: end time in second of the event kept from 'energetic_detector'
% 3rd col: number of the event (kept from 'energetic_detector')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detector parameters  - INPUTS

% sig_acou: acoustic signal
% fs: sampling rate of the acoustic signal

% dt_aPSD : aPSD resoloution => time window over which thre PSD is estimated to then average between the bandwidth frequency call of interest:
% bandwidth_call_freq=[min tonal frequency; max tonal frequency]

% window_call_durations=[min tonal duraitn; max tonal duration]
% space_btw_diff_call= duration in sec between two energetic events to consider new call

% thresh_sigma= best to set it at 2 or 3,
% for instance thresh_sigma= 3 means  [−3 σ; 3 σ] with σ being the standard deviation of the distribution
% here, for a random variable that follows a normal distribution, 99.74% of the values (corresponding
% to ambient noise within our framework) belong to the [−3 σ; 3 σ] interval in comparison with the mean

% appli_fit_type : can be set to 'no' to apply only the energetic detector,
% otehrwise can be set to 'yes' or left empty
% fit_type = 'poly3','exp1', 'exp2'... see "help fit"
% thresh_Rsquare: minimum R² to consider fit ok

% 'threshold_detect' is estimated usin the mode, but can be change to median or mean if
% requires, but mode seemed to be the best approach.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example of inputs

% dt_aPSD=0.1; % integration time for aPSD => resolution
% window_call_durations=[1 120]; % in sec, ex for bearded seals : at least one sec and could reach 2minutes
% space_btw_diff_call=1;% difference in sec to consider new call
% bandwidth_call_freq=[500 5000]; % in Hz, ex for bearded seals vocalisations
% thresh_sigma=3

% appli_fit_type='yes'
% fit_type='poly3' % poly3 enables for oscillations in the tonal
% thresh_Rsquare=0.5 % at least 50% of variance explained by the fit

% %% run function:
% [sig_num, fs] = audioread('file.wav');
% % best to convert the signal in dB:
% sig_acou = 10^(-hydro_sensitivity/20).*10^(-gain/20)*dynamic*sig_num;
% [tonals] = func_tonal_detector(sig_acou,fs,dt_aPSD,bandwidth_call_freq,window_call_durations,space_btw_diff_call,thresh_sigma,appli_fit_type,fit_type,thresh_Rsquare)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off');
min_call_freq=bandwidth_call_freq(1);
max_call_freq=bandwidth_call_freq(2);
min_call_duration=window_call_durations(1);
max_call_duration=window_call_durations(2);
appli_fit_type=categorical(string(appli_fit_type));
hanning_window = 1024;
[S,Fq,Tps,P] = spectrogram(sig_acou,hann(hanning_window),0.5*hanning_window,hanning_window,fs);
hanning_window_2 = 128;  %  lower frequency resolution but higer temporal resolution
[S2,Fq2,Tps2,P2] = spectrogram(sig_acou,hann(hanning_window_2),0.5*hanning_window_2,hanning_window_2,fs);



time_to_test=1:(dt_aPSD*fs):((length(sig_acou))-(dt_aPSD*fs));
time_aPSD=(1:length(time_to_test))*dt_aPSD;

for nter=1:length(time_to_test)
    tdeb=floor(time_to_test(nter));
    tfin=ceil(tdeb+(dt_aPSD*fs));
    s_analyse2=sig_acou(((tdeb)):(tfin),1);
    N=length(s_analyse2);
    nfft=N;
    window=ones(size(s_analyse2));
    [psd_period,f_psd]=periodogram((s_analyse2-mean(s_analyse2)),window,nfft, fs,'psd'); % power spectral density
    psd_period2= psd_period(and(f_psd>min_call_freq,f_psd<max_call_freq));
    aPSD(nter)=mean(10*log10(psd_period2));

end



% Sliding window of a specified size (in samples) to calculate the standard
% deviation of the 5 values in the vector of averaged DSPs (aPSD).
% The window moves forward one sample at a time
%  -> this helps to better highlight the energy peaks
%
% Std = vector with standard deviations avec les écarts-types
% Std(i) = standard deviation of 5 values from aPSD : i-2, i-1, i, i+1, i+2

window_size = 5;  % must be an integer
for nter=1:length(time_to_test)
    if nter <= floor(window_size/2)
        tdeb=floor(time_to_test(1));
        tfin=floor(time_to_test(window_size));
        s_analyse2=sig_acou(((tdeb)):(tfin),1);
        N=length(s_analyse2);
        nfft=N;
        window=ones(size(s_analyse2));
        [psd_period,f_psd]=periodogram((s_analyse2-mean(s_analyse2)),window,nfft, fs,'psd'); % power spectral density
        psd_period2= psd_period(and(f_psd>min_call_freq,f_psd<max_call_freq));
        psd_period2 = 10*log10(psd_period2);
        Std(nter) = std(psd_period2);

    elseif nter >= (length(time_to_test)-round(window_size/2))
        tdeb=floor(time_to_test(length(time_to_test)-window_size));
        tfin=floor(time_to_test(length(time_to_test)));
        s_analyse2=sig_acou(((tdeb)):(tfin),1);
        N=length(s_analyse2);
        nfft=N;
        window=ones(size(s_analyse2));
        [psd_period,f_psd]=periodogram((s_analyse2-mean(s_analyse2)),window,nfft, fs,'psd'); % power spectral density
        psd_period2= psd_period(and(f_psd>min_call_freq,f_psd<max_call_freq));
        psd_period2 = 10*log10(psd_period2);
        Std(nter) = std(psd_period2);
    else
        tdeb=floor(time_to_test(nter-floor(window_size/2)));
        tfin=floor(time_to_test(nter+floor(window_size/2)+rem(window_size,2)));
        s_analyse2=sig_acou(((tdeb)):(tfin),1);
        N=length(s_analyse2);
        nfft=N;
        window=ones(size(s_analyse2));
        [psd_period,f_psd]=periodogram((s_analyse2-mean(s_analyse2)),window,nfft, fs,'psd'); % power spectral density
        psd_period2= psd_period(and(f_psd>min_call_freq,f_psd<max_call_freq));
        psd_period2= 10*log10(psd_period2);
        Std(nter)=std(psd_period2);
    end
end

threshold_detect = mode(round(Std,1))+thresh_sigma*std(Std);

clear pk_apsd pk_apsd2 tpeak tpeak2
Idx=find(Std>threshold_detect);
if length(Idx)>2
    pk_apsd(:,1)=Idx;% index where it happens

    % Count number of peaks
    k=1;
    pk_apsd(1,2)=k;
    for i=2:length(pk_apsd)
        if Idx(i)-Idx(i-1)<=(space_btw_diff_call/dt_aPSD) ;%% if there is a peak at (i-1) and also (i) thus still the same event
            pk_apsd(i,2)=k;
        else
            k=k+1; % otherwise new event => Ipca_raw(i)-Ipca_raw(i-1)>1
            pk_apsd(i,2)=k;
        end
    end

    events= unique(pk_apsd(:,2));
    k=0;
    % add duration condition on peak of aPSD
    for n=1:length(events)
        F=find(pk_apsd(:,2)==events(n));                             % position for event n
        tpeak(n)=time_aPSD(pk_apsd(F(round(end/2)),1));              % center time for the event
        duration_event=pk_apsd(F(end),1)-pk_apsd(F(1),1);           % samples
        if duration_event>=1/dt_aPSD*min_call_duration&duration_event<=1/dt_aPSD*max_call_duration ;
            k=k+1;
            F2=find(pk_apsd(:,2)==events(n));
            pk_apsd2(k,:)=pk_apsd(F2(round(end/2)),1:2);
            tpeak2(k,1)=pk_apsd(F2(1),1)*dt_aPSD;
            tpeak2(k,2)=pk_apsd(F2(end),1)*dt_aPSD;

        else
            k=k+1;
            pk_apsd2(k,1:2)=0;
            tpeak2(k,1)=0;
            tpeak2(k,2)=0;
        end
    end


    pk_apsd2=pk_apsd2(pk_apsd2(:,1)>0,:);
    tpeak2=tpeak2(tpeak2(:,1)>0,:);

    energetic_detector=horzcat(time_aPSD(pk_apsd(:,1))',pk_apsd(:,2)); % detector of the energetic peaks
    energetic_duration_detector=horzcat(tpeak2,pk_apsd2(:,2)); % detector of the energetic peaks with the time condition
    %% Contouring + fit model
    c=0;
    raw_tonals={};
    smooth_tonals={};
    for k=1:length(pk_apsd2(:,1))
        if appli_fit_type == 'no'
            Tcall=Tps2(Tps2>=tpeak2(k,1)&Tps2<=tpeak2(k,2));
            Pcall=P2(:,(Tps2>=tpeak2(k,1)&Tps2<=tpeak2(k,2)));
        else
            Tcall=Tps(Tps>=tpeak2(k,1)&Tps<=tpeak2(k,2));
            Pcall=P(:,(Tps>=tpeak2(k,1)&Tps<=tpeak2(k,2)));
        end
        Pcall=10*log10(abs(Pcall));
        clear A
        for jj=1:length(Tcall)
            A(jj) = (find(Pcall(:,jj)==max(Pcall(min(find(Fq>=min_call_freq)):end,jj)))); % above the minimum frequncy band

        end

        if appli_fit_type == 'no'
            c=c+1;
            raw_tonals{c}=[Tcall' Fq(A)];
        else
            [f,gof] = fit(Tcall',Fq(A),fit_type);
            if(round(gof.rsquare,1)>=thresh_Rsquare)
                c=c+1;
                xq=min(Tcall):0.2:max(Tcall); % smoothing the detection to 0.2s
                vq = interp1(Tcall,Fq(A),xq,'pchip');
                raw_tonals{c}=[Tcall' Fq(A)];
                smooth_tonals{c}=[xq' vq'];

            end
        end

    end

else
    c=0;
    raw_tonals={};
end