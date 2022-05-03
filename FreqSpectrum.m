function [IBI_FFT, IBI_Fp, tDurationAnalyzed, Nwindow_out, Noverlap_out] = FreqSpectrum(t_finalshift, IBI_final, Nwindow ,Noverlap , Nfft ,Fs_interp, Detrend_vec, ARorder, Filter, FreqMethod); 
% [IBI_FFT, IBI_Fp, tDurationAnalyzed, Nwindow_out, Noverlap_out] = FreqSpectrum(t_finalshift, IBI_final, Nwindow ,Noverlap , Nfft ,Fs_interp, Detrend_vec, ARorder, Filter, FreqMethod); 
% Compute the power spectral density estimate of IBI_Final using various 
% methods (Welch, AR, simple FFT) with various options for windows filtering 
% (Filter), overlap (Overlap), padding increase resolution (NFFT) and
% detrending the signal (polynomial or smoothing; Detrend_vec)
% POwer spectral density method is given by FreqMethod: 
%       1. (used to be 1, 2 or 3). Classic pwelch using overlap, padding (NFFT)
%       2 (was 4). Lomb-Scargle No interp
%       3 (was 5). Simple Fourier, no windows (rect window)/overlap
%       4 (was 7). Simple fourier with windows filtering of size Nt  
%       5 (was 6). Wavelet (not actually doing anything as it is implemented in SignalBrowser)
%       6 (was 8). Autoregressive using Burg method 
%       7 (was 9). Autoregressive using Yule-Walker method  
%       8 (was 10). Autoregressive using Covariance method 
% 
% Window filtering option given by Filter
%       1. None
%       2. Hanning
%       3. Hamming
%       4. Bartlett
%       5. Blackman
%
% Didier: 4/26/2022 (with SignalBrowser version 4.6)
% Reduce, cleaned up and fixed the number of psd choices

FreqMethodString = {'pwelch','Lomb-Scargle','FFT with no window filtering','FFT with window filtering',...
                    'Morlet', 'AR - Burg','AR - Yule-Walker', 'AR - Covariance' };
Nwindow_out = Nwindow;
Noverlap_out = Noverlap;
Nt = length(t_finalshift);
Nseg = 1;
Norder = Detrend_vec(1); % Order of the polynomial detrending <3
smoothN = Detrend_vec(2); % otherwise uses smoothing function

fprintf('psd method %d (%s): \n', FreqMethod,FreqMethodString{FreqMethod});


if FreqMethod == 2 %Lomb-Scargle No interp
    if Norder <3
        [P, S ] = polyfit(t_finalshift, IBI_final, Norder);
        IBI_detrend = IBI_final - polyval(P,t_finalshift);
    else
        IBI_detrend = IBI_final - smooth(IBI_final,smoothN);
    end
    
    if isempty(Fs_interp) 
        Fs_interp = 1/(t_finalshift(2) - t_finalshift(1));
    end
    %     IBI_Fp = (0:Nfft)'*Fs_interp/Nfft;
    %     [IBI_FFT,Prob] = lomb(t_finalshift,IBI_detrend,IBI_Fp);
    [IBI_FFT,IBI_Fp] = plomb(IBI_detrend,t_finalshift);
    Nwindow_out = Inf;
    Noverlap_out = 0;
    tDurationAnalyzed = t_finalshift(end)-t_finalshift(1);
else 
    %interpolating first if needed as signal must be evenly sampled 
    if ~isempty(Fs_interp) % = [] ;%1/Fs_interp ~= t_finalshift(2) - t_finalshift(1) 
        % interpolating 
        tR_interp = (t_finalshift(1):1/Fs_interp:t_finalshift(Nt))';
        IBI_interp = interp1(t_finalshift, IBI_final,tR_interp,'spline');
    else %no need for interpolation
       tR_interp = t_finalshift;
       IBI_interp = IBI_final;
       Fs_interp = 1/(t_finalshift(2) - t_finalshift(1));
   end
    %rg = find(tR_interp>=tIBI_final(
    % Detrending
    %         [P, S ] = polyfit(tR_interp(rg), tIBI2(rg),Norder);
    %         tIBI_detrend = tIBI2(rg) - polyval(P,tR_interp(rg));
    
    if Norder <3
        [P, S ] = polyfit(tR_interp, IBI_interp, Norder);
        IBI_detrend = IBI_interp - polyval(P,tR_interp);
    else
        IBI_detrend = IBI_interp - smooth(IBI_interp,smoothN);
    end
   
   
   if Nwindow == 0 %(default pwelch values of 8 with 50% overlaping windows)
        Nwindow = [];
        Noverlap = [];
   elseif Nwindow == inf %full signal length
       Nwindow = length(tR_interp);
       Nwindow_out = Inf;
   %elseif length(tR_interp)<(Nwindow*2 - Noverlap) || Nwindow == Inf %#ok<OR2>
   elseif length(tR_interp)<(Nwindow) || Nwindow == Inf %#ok<OR2>
        %if FreqMethod == 7
            Nwindow = length(tR_interp);%length(rg);
            Nwindow_out = Inf;
        %end
        Noverlap = 0;
        Noverlap_out = 0;
    end
    Dt = max(tR_interp) - tR_interp(1);
    
    fprintf(' Dt = %0.1f; N = %d;  Win = %d; Over = %d; \n',Dt, length(tR_interp), Nwindow, Noverlap);
    
    if FreqMethod ~= 1 
       %force  Nwindow to be length of tR_interp
       % pwelch is only used for FreqMethod = 1 
       Nwindow = length(tR_interp);
       Nwindow_out = Inf;
       Noverlap = 0;
       Noverlap_out = 0;
    end
    
    if Filter == 2; win = hamming(Nwindow);
    elseif (Filter == 3); win = hanning(Nwindow);
    elseif (Filter == 4); win = bartlett(Nwindow);
    elseif (Filter == 5); win = blackman(Nwindow); 
        %elseif (Filter == 6) xn = (1:Nwindow)'; x=(xn-floor(Nwindow/2))/(Nwindow-1)*2*pi; win =(1+cos(x))/2;
    else win = ones(Nwindow,1); %rectangular i.e. no windows
    end;
    
    % Duration of Signal Analyzed
    if FreqMethod <= 3
        % Compute the number of segments
        Nseg = (length(IBI_detrend)-Noverlap)./(Nwindow-Noverlap);
        %%%%%%%%%%%%%%%%%% Didier  - SignalBrowser V4.2 4/27/2010 %%%%%%%%%%%%%
        %tDurationAnalyzed = Dt - (Nseg - fix(Nseg))*Fs_interp;
        tDurationAnalyzed = ((Nwindow-Noverlap)*fix(Nseg) + Noverlap-1)/Fs_interp;
        %%%%%%%%%%%%%%%%%% Didier 4/27/2010 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        tDurationAnalyzed = Dt;
    end
    fprintf('; Nwin = %d, DtAnalyzed = %0.1f \n',fix(Nseg), tDurationAnalyzed)
    
    %%%%%%%%%%%%%%%%%% Didier  - SignalBrowser V4.6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Removed Choice 2 and 3; Redundant to method 1 after correcting minor
    % issues which were differentiating them 
    % revert to computing psd as "onesided" to include all power

    if FreqMethod == 1 % Welch 1; Uses most recent Matlab pwelch function
        [IBI_FFT,IBI_Fp] = pwelch(IBI_detrend,win,Noverlap,Nfft,Fs_interp);
        IBI_FFT = IBI_FFT;%/2;
%   Keep the following commented line as information for prior
%   SignalBrowser versions <4.6
%     elseif  FreqMethod == 2 % Welch 2; uses a very old pwelch function??? 
%         % Since Filter is used in computing win, I am not sure if Filter as
%         % input is used (??)
%         [IBI_FFT,IBI_Fp] = pwelch2(IBI_detrend,win,Noverlap,Nfft,Fs_interp,Filter);
%         IBI_FFT = IBI_FFT/Fs_interp;
%     elseif FreqMethod == 3 % Welch 3; used an old pwelch function?? 
%         % Since Filter is used in computing win, I am not sure if Filter as
%         % input is used (??)
%         %[IBI_FFT1,IBI_Fp1] = pwelch3(IBI_detrend/1000,Nwindow,Noverlap,Nfft,Fs_interp,[],3);
%         [IBI_FFT,IBI_Fp] = pwelch3(IBI_detrend,win,Noverlap,Nfft,Fs_interp,'twosided',Filter);
%         IBI_FFT = IBI_FFT;%/2; 
    elseif FreqMethod == 3 || FreqMethod == 4 || FreqMethod == 5 
                            % FFT window; Uses most recent Matlab pWelch function 
                            %  with no overlaping, Nwindow = Nt
                            % For FreqMethod = 5, for IBI it compute the
                            % wavelet psd within SignalBrowser (FreqMethod
                            % is not called). However, for any other signal
                            % (respiration), FreqSpectrum is called and will 
                            % compute the simple FFT psd with window filtering as
                            % wavelet decomposition is not computed 
        Noverlap = 0;
        if FreqMethod == 4 % no window filtering 
            win = ones(Nwindow,1); % rect window
        end
        % using pwelch (with rect win) of size Nt is similar to fft psd
        [IBI_FFT,IBI_Fp] = pwelch(IBI_detrend,win,Noverlap,Nfft,Fs_interp);
        
    elseif FreqMethod == 5 % Wavelet 
        % do nothing as it is computed by SignalBrowser

    elseif FreqMethod == 6 % AR method - Burg  
        [IBI_FFT,IBI_Fp] = pburg(IBI_detrend.*win,ARorder,Nfft,Fs_interp);
        % Using spectrum function is no longer recommended   
        %         Hburg = spectrum.burg(ARorder);
        %         hpsd = psd(Hburg,IBI_detrend,'Fs',Fs_interp,'NFFT',Nfft,'SpectrumType','twosided');
        %         IBI_FFT = hpsd.Data;
        %         IBI_Fp = hpsd.Frequencies;
    elseif FreqMethod == 7 % AR method - Yule-Walker  
        [IBI_FFT,IBI_Fp] = pyulear(IBI_detrend.*win,ARorder,Nfft,Fs_interp);
        % Using spectrum function is no longer recommended   
    elseif FreqMethod == 8 % AR method - Cov  
        [IBI_FFT,IBI_Fp] = pcov(IBI_detrend.*win,ARorder,Nfft,Fs_interp);
        % Using Spectrum function is no longer recommended   

    end
    
    %%%%%%%%%%%%%%%%%% Didier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
