function [RRpos, RRout,tpos] = findRpeaks(ECGsig, tsig, Fs, RRmin, MODEin, RRpos_old) 
% function [RRpos, RRout,tpos] = findRpeaks(ECGsig, tsig, Fs, RRmin, MODEin, RRpos_old)
% 
% Compute Rpeaks in the ECG signal given by ECGsig and corresponding time
% instance tsig, sampling frequency Fs. 
% returns indices of Rpeaks in RRpos, corresponding time in tpos
% RRout is no longer computed; returns empty vector
% RRmin defines the smallest RR intervals
% MODEin select the method to be used:
%   MODE = 1 Hilbert transform
%   MODE = 2 Regular threshold detection
%   MODE = 3 Filter and use envelope of the e.c.g to delineate the QRS complex 
%               (see  processing2.m function from the BIOSIG toolbox)
%   MODE = 4 uses prior Rpeaks indices in RRpos_old to recompute slightly 
%           displaced prior R peaks(due to filtering , etc...) within +/- 6
%           data samples 

% Afilt = [0.028  0.053 0.071  0.053 0.028];
% Bfilt = [1.000 -2.026 2.148 -1.159 0.279];

h = gca;
m = length(ECGsig);
fprintf('t = %5.1f \n',m/Fs)

RRout = [];
tpos = [];
RRpos = [];
kmin = floor(RRmin*Fs);
RRpos1 = -2*kmin;

MODE = 3; % default mode

if nargin >= 5
    MODE = MODEin;
end
fprintf('Finding Rpeaks MODE %d\n',MODE)

if MODE == 4 & nargin < 6
    fprintf('You must specify old RR position for MODE 4\n')
    return
end


if MODE == 1 | MODE == 2 | MODE == 3
    if MODE == 3 %using processing2 BIOSIG toolbox
        %y_filt = processing2({'ECG_envelope',Fs},(ECGsig>0).*ECGsig);
        y_filt = processing2({'ECG_envelope',Fs},ECGsig);
        m_y = mean(y_filt(find(~isnan(y_filt))));
        rmsv = norm(y_filt-m_y)/sqrt(length(y_filt));
        y_filt = y_filt - m_y+rmsv/2; %"keep y_filt just above y =0, remove DC component"
    end
    
    i1 = 1; i2 = 1024*2; 
     
    while i1 + i2 + 1 <m
        rg1 = i1+(1:i2)';
        if MODE == 1; % use Hilbert transform to detects R peaks
            y = (ECGsig(rg1+1)-ECGsig(rg1-1))/2*Fs;
            yfft = fft(y);
            yfft(1) = 0;
            yhilb = real(ifft([-j*yfft(1:i2/2) ; j*yfft(i2/2+1:i2)]));
            y_RR = yhilb;
        elseif MODE == 2; % use ECG sig to detect R peaks
            y_RR = ECGsig(rg1);
        else MODE == 3; % use filtered signal
            y_RR = y_filt(rg1);
        end
%             figure(2);clf
%             plot(tsig(rg1+1),y,'k',tsig(rg1+1),yhilb,'r');
%             figure(1); clf
%             plot(tsig(rg1+1),ECGsig(rg1+1),'b',tsig(rg1+1),yhilb/500,'r')
% 
%             figure(2);clf
%             plot(tsig(rg1+1),ECGsig(rg1+1),'b',tsig(rg1+1),y_RR,'r');
%             pause
        
        % calculate threshold
        rmsv = norm(y_RR)/sqrt(length(y_RR));
%         fprintf('r = %2.2g -\n',100*rmsv/max(y_RR));%
        if rmsv < 0.45*max(y_RR)
            if rmsv > 0.1*max(y_RR)
                th = 0.4*max(y_RR);
            else
                th = 0.2* max(y_RR); %2*rmsv;
            end
            
            k1 = min(find(y_RR>th));
            rgk = k1:min(i2,k1+kmin);
            kmax = k1-1+min(find(y_RR(rgk)==max(y_RR(rgk))));
            %fprintf('kmax=%2.3g \n',kmax);
            
            if isempty(RRpos); RRmean = max(ECGsig(rg1)); 
            else
                RRmean = mean(ECGsig(RRpos)) ;
            end
            RRpos0 = i1 -1 + kmax;
            y_loc =  ECGsig( max(1,RRpos0-kmin):min(m,RRpos0+kmin) );
            rms_loc = norm(y_loc.*(y_loc>mean(y_loc)))/sqrt(length(y_loc)); %calculate only for y_loc>mean(y_loc) since
                                                                            % some ECG have bad QRS complex [symetrique]
%             fprintf('i1 = %d, %d : %d - %d - %d (rmsloc = %2.1g-%2.1g//max=%2.1g) \n',i1, kmax,kmax<min(i2,k1+kmin) & kmax > 1, ... 
%                     RRpos0 - RRpos1 > kmin ... 
%                     , max(y_loc)>3*rms_loc,rms_loc,norm(y_loc)/sqrt(length(y_loc)) , max(y_loc))   
%                     figure(3); plot(y_RR)

            if kmax<min(i2,k1+kmin) & kmax > 1 ... %make sure it is a peak and not a boundary maxima
                    & RRpos0 - RRpos1 > kmin ... %accept only RR larger than RRmin
                    & (max(y_loc)>3*rms_loc) %exclude peaks in very locally noisy signal
                     %& ECGsig(RRpos0) > min(0.3*ECGsig(max(RRpos1,1)),RRmean)... % exclude peaks too high (noise/artefact)
%                 fprintf('OK1 \n')           
                RRpos1 = RRpos0;
                %find the position of the maximum ECGsig starting from 
                %estimated location using Hilbert transform
                
                nw = 5; % width of where to look for maximum
                RRpos2 = RRpos1;
                
                if MODE == 1
                    nw1 = 5;
                    rgf = max(1,RRpos1-nw1):min(m,RRpos1+nw1);
                    RRpos1 = rgf(1) -1 + min( find(ECGsig(rgf) == max(ECGsig(rgf))) );
                elseif MODE == 3
                    nw1 = floor(RRmin*Fs);
                    rgf = max(1,RRpos1-nw1):min(m,RRpos1+1);
                    if rgf(1)==1; rgf = 1:10; end
                    if rgf(length(rgf)) == m; rgf = m-10:m; end
                    RRpos1 = rgf(1) -1 + min(find(ECGsig(rgf) == max(ECGsig(rgf))));
                end
                
                if ECGsig(RRpos1) < min(0.3*ECGsig(max(RRpos2,1)),RRmean) % exclude peaks too high (noise/artefact)
                    RRpos1 = RRpos2;
                    i1 = i1 + kmax + floor(kmin/2);
%                     fprintf('jump')
                    continue; %jump to next iteration of while loop
                end
                % Find a more precise value in time using interpolation
                rgf = max(1,RRpos1-nw):min(m,RRpos1+nw);
                if rgf(1)==1; rgf = 1:10; end
                if rgf(length(rgf)) == m; rgf = m-10:m; end

                ECGsig_interp = interp(ECGsig(rgf),10);
                npos = min(find(ECGsig_interp == max(ECGsig_interp)));
                tpos2 = tsig(RRpos1) + (npos-nw*10-1)/10/Fs; 
                %---------PLOT -------------
%                 figure(2); plot(tsig(RRpos1)+(-nw*10:(nw*10+9))'/10/Fs, ECGsig_interp,'.b',...
%                                 tsig(RRpos1+(-100:100)'),ECGsig(RRpos1+(-100:100)'),'b',...
%                                 tpos2,ECGsig_interp(npos),'+b',...
%                                 tsig(RRpos1), ECGsig(RRpos1),'+r', ...
%                                 tsig(RRpos2), ECGsig(RRpos2),'or',...
%                                 tsig(RRpos2 - kmax + (max(kmax-140,1):min(kmax+40,i2)) ),y_RR( max(kmax-140,1):min(kmax+40,i2)) /(1000*(MODE==1)+0.2*(MODE == 3)),'r')
%                 pause          
%                 pause(0.01)
                tpos = [tpos;tpos2];
                RRpos = [RRpos; RRpos1];
            end
            i1 = i1 + kmax + floor(kmin/2); 
        else
            %i1 = i1 + kmax + 1;
            i1 = i1 + kmin + 1;
        end
        %pause
        if rem(length(RRpos),5)==0
            fprintf('%d (%5.1f) ',length(RRpos),RRpos1/Fs);
        end
        if rem(length(RRpos),1000)==0
            fprintf('\n');
        end
        %pause
        %       
        %     
        %     figure(2);clf
        %     plot(tsig(rg1+1),y,'k',tsig(rg1+1),y_RR,'r');
%             figure(1); clf
%             rgp = (max(i1-4000,1):min(i1+2000,m))';
%             rgp2 = find(RRpos>i1-4000 & RRpos<i1+2000);
%             %plot(tsig(rgp),ECGsig(rgp),'k',tsig(rg1),y_RR/(1000*(MODE==1)+0.2*(MODE == 3)),'b',tsig(RRpos(rgp2)),ECGsig(RRpos(rgp2)),'.r');
%             plot(tsig(rgp),ECGsig(rgp),'k',tsig(rg1),y_RR/(1000*(MODE==1)+(MODE == 3)),'b',tsig(RRpos(rgp2)),ECGsig(RRpos(rgp2)),'.r');
%            pause(0.01)
    end
    

%     figure(1); clf
%     plot(tsig,ECGsig,'k',tsig(RRpos),ECGsig(RRpos),'.r',tsig(RRout),ECGsig(RRout),'bo');
    
end
    
if MODE == 4 
    mR = length(RRpos_old);
    RRpos = [];
    tpos = [];
    for i = 1:mR
        rg = max(1,min(RRpos_old(i)-6,m)):min(max(RRpos_old(i)+6,1),m);
        RRpos1 = rg(1) - 1 + find(ECGsig(rg) == max(ECGsig(rg)));
        if RRpos1 > rg(1) & RRpos1 < rg(length(rg))
            RRpos = [RRpos; RRpos1];
            % Find a more precise value in time using interpolation
            nw = 5;
            rgf = max(1,RRpos1-nw):min(m,RRpos1+nw);
            if rgf(1)==1; rgf = 1:10; end
            if rgf(length(rgf)) == m; rgf = m-10:m; end
            
            ECGsig_interp = interp(ECGsig(rgf),10);
            npos = min(find(ECGsig_interp == max(ECGsig_interp)));
            tpos2 = tsig(RRpos1) + (npos-nw*10-1)/10/Fs; 
            tpos = [tpos;tpos2];
        end
    end
    
end

if ~isempty(RRpos)
    %find outsider
    %case of noise
%     nout = find( ECGsig(RRpos) > 3*RRmean);
%     mR = length(RRpos)
%     IBI = ECGsig(2:mR)-ECGsig(1:mR-1);
%     %distIBI = sqrt((IBI(2:mR-2)-IBI(1:mR-3)).^2 + (IBI(2:mR-2)-IBI(3:mR-1)).^2)/2;
%     nIBI = find((IBI(2:mR-2) >= 1.2*IBI(1:mR-1) & IBI(2:mR-2) >= 1.2*IBI(1:mR+1)) | ...
%         (IBI(2:mR-2) <= 0.8*IBI(1:mR-1) & IBI(2:mR-2) <= 0.8*IBI(1:mR+1)) );
%     nout = [nout; RRpos(nIBI+2)];
%     RRout = RRpos(nout);
%     size(RRpos)
end

axes(gca)    
    
    
    