function [tEDR, EDRsig] = EDR(tsig, ECGsig, nR, dn_QR, dn_RS, Fs1);

% function [tEDR, EDRsig] = EDR(tsig, ECGsig, nR, dn_QR, dn_RS, Fs1);
% compute the ECG derived respiration signal based on amplitude change in 
% QRS with chest wall movement. 
% EDR takes the ECG and corresponding time signal ECGsig and tsig, 
% indices of Rpeaks in nR and sampling frequency Fs1. 
% dn_QR and dn_SR provides an estimate of the distance from R to the Q and S 
% complex (typically dn_QR = 40 and dn_QS = 30 samples at 1000Hz) to compute  
% the area under the QRS curve (METHOD=1) or remove drift (METHOD=2) 




% METHOD = 1;  %METHOD = 1: use area under curve...
% METHOD = 2;  %METHOD = 2: use ECG ampl...
METHOD = 2;  %METHOD = 2: use ECG ampl... (this methods works better)

m = length(tsig(nR));
m2 = length(tsig);
tR = tsig(nR);
tIBI = 1000*(tR(2:m)-tR(1:m-1));

%Rsig = ECGsig(nR);
% %make sure that the extreme QRS complexes are whole in
N_QS = length(dn_QR); %N_QS = 1 for single value
Rsig1 = zeros(m,1);
Rsig2 = zeros(m,1);
for i = 1:m %remove extremities
    dn_QR = dn_QR(min(i,N_QS)); % set dn_QR = 
    % define index
    n1 = max(1,nR(i) - round(dn_QR*Fs1/1000));
    n2 = min(m2,nR(i) + round(dn_RS*Fs1/1000));
    p = polyfit([tsig(n1)   tsig(n2)], [ECGsig(n1)   ECGsig(n2)] ,1);

    y_baseline = polyval(p,tsig(n1:n2));
    y = ECGsig(n1:n2) - y_baseline;
    %Rsig(i-1) = (tsig(2)-tsig(1))*( 0.5*y(1) + 0.5*y(n2-n1+1) + sum(y(2:n2-n1)) ) ; 
    Rsig1(i) = (tsig(2)-tsig(1))*trapz(y(1:(n2-n1+1))); % area under the curve
    Rsig2(i) = ECGsig(nR(i))-(ECGsig(n1) +  ECGsig(n2))/2; %take amplitude of ECG and remove drift 
end
%Rsigplot = mean(tIBI) + 2*(Rsig-mean(Rsig))/std(Rsig)*std(tIBI);
Rsigplot = 2*(Rsig1-mean(Rsig1))/std(Rsig1);
%Rsigplot = 2*(Rsig-mean(Rsig1))/std(Rsig1);

%Rsig2 = ECGsig(nR);
Rsig2 = Rsig2 - mean(Rsig2); Rsig2 = Rsig2/std(Rsig2);

% figure(2); cla
% hold on
% length(tR)
% length(Rsigplot)
% tIBIm = mean(tIBI);
% uda.hRsig = plot(tR(1:m-1), (tIBI-tIBIm)/std(tIBI-tIBIm),'b',...%tsig, ECGsig-mean(ECGsig),'k',
%      tR(3:m-1),Rsigplot,'-r',...
%     tR(3:m-1),Rsigplot,'.r',tR(3:m-1),Rsig2,'m');
% axis auto
% %axis([tsig(1) tsig(m2) 0.9*min(tIBI) 1.1*max(tIBI)])
% hold off;

if METHOD == 1  %use area under curve...
    EDRsig = Rsig1;
    %tEDR = (tR(2):1/Fs1:tR(m-1))';
else %use amplitude of ECG
    EDRsig = Rsig2;
    %tEDR = (tR(2):1/Fs1:tR(m-1))';
end
tEDR = (tsig(1):1/Fs1:max(tR(m)))';;
%---- Extract EDR signal
% resample at 4 Hz
Fs1 = 4;


EDRsig1 = (EDRsig - mean(EDRsig))/std(EDRsig); tEDR1 = tR(1:m); % unused; used to check processing

% Interpolate before filtering
EDRsig = interp1(tR(1:m), EDRsig,tEDR,'spline');

EDRsig2 = (EDRsig - mean(EDRsig))/std(EDRsig); tEDR2 = tEDR; % unused; used to check processing

% build the band-pass filter 1/30Hz - 1Hz
 fp1 = 0.05; fp2 = 0.5; 
 [b,a]=butter(4,2*[fp1 fp2]/Fs1);
% filter
EDRsig=filter(b,a,EDRsig);
  
% normalize and center again
EDRsig = (EDRsig - mean(EDRsig))/std(EDRsig);


% figure(4); cla
% 
% tEDR3 = tEDR;
% EDRsig3 = smooth(EDRsig,10); EDRsig3 = EDRsig3-smooth(EDRsig,200);EDRsig3 = EDRsig3/std(EDRsig3);
% 
% 
% plot(tEDR,EDRsig,'-b',tEDR1,EDRsig1,'-k',tEDR2,EDRsig2,'-r',tEDR3, EDRsig3,'m',tR(1:m-1), (tIBI-tIBIm)/std(tIBI-tIBIm),'k--');
% axis auto
% hold off;



%EDRtest


