function [N_SCR, nSCR] = Compute_SCR(t,y)
% function [N_SCR, nSCR] = Compute_SCR(t,y) computes the Skin Conductance
% responses or so called electrodermal response) from the continuous skin 
% conductance vector y and associated time samples.
% Returns the number of responses N_SCR and sample indices of y of those
% responses nSCR
% Local peak responses must be greater than EDRamp=0.5us and at least
% spaced 1s apart


% making sure that y is a vertical vector
[nc,nr] = size(y);
if nr>nc
    y=y'; t=t';
end

%---------- Didier v3.3 --------------
N_SCR = 0;
Trisetime_min = 1;  %1 sec latency
Trisetime_max = 10; %10 sec max latency
EDRamp = 0.05;      % 0.5uS EDR ampl
Ns = length(y);
% find SCL/EDR peaks
npeak = find( (y(2:Ns-2)>y(1:Ns-3)) & ...
    ( (y(2:Ns-2)>y(3:Ns-1)) | ...
    ( (y(2:Ns-2)== y(3:Ns-1)) & (y(3:Ns-1)>y(4:Ns))  ) ...
    ) )+1;

nSCR = [];
for ns = npeak'
    ns1 = ns;
    ns2 = ns;
    % go backward to find preceding index ns1 of the lowest y value prior to
    % peak ns within the past Trisetime_max=10s 
    % this is to find the local amplitude of the EDR peak (y(ns2)-y(ns1))
    while y(ns1)>=y(ns1-1) & ns1>2 & (t(ns2) - t(ns1))<Trisetime_max+1
        ns1 = ns1 - 1;
    end
    % Check that the local peak amplitude is greater than EDRamp and no
    % closer than Trisetime_min/2 to consider it an electrodermal response 
    if (y(ns2) - y(ns1)) >= EDRamp & (t(ns2) - t(ns1))>Trisetime_min/2 & (t(ns2) - t(ns1))<=Trisetime_max
        N_SCR = N_SCR + 1;
        nSCR = [nSCR;ns];
    end
end   
% fig = figure(gcf);
% figure(5); plot(t,y,'b',t(npeak),y(npeak),'.b',t(nSCR),y(nSCR),'or')
% figure(fig)