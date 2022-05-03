function [f_mean, f_median, f_std] = geometric_measures(f,pwr,METHOD)
% function [f_mean, f_median, f_std] = geometric_measures(f,pwr,METHOD)
% 
% Compute and return frequency related geometric measures (mean, median, std)
% of the "pwr" distribution. pwr is a vector containing the power
% spectrum values at frequency f
% 

% This is to validate the computation
% f = [1 2 3 4 5 6 7]';
% pwr = [0.1 0.3 0.4 0.2 0.4 0.05 0.01]';
% f = 0:0.01:10;
% sigma = 0.1; fc = 6.3;
% pwr = exp(-(sqrt(f)-sqrt(fc)).^2/sigma); %+exp(-(f-3).^2/2);

% plot(f,pwr);
% pause

if nargin <3
    METHOD =1;
end
area_all = abs(trapz(f,pwr));
pwr_sum = sum(pwr);

[n, nx] = size(pwr);

t1 = cputime;

if METHOD == 1
    % calculate mean median, std of freq distribution at 1% accuracy using
    % integration method

    f_mean = (f'*pwr)./pwr_sum;
    for j = 1:nx
        f_std(j) = sqrt( (f-f_mean(j))'.^2*pwr(:,j)./pwr_sum(j) );
    end

    % calculating median
    n = length(f);
    for j = 1:nx
        area = abs(trapz(f(1:round(n/2)), pwr(1:round(n/2),j) ));
        if area<area_all(j)/2

            for i = round(n/2)-1:length(f);
                area_old = area;
                area = abs(trapz(f(1:i),pwr(1:i,j)));
                if area> area_all(j)/2 ;
                    i_f = i; break;
                end;
            end
            f_median(j) = f(i_f-1)+(f(i_f)-f(i_f-1))*(area_all(j)/2-area_old)/(area - area_old);
            %         catch
            %             1==1;
            %         end
        else
            for i = round(n/2)+1:-1:1;
                area_old = area;
                try
                    area = abs(trapz(f(1:i),pwr(1:i,j)));
                catch
                    1==1;
                end
                if area< area_all(j)/2
                    i_f = i; break;
                end;
            end
            f_median(j) = f(i_f+1)+(f(i_f)-f(i_f+1))*(area_all(j)/2-area_old)/(area - area_old);
        end
    end
else
    % calculate mean median, std of freq distribution at 1% accuracy using
    % bin method
    Npoint = 100;
    pwr = FFTpwr(rg,nx)/min(min(FFTpwr(find(FFTpwr(rg,nx)>0))));
    pwr = pwr/max(max(pwr))*Npoint;

    for j = nx
        f_all = [];
        for i = 1:max(1,round(NFFT/Npoint/10)):NFFT
            f_all = [f_all; repmat(f(i),round(pwr(i,j)),1)];
        end
        f_median(j) = median(f_all);
        f_mean(j) = mean(f_all);
        f_std(j) = std(f_all);
    end
end





%[f_mean f_mediana f_std]
% cputime-t1
% plot(f,pwr,[f_mean f_mean],[0 max(pwr)],'r',[f_mediana f_mediana],[0 max(pwr)],'b',[f_mean-f_std/2 f_mean+f_std/2],[max(pwr)/2 max(pwr)/2] )
%
% t2 = cputime
% f_all = [];
% Npoint = 1000;
% pwrN = pwr/min(pwr(find(pwr>0)));
% pwrN = pwrN/max(pwrN)*Npoint;
% for i = 1:max(1,round(n/Npoint)):n
%     f_all = [f_all; repmat(f(i),round(pwrN(i)),1)];
% end
% [mean(f_all) median(f_all) std(f_all)]
%
% cputime - t2
%
% [SpectrumPwr, SpectrumMax, f_max, f_mean, f_median, f_std] = HRV_measures(f,[pwr;pwr],[],min(f),max(f),1,2)