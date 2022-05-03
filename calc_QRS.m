function [nQR, nRS] = calc_QRS(x,y,tR,METHOD) 
nR = min(find(x>tR))
tR
dn1 = 90; dn2 = 30; m=length(y);
n1 = max(nR-dn1,1); n2 = min(m,nR+dn2);
nR = find( y==max( y(n1:n2) )  );

%METHOD = 2; % choose method

% calculate derivative
dydx = zeros(m,1);
dydx = (y(2:m)-y(1:m-1))./(x(2:m)-x(1:m-1));size(dydx)
dydx(m)=dydx(m-1);
% 5 points derivative
% dydx(5:m) = 2*y(5:m) + y(4:m-1) - y(2:m-3) - 2* y(1:m-4);
% dydx(1:4)=dydx(5);
dydx(3:m-2) = 2*y(5:m) + y(4:m-1) - y(2:m-3) - 2* y(1:m-4);
dydx(1:2)=dydx(3);
dydx(m-1:m)=dydx(m-2);


switch METHOD
    case 1
        %Method 1
        % based on first derivative 
        % finding Q point first minimum for derivative before first derivative zero crossing prior to R point 
        d2yd2x = zeros(m,1);
        d2yd2x(3:m-2) = 2*dydx(5:m) + dydx(4:m-1) - dydx(2:m-3) - 2* dydx(1:m-4);
        %d2yd2x = dydx(2:m) - dydx(1:m-1);
        
        
        %%%% Defining QRS onset defined as nQ
        nQ=nR;
        nQ1 = n1-1 + max(find(dydx(n1:nR-5)<0)); %first minima for y 
        n3 = max(1, nQ1-100);
        nQa = n3-1 + max(find(dydx(n3:nQ1-5)>=0));   %first derivative crossing zero
        %nQ = nQa;
        % Q point is defined as the first derivative minima (inflexion
        % point for ECG) in case there is no zero crossing or it is way too
        % far
        nQ2 = n3-1 + max(find(d2yd2x(n3:nQ1-1)<0)); % use inflexion point (minimum for derivative)
        nQb = n3-1 + max(find(d2yd2x(n3:nQ2-1)>=0)); 
        nQ = max(nQa,nQb)
        if isempty(nQa); nQ = nQb; end
        if x(nR) - x(nQa) > 0.08
          %nQ = nQb;  
        end

        %%%%%%%%%%%% Defining QRS defined as nS
        % finding S point first minimum for derivative before first derivative zero crossing prior to R point 

        nS=nR;
        n4 = min(m,nR+100);
        nS1 = nR+4 + min(find(dydx(nR+5:n4)>0)); %first minima for y 
        n5 = min(m,nS1+80);
        nSa = nS1 + 4 + min(find(dydx(nS1+5:n5)<=0)); %first derivative crossing zero
        
        % S point is defined as the first derivative minima (inflexion
        % point for ECG) in case there is no zero crossing or it is way too
        % far
        nS2 = nS1+3 + min(find(d2yd2x(nS1+4:n5)<=0)); % find 
        nSb = nS2+3 + min(find(d2yd2x(nS2+4:n5)>=0)); 
        nS = min(nSa,nSb);
        if isempty(nSa); nS = nSb; end
        if x(nSa) - x(nR) > 0.08
          %nS = nSb;  
        end
        

%         figure(1); plot(x,y,'b',x,dydx/max(abs(dydx)),'r-',x,d2yd2x/max(abs(d2yd2x)),'r--',x(nQ),y(nQ),'xk',x(nS),y(nS),'xk')

    case 2
        
        %Method 2
        % based on 
        %"Pan & Tompkins, 1985, IEEE bme, A real-time qrs detection algorithm"
        % Moving Averaging Window Integration
        Wt = 0.180;% ms
        N = round(0.180/(x(2)-x(1))); LM = zeros(m,1);
        for i = N+5:m
            LM(i) = 1/N*sum(dydx(i-N+1:i).^2);
            %LM(i) = 1/N*sum((dydx(i-N+1:i)));
        end
        LM = LM/min(max(LM));
        nQ = max(find(LM(1:nR)<=0.01));
        %nS = nR -1 + min(find(LM(nR:m-N+1)<=0.01));
        nTotal = nR -1 + min(find(LM(nR:m-N+1)<=0.01));
        nS = nTotal - round(Wt/(x(2)-x(1)));
        %nS = nQ + nQS; %min(find(x>= round(x(nTotal) - Wt)));
        figure(2); plot(x,y,'b',x,dydx.^2/max(dydx.^2),'r-', x,LM,'r',x(nQ),y(nQ),'xk',x(nS),y(nS),'xk')
        
    case 3
        
        % Method 3
        % Hayn, 2006, Comp in Cardiol, Automated QT Interval Measurement from Multilead ECG Signals
        % calculate the range curve for 40ms window or 40 data point
        figure(3); clf;
        f_range = zeros(m,1);
        for i = 2:m-1
            n1 = max(1,i-19);
            n2 = min(m,i+19);
            f_range(i) = max(y(n1:n2))-min(y(n1:n2));
        end
        f_range(1) = f_range(2); f_range(m) = f_range(m-1);
        Th = 0.1*(max(y)-min(y));% set Threshold value
        nr_sub = find(f_range(1:nR-1)<Th); % subthrehold range
        nS1 = max(nr_sub);
        figure(2); 
        Th2 = Th; nS = nS1
        Drg_meanS = (f_range(nS1) + f_range(nS1-1))/(f_range(nS1) + f_range(nS1+1));
        plot(x,y,'b',x,f_range,'r',x(nr_sub), y(nr_sub),'.b',x(nS1), y(nS1),'xk');
        while length(nr_sub)>1
            %Th2 = Th2-min([0.1*Th (max(f_range(nr_sub)) - min(f_range(nr_sub)))/2 ]) %Decreasing the Threshold
            %[Th2-0.1*Th mean(f_range(nr_sub))]
            nrg_max = nr_sub(1) - 1 + find(f_range(nr_sub) ~= max(f_range(nr_sub))); %find max
            nrg_max2 = nrg_max(1) - 1 + max( find(f_range(nrg_max) == max(f_range(nrg_max))) ); %find next max
            [max(f_range(nr_sub)) f_range(nrg_max2)]
            Th2 = mean([max(f_range(nr_sub)) f_range(nrg_max2)]); %Decreasing the Threshold
            
            nr_sub = nS-41 + find(f_range(nS-40:nS-1)<Th2); % subthrehold range
            nS1 = max(nr_sub);
            Drg_mean = (f_range(nS1) + f_range(nS1-1))/(f_range(nS1) + f_range(nS1+1));
            if Drg_mean < Drg_meanS
                nS = nS1; Drg_meanS = Drg_mean;
            end
            hold off; plot(x,y,'b',x,f_range,'r',x(nr_sub), y(nr_sub),'.b',x(nS1), y(nS1),'xk',x(nS), y(nS),'ok',x,Th2,'--k');
            hold on; plot(x(nr_sub), f_range(nr_sub),'.b');
            rg = nr_sub(1) - 20:nS1 + 20;
            axis([x(nr_sub(1) - 20) x(nS1 + 20) min([y(rg); f_range(rg)]) max([y(rg); f_range(rg)])])
            drawnow; pause(0.1) 
        end
        
end

nQR = nR-nQ; nRS = nS-nR;
