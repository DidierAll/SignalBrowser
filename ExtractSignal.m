%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                        %%%%%%%%%%
%%%%%%%%%%                ExtractSignal Function                  %%%%%%%%%%
%%%%%%%%%%                                                        %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sout , tout, toutshift] = ExtractSignal(Sin,tin,tS1,tS2,t1,t2);%,SHIFT)

% [Sout , tout, toutshift] = ExtractSignal(Sin,tin,tS1,tS2,t1,t2)
% Remove data segments with artifacts from Sin using artifacts "Signal" segments
% defined by vectors of initial and end time tS1 and tS2 (each index in 
% the pair of data tS1 and tS2 defines an artifact segment between time tS1
% and tS2)
% ExtractSignal returns the truncated signal values Sout and corresponding 
% original timepoints tout after the artifact segments are removed (time 
% discontinuities are present everytime a segment is removed). 
% toutshift returns the new "shifted" time with the time discontinuities 
% in tout removed.

if ~isempty(tS1) && ~isempty(Sin)
    [r c] = size(Sin);
    if r<c
        Sin = Sin';
        tin = tin';
    end
    
    N = length(Sin);
    
    % if narging <=6
    %     SHIFT = 0;
    % end
    if nargin <= 4 || (isempty(t1) && isempty(t2))
        t1 = tin(1); t2 = tin(N); 
    end
    nS_vec = find(tS2>=t1 & tS1<=t2);
    NS = length(nS_vec);
%     nmax = min(find(tin>=t1));
%     nmin = max(find(tin<=t2));
    nmin = min(find(tin>=t1));
    nmax = max(find(tin<=t2));
    %rg = find( tin<=max(t1,tS1(nS_vec(1))) );
    if ~isempty(nS_vec)
        rg = find( tin >= t1 & tin <= tS1(nS_vec(1) ));
    else
        rg = 1:N;
    end
    tout = tin(rg);
    toutshift = tin(rg);
    Sout = Sin(rg);
    if ~isempty(rg)
        n1 = rg(length(rg)); S1 = Sin(n1); tshift = 0;
    else
        n1 = 0; tshift = 0;
    end
    METHOD = 2;
    if METHOD == 1 
        %METHOD 1----------------------------------------------
        for nS = nS_vec
            %extract range
            if nS < nS_vec(NS)
                rg = find(tin>=tS2(nS) & tin<=tS1(nS+1));
            else
                %rg = find(tin>=min(t2,tS2(nS)));
                rg = find(tin>=tS2(nS) & tin<=t2 );
            end
            
            if nS > max(nS_vec)-4
                'ok';
            end
            nout = length(tout);
            if ~isempty(rg)
                n2 = rg(1); S2 = Sin(n2);
                
                if ~isempty(rg) & n2>n1+1 %make sure that there is a discontinuity
                    if n1>nmin & n2 <nmax %avoid adding values for first and last datapoint
                        S_between = (S1+S2)/2;
                        % shift due to introducing a new point (uneven sampling)
                        % technically t(n1) = t(n1-1) + S(n1)
                        %tshift_between = (tin(n1)-tin(n1-1) + tin(n2)-tin(n2-1))/2;
                        tshift_between = (tin(n1+1)-tin(n1) + tin(n2+1)-tin(n2))/2;
                        t_between = tout(nout)+ tshift_between;
                        t_between = tout(nout)+ S_between/1000;
                    else
                        t_between = []; S_between = []; tshift_between = 0;%[];
                    end
                    
                    tout = [tout;t_between; tin(rg)]; % noshift
                    
                    %negative shift due to deleting a  huge shunk
                    %tshift_new = tin(n2)-tin(n1)- (n1~=0)*(tin(n2)-tin(n2-1));
                    try
                        tshift_new = tin(n2)-tin(n1+1);% + (n1==0)*(tin(n1+1)-tin(n1));
                    catch
                        pause(1)
                    end
                    %total shift
                    tshiftold = tshift;
                    tshift = tshift + tshift_between - tshift_new;
                    toutshift = [toutshift ; t_between + tshiftold ; tin(rg)+tshift];
                    
                    %add a value between discontinuity
                    Sout = [Sout; S_between;Sin(rg)];
                end
                n1 = rg(length(rg)); S1 = Sin(n1);
            end
            nS=nS+1;
        end
    else
        %METHOD 2----------------------------------------------
        for nS = nS_vec
            %extract range
            if nS < nS_vec(NS)
                rg = find(tin>=tS2(nS) & tin<=tS1(nS+1));
            else
                %rg = find(tin>=min(t2,tS2(nS)));
                rg = find(tin>=tS2(nS) & tin<=t2 );
                
            end
            if nS > max(nS_vec)-4
                'ok';
            end
            nout = length(tout);
            if ~isempty(rg)
                n2 = rg(1); S2 = Sin(n2);
                
                if ~isempty(rg) & n2>n1+1 %make sure that there is a discontinuity
                    if n1>nmin & n2 <nmax %add values only for other than first and last datapoint
                        S_between = (S1+S2)/2;
                        % shift due to introducing a new point (uneven sampling)
                        % technically t(n1) = t(n1-1) + S(n1)
                        % tshift_between = (tin(n1)-tin(n1-1) + tin(n2)-tin(n2-1))/2;
                        t_between = tout(nout)+ S_between/1000; %S_between = interbeat intervals
                    else
                        t_between = []; S_between = []; 
                    end
                    
                    tout = [tout;t_between; tin(rg)]; % noshift & add in between point
                    
                    %add a value between discontinuity
                    Sout = [Sout; S_between;Sin(rg)];
                end
                n1 = rg(length(rg)); S1 = Sin(n1);
            end
            nS=nS+1;
        end
        
        
        
        % recompute toutshift based on Sout (interbeat intervals)
        toutshift = tout*0;
        toutshift(1) = tout(1);
        for i = 2:length(tout);
            toutshift(i) = toutshift(i-1) + Sout(i)/1000; 
        end
   end
   if r<c
        Sout = Sout';
        tout = tout';
        toutshift = toutshift';
   end
else
    Sout = Sin;
    tout = tin;
    toutshift = tin;
end