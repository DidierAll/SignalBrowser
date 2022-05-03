%%%%%%%% Default parameter Values  %%%%%%%%%%%

fprintf('Loading Parameters Values \n')
        ud.RRmin = 0.3; %minimum IBI interval for human (mice RRmin = 0.05/0.2)
        %HRV bands 
        ud.VLFmin = 0; 
        ud.VLFmax = 0.04; 
        ud.LFmin = 0.04; 
        ud.LFmax = 0.15; 
        ud.HFmin = 0.15; 
        ud.HFmax = 0.6;
        
        % filters 
        ud.ECGfilter_freq = [50];
        ud.ECGfilter_details = 'low pass butter, 4th order';
        ud.RESPfilter_freq = [0.05 0.5];
        ud.RESPfilter_details = 'band pass butter, 4th order';
        ud.SCLfilter_freq = [1];
        ud.SCLfilter_details = 'low pass butter, 4th order';
        ud.ECGfiltered = 'yes';
        
        
        % HRV values
        %interpolation Frequency = 4Hz
        ud.Fs_interp = 4; 
        %Detrend 2nd order
        ud.detrend = 2;
        ud.NFFT = 2048;
        ud.Nwindow = 128*ud.Fs_interp;
        ud.Noverlap = ud.Nwindow/2;
        ud.FreqMethod = 1; %Welch methods
        ud.Filter = 2; %2 = Hamming, 3 = hanning, 4 = barlett, 5 = blackman, [] = rectangular

        % Spectrogram Values
        ud.SpectrogramMethod = 4;

        
        % EDR
        ud.nQR = 43; ud.nRS = 34;

        % Samples informations
        % this is used to compute metrics on specific data samples
        % segments given a specific experimental protocol . 
        % This can also be defined via the GUI/Toolbox  
        ud.Samplet1 = [];
        ud.Samplet2 = [];
        ud.Nber_of_Samples = 0;
        ud.Sample_Nber = [];
        ud.SampleLabel = {};

%         % examples of data sample computing metrics for the first 5 min  
%         % [0 300] and next 5 min [300 600]
%         ud.Samplet1 = [0 300];
%         ud.Samplet2 = [300 600];
%         ud.Nber_of_Samples = 2;
%         ud.Sample_Nber = 1;
%         ud.SampleLabel = {'PRE','POST'};

        
        %Signal Excluded Samples
        if ~isfield(ud,'Signalt1') || isempty(ud.Signalt1)
            %Signal Excluded Samples
            ud.Signalt1 = [];
            ud.Signalt2 = [];
            ud.Nber_of_Signals = 0;
            ud.Signal_Nber = [];
        end
        
        
                %updating fields
        set(h.VLFmin,'String',sprintf('%4.2f',ud.VLFmin));
        set(h.VLFmax,'String',sprintf('%4.2f',ud.VLFmax));
        set(h.LFmin,'String',sprintf('%4.2f',ud.LFmin));
        set(h.LFmax,'String',sprintf('%4.2f',ud.LFmax));
        set(h.HFmin,'String',sprintf('%4.2f',ud.HFmin));
        set(h.HFmax,'String',sprintf('%4.2f',ud.HFmax));
        
        set(h.interpol_edit,'String',sprintf('%d',ud.Fs_interp));
        set(h.detrend_menu,'Value',ud.detrend);
        set(h.NFFT_edit,'String',sprintf('%d',ud.NFFT));
        set(h.WindowLength_edit,'String',sprintf('%d',ud.Nwindow));
        set(h.WindowOverlap_edit,'String',sprintf('%d',ud.Noverlap));
        set(h.FreqMethod_menu,'Value',ud.FreqMethod);
        
        % Update Sample Edit
        set(h.SampleLength_edit,'String',[]);
        set(h.SampleStart_edit,'String',[]);
        set(h.SampleEnd_edit,'String',[]);
        set(h.Sample_Nber,'String',[]);
        set(h.Nber_of_Samples,'String',[]);
        
        set([ud.SamplePushbuttonHandle;ud.SampleEditHandle],'Enable','off');
        set(h.Nber_of_Samples,'Enable','on');
        
        % Update Measure Edit
        set(h.HRVMeasurePeak_text,'string',[]);
        set(h.HRVMeasurePwr_ms_text,'string',[]);
        set(h.HRVMeasurePwr_perc_text,'string',[]);
        set(h.HRVMeasurePwr_nu_text,'string',[]);
        set(h.TimeMeasure_text,'string',[]);

        set(h.SCLMeasure_text,'string',[])
        set(h.RespirationMeasure_text,'string',[])

