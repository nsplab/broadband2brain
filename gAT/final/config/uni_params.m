function [args, T, window, hp_handle, diode_handle, t0] = uni_params(method_name, paramID)
% Returns default parameters to be used for each method

t0 = 0;  % starting point for first interval

args = struct();

if strcmp(method_name, 'int')

    if paramID == 1

        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        % successive integrals
        T = 0.025;
        args.L = 5;
        args.K = 1;
        args.delta_t = get_dt();
        args.numiter = 3;
        args.t_resolution = T / 100;
        args.delta = 0;
        window = 0.005;

    elseif paramID == 2

        % Drew's filter
        [A, B] = butter_hp(1000, 2);
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        % successive integrals
        T = 0.002;
        args.L = 2;
        args.K = 1;
        args.delta_t = get_dt();
        args.numiter = 3;
        args.t_resolution = T / 100;
        args.delta = 0.0025;
        window = 0.005;

    end

elseif strcmp(method_name, 'ker')

    if paramID == 1

        % RC kernel with Drew's filter

        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        T = 0.03;
        args.L = 15;
        args.RC = 0.004;
        args.K = 1;
        args.delta_t = get_dt();
        args.numiter = 3;
        args.t_resolution = T / 100;
        args.h = @RC_lowpass;
        args.delta = 0;
        window = 0.1;

    elseif paramID == 2

        % RLC kernel with hp filter

        % Drew's filter
        [A, B] = butter_hp(1000, 2);
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        T = 0.03;
        args.L = 30; % 1 kHz
        args.s1 = -1000; %-400;
        args.s2 = -800; %-500;
        args.K = 2;
        args.delta_t = 0.0005;
        args.numiter = 3;
        args.t_resolution = T / 100;
        args.h = @RLC;
        args.delta = 0.0025;
        window = 0.01;
        
    elseif paramID >= 10 && paramID < 20
        
        % increasing sampling rate
        sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
        
        % RLC kernel with hp filter

        % Drew's filter
        [A, B] = butter_hp(1000, 2);
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        T = 0.03;
        args.L = round(T * sampling_rates(paramID - 9));
        args.s1 = -1000; %-400;
        args.s2 = -800; %-500;
        args.K = 2;
        args.delta_t = 0.0005;
        args.numiter = 3;
        args.t_resolution = T / 100;
        args.h = @RLC;
        args.delta = 0.0025;
        window = 0.01;
    
    elseif paramID >= 20 && paramID < 30
        
        % variety of sampling rates with changing kernel width
        
        % increasing sampling rate
        sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
        
        % RLC kernel with hp filter

        % Drew's filter
        [A, B] = butter_hp(1000, 2);
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        T = 0.03;
        args.L = round(T * sampling_rates(paramID - 19));
        
        % scaling factor for kernel width
        sc = args.L / (800 * T);
        
        args.s1 = -1000 * sc; %-400;
        args.s2 = -800 * sc; %-500;
        args.K = 2;
        args.delta_t = 0.0005;
        args.numiter = 3;
        args.t_resolution = T / 100;
        args.h = @RLC;
        args.delta = 0.0025;
        window = 0.01;
        
    elseif paramID == 30
        
        % high sampling rate, T decreased, L kept constant
        % for comparison with paramID = 14
        
        % Drew's filter
        [A, B] = butter_hp(1000, 2);
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        T = 0.03;
        args.L = round(T * 1000);
        
        sampling_rate = 3000;
        T = args.L / sampling_rate;
        
        args.s1 = -1000; %-400;
        args.s2 = -800; %-500;
        args.K = 2;
        args.delta_t = 0.0005;
        args.numiter = 3;
        args.t_resolution = T / 100;
        args.h = @RLC;
        args.delta = 0.0025;
        window = 0.01;

    elseif paramID >= 31 && paramID < 41
        
        % increasing sampling rate
        sampling_rates = [20, 40, 60, 80, 100, 120, 150, 200, 250, 300];
        
        % RLC kernel with hp filter

        % Drew's filter
        [A, B] = butter_hp(1000, 2);
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        T = 0.03;
        args.L = round(T * sampling_rates(paramID - 30));
        args.s1 = -1000; %-400;
        args.s2 = -800; %-500;
        args.K = 2;
        args.delta_t = 0.0005;
        args.numiter = 3;
        args.t_resolution = T / 100;
        args.h = @RLC;
        args.delta = 0.0025;
        window = 0.01;
    end
    
elseif strcmp(method_name, 'gibbs')
    
    % gibbs sampling
    
    if paramID >= 10 && paramID < 20
        
        % standard params, pretty low numiter
        
        % increasing sampling rate
        sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];

        % Drew's filter
        [A, B] = butter_hp(1000, 2);  % TODO: choose filter
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;  % TODO: choose diode

        T = 0.03;
        args.L = round(T * sampling_rates(paramID - 9));
        %args.s1 = -1000; %-400;
        %args.s2 = -800; %-500;
        args.K = 2;
        %args.delta_t = 0.0005;
        args.numiter = 5;
        args.mean_iter = 2;  % should really be like 50 and 10 but need performance
        %args.t_resolution = T / 100;  % TODO: add this variable?
        %args.h = @RLC;
        %args.delta = 0.0025;
        window = 0.01;
        args.sigmah = 0.001;  % gaussian kernel width
        
    elseif paramID >= 20 && paramID < 30
        
        % gibbs with no rectification, still small numiter
        
        % increasing sampling rate
        sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
        
        % gibbs sampling

        % Drew's filter
        [A, B] = butter_hp(1000, 2);  % TODO: choose filter
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = 0; %@ideal_diode;  % TODO: choose diode

        T = 0.03;
        args.L = round(T * sampling_rates(paramID - 19));
        %args.s1 = -1000; %-400;
        %args.s2 = -800; %-500;
        args.K = 2;
        %args.delta_t = 0.0005;
        args.numiter = 5;
        args.mean_iter = 2;  % should really be like 50 and 10 but need performance
        %args.t_resolution = T / 100;  % TODO: add this variable?
        %args.h = @RLC;
        %args.delta = 0.0025;
        window = 0.01;
        args.sigmah = 0.001;  % gaussian kernel width
        
    elseif paramID >= 30 && paramID < 40
        
        % gibbs with rectification, more iterations
        
        % increasing sampling rate
        sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
        
        % gibbs sampling

        % Drew's filter
        [A, B] = butter_hp(1000, 2);  % TODO: choose filter
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;  % TODO: choose diode

        T = 0.03;
        args.L = round(T * sampling_rates(paramID - 29));
        %args.s1 = -1000; %-400;
        %args.s2 = -800; %-500;
        args.K = 2;
        %args.delta_t = 0.0005;
        args.numiter = 20;
        args.mean_iter = 6;  % should really be like 50 and 10 but need performance
        %args.t_resolution = T / 100;  % TODO: add this variable?
        %args.h = @RLC;
        %args.delta = 0.0025;
        window = 0.01;
        args.sigmah = 0.001;  % gaussian kernel width
        
    elseif paramID >= 40 && paramID < 50
        
        % conference: gibbs with no rectification, high iter
        
        % increasing sampling rate
        sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
        
        % gibbs sampling

        % Drew's filter
        [A, B] = butter_hp(1000, 2);  % TODO: choose filter
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = 0; %@ideal_diode;  % TODO: choose diode

        T = 0.03;
        args.L = round(T * sampling_rates(paramID - 39));
        %args.s1 = -1000; %-400;
        %args.s2 = -800; %-500;
        args.K = 2;
        %args.delta_t = 0.0005;
        args.numiter = 20;
        args.mean_iter = 6;  % should really be like 50 and 10 but need performance
        %args.t_resolution = T / 100;  % TODO: add this variable?
        %args.h = @RLC;
        %args.delta = 0.0025;
        window = 0.01;
        args.sigmah = 0.001;%0.001;  % gaussian kernel width
        
    elseif paramID >= 50 && paramID < 60
        
        % gibbs with rectification, more iterations
        
        % increasing sampling rate
        sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
        
        % gibbs sampling

        % Drew's filter
        [A, B] = butter_hp(1000, 2);  % TODO: choose filter
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;  % TODO: choose diode

        T = 0.03;
        args.L = round(T * sampling_rates(paramID - 49));
        %args.s1 = -1000; %-400;
        %args.s2 = -800; %-500;
        args.K = 2;
        %args.delta_t = 0.0005;
        args.numiter = 50;
        args.mean_iter = 25;
        %args.t_resolution = T / 100;  % TODO: add this variable?
        %args.h = @RLC;
        %args.delta = 0.0025;
        window = 0.01;
        args.sigmah = 0.001;  % gaussian kernel width
        
    end
    
elseif strcmp(method_name, 'RMSE') || strcmp(method_name, 'RMSE_cal')
    
    if paramID >= 10 && paramID < 20
        
        % increasing sampling rate
        sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
        
        % RLC kernel with hp filter

        % Drew's filter
        [A, B] = butter_hp(1000, 2);
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        T = 0.03;
        args.L = round(T * sampling_rates(paramID - 9));
        args.s1 = -1000; %-400;
        args.s2 = -800; %-500;
        args.K = 2;
        args.delta_t = 0.0005;
        args.numiter = 3;
        args.t_resolution = T / 100;
        args.h = @RLC;
        args.delta = 0.0025;
        window = 0.01;
        
    elseif paramID >= 20 && paramID < 30
        
        % no rectification
        
        % increasing sampling rate
        sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
        
        % RLC kernel with hp filter

        % Drew's filter
        [A, B] = butter_hp(1000, 2);
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = 0; %@ideal_diode;

        T = 0.03;
        args.L = round(T * sampling_rates(paramID - 19));
        args.s1 = -1000; %-400;
        args.s2 = -800; %-500;
        args.K = 2;
        args.delta_t = 0.0005;
        args.numiter = 3;
        args.t_resolution = T / 100;
        args.h = @RLC;
        args.delta = 0.0025;
        window = 0.01;
        
    elseif paramID >= 30 && paramID < 34
        
        % increasing number of iterations
        numiters = [3, 5, 10, 20];
        sampling_rate = 1000;
        
        % RLC kernel with hp filter

        % Drew's filter
        [A, B] = butter_hp(1000, 2);
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        T = 0.03;
        args.L = round(T * sampling_rate);
        args.s1 = -1000; %-400;
        args.s2 = -800; %-500;
        args.K = 2;
        args.delta_t = 0.0005;
        args.numiter = numiters(paramID - 29);
        args.t_resolution = T / 100;
        args.h = @RLC;
        args.delta = 0.0025;
        window = 0.01;
                
    end

elseif strcmp(method_name, 'ML')
        
    if paramID >= 10 && paramID < 20
        
        % increasing sampling rate
        sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
        
        % RLC kernel with hp filter

        % Drew's filter
        [A, B] = butter_hp(1000, 2);
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        T = 0.03;
        args.L = round(T * sampling_rates(paramID - 9));
        args.s1 = -1000; %-400;
        args.s2 = -800; %-500;
        args.K = 2;
        args.delta_t = 0.0005;
        args.numiter = 3;
        args.t_resolution = T / 100;
        args.h = @RLC;
        args.delta = 0.0025;
        window = 0.01;
        
    end
    
elseif strcmp(method_name, 'deconv')
        
    if paramID >= 10 && paramID < 20
        
        % increasing sampling rate
        sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
        
        % RLC kernel with hp filter

        % Drew's filter
        [A, B] = butter_hp(1000, 2);
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        T = 0.03;
        args.L = round(T * sampling_rates(paramID - 9));
        args.K = 2;
        args.h = @RC_lowpass;
        args.RC = 0.001;
        window = 0.01;
    
    end
    
elseif strcmp(method_name, 'analog')  % gAT (new in paper2)
        
    if paramID == 1
    
        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);  % highpass
        args.A_filter = A;
        args.B_filter = B;
        %[A, B] = butter_bp(10, 2000, 2);  % lowpass
        %args.A_lpfilter = A;
        %args.B_lpfilter = B;
        hp_handle = @filter_generic;  %@highpass_lowpass;
        diode_handle = 0;

        T = 0.01;
        args.L = 2;  % number of samples
        args.K = 1;  % number of spikes
        window = 0.02;
        
    elseif paramID == 2
        
        % experimental -- feel free to edit
        
        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);  % highpass
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  %@highpass_lowpass;
        diode_handle = 0;

        T = 0.02;
        args.L = 2;
        args.K = 1;
        window = 0.02;
        
    %elseif paramID >= 10 && paramID <= 20
    elseif paramID >= 10 && paramID <= 25
        
        % different T values
        
        %T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04 0.05];  % Sampling periods
        T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04:0.01:0.10];  % Sampling periods
        T = T_vec(paramID-9);
    
        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);  % highpass
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  %@highpass_lowpass;
        diode_handle = 0;
        
        args.L = 2;
        args.K = 1;
        window = 0.02;
        
    end
    
elseif strcmp(method_name, 'AT')  % analog thresholding
    
    %if paramID >= 10 && paramID <= 20
    if paramID >= 10 && paramID <= 25
        
        % different T values
        
        %T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04 0.05];  % Sampling periods
        T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04:0.01:0.10];  % Sampling periods
        T = T_vec(paramID-9);
    
        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);  % highpass
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  %@highpass_lowpass;
        diode_handle = 0;
        
        args.L = 1;
        args.K = 1;
        window = 0.02;
        
    end
    
elseif strcmp(method_name, 'ann')  % annihilating filer

    if paramID == 1

        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        % successive integrals
        T = 0.01;%0.025;
        args.L = 3;
        args.K = 1;
        %args.delta_t = get_dt();
        %args.numiter = 3;
        %args.t_resolution = T / 100;
        %args.delta = 0;
        window = 0.02;
        
    elseif paramID == 2

        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = 0;

        % successive integrals
        T = 0.01;%0.025;
        args.L = 3;
        args.K = 1;
        %args.delta_t = get_dt();
        %args.numiter = 3;
        %args.t_resolution = T / 100;
        %args.delta = 0;
        window = 0.02;
        
    elseif paramID == 3

        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        % successive integrals
        T = 0.03;%0.025;
        args.L = 5;
        args.K = 2;
        %args.delta_t = get_dt();
        %args.numiter = 3;
        %args.t_resolution = T / 100;
        %args.delta = 0;
        window = 0.02;
        
	elseif paramID == 4

        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = 0;

        % successive integrals
        T = 0.03;%0.025;
        args.L = 5;
        args.K = 2;
        %args.delta_t = get_dt();
        %args.numiter = 3;
        %args.t_resolution = T / 100;
        %args.delta = 0;
        window = 0.02;

    elseif paramID == 5

        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = 0;

        % successive integrals
        T = 0.005;%0.025;
        args.L = 5;
        args.K = 2;
        %args.delta_t = get_dt();
        %args.numiter = 3;
        %args.t_resolution = T / 100;
        %args.delta = 0;
        window = 0.02;
        
	elseif paramID >= 10 && paramID < 20
        
        % increasing sampling rate
        sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
        
        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        % successive integrals
        args.K = 1;
        args.L = 2*args.K+1;
        T = args.L/sampling_rates(paramID-9);
        %args.delta_t = get_dt();
        %args.numiter = 3;
        %args.t_resolution = T / 100;
        %args.delta = 0;
        window = 0.02;
        
	elseif paramID >= 20 && paramID < 23
        
        % lower sampling rate
        
        % increasing sampling rate
        sampling_rates = [1, 25, 100];
        
        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        % successive integrals
        args.K = 1;
        args.L = 2*args.K+1;
        T = args.L/sampling_rates(paramID-19);
        %args.delta_t = get_dt();
        %args.numiter = 3;
        %args.t_resolution = T / 100;
        %args.delta = 0;
        window = 0.02;
        
    end

elseif strcmp(method_name, 'twodelta')

    % Note from Bryan (4/15/13): this part is copied from 'int'
    % I don't know if this is actually a good setup or not.
    %
    % (4/20/13): now based on 'analog'
    % I think that almost all of these values are ignored anyways
    if paramID == 1
    
        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);  % highpass
        args.A_filter = A;
        args.B_filter = B;
        %[A, B] = butter_bp(10, 2000, 2);  % lowpass
        %args.A_lpfilter = A;
        %args.B_lpfilter = B;
        hp_handle = @filter_generic;  %@highpass_lowpass;
        diode_handle = 0;

        T = 0.01;
        args.L = 2;  % number of samples
        args.K = 1;  % number of spikes
        window = 0.02;

    % (4/21/13): also copied from analog
    elseif paramID == 2
        
        % experimental -- feel free to edit
        
        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);  % highpass
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  %@highpass_lowpass;
        diode_handle = 0;

        T = 0.02;
        args.L = 2;
        args.K = 1;
        window = 0.02;
        
    %elseif paramID >= 10 && paramID <= 20
    elseif paramID >= 10 && paramID <= 25
        
        % different T values
        
        %T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04 0.05];  % Sampling periods
        T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04:0.01:0.10];  % Sampling periods
        T = T_vec(paramID-9);
    
        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);  % highpass
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  %@highpass_lowpass;
        diode_handle = 0;
        
        args.L = 2;
        args.K = 1;
        window = 0.02;
        
    end

end
