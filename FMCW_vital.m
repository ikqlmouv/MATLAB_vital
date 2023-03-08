%clc
%clear all
%close all

% Initialization
chirp_sampling_frequency = 170000; % sampling frequency should be much higher while simuation (Acc to Nyquist Criteria Fs >= 2.fc)
number_of_samples = 128;
resolution = chirp_sampling_frequency / number_of_samples; % resolution = Fs/N, very important parameter to decide unknown frequencies
result_sampling_rate = 20;

number_of_chrip = size(I_data);

% remove DC
for i=1:number_of_chrip(1)

    signal = I_data(i,:);
    signal = signal - mean(signal);
    I_data(i,:) = signal;
    
    signal = Q_data(i,:);
    signal = signal - mean(signal);
    Q_data(i,:) = signal;    

end

number_of_chrip = size(I_data);
tic
rp = plot_phase_data(I_data,Q_data,number_of_chrip(1),chirp_sampling_frequency,1024,result_sampling_rate);
toc

% tic
% [status, cmdout] = system('python segmentation.py data.csv 20')
% cmdout = deblank(cmdout);
% cmdout = regexp(cmdout, '\s+', 'split')
% % HR_ans = cell2mat(cmdout(1))
% % BR_ans = cell2mat(cmdout(2))
% toc


   

function [rp] = plot_phase_data(I_signal_2d,Q_signal_2d,chirp_num,fs,N,result_smapling_rate)

    range_coeff = 3e8 / (4*pi*24.125e9) / 1.5e-3;
    phase_to_plot = zeros(chirp_num,1);
    check_phase = zeros(chirp_num,1);
    indeces_to_plot = zeros(512,1);
    I_range_profile_matrix=zeros(chirp_num,512,1);
    Q_range_profile_matrix=zeros(chirp_num,512,1);
    
         for i=1:chirp_num % for chirp 
            
            I_signal = I_signal_2d(i,:);
            Q_signal = Q_signal_2d(i,:);

            I_signal(N) = 0;
            Q_signal(N) = 0;
            
    %       Compute the FFT signal 
            y = abs(fft(I_signal));   
            f = (0:length(y)-1)*fs/length(y);
            y = y(1:N/2); %remove half len
            f = f(1:N/2); %remove half len
            I_range_profile_matrix(i,:) = y;

            y = abs(fft(Q_signal));   
            f = (0:length(y)-1)*fs/length(y);
            y = y(1:N/2); %remove half len
            f = f(1:N/2); %remove half len
            Q_range_profile_matrix(i,:) = y;
    
    %       Extrack the trrget bin   
            [maxpeak, maxpeakindes] = max(y);
            indeces_to_plot(maxpeakindes) = indeces_to_plot(maxpeakindes) + 1;
    
         end

% %plotting range-profile matrix
% 
%         figure
%         pcolor(I_range_profile_matrix);
%         title('range-profile matrix');
%         ylabel('ramp-num');
%         xlabel('range-bin');

%Check the initial phase
        check_phase = angle(I_range_profile_matrix(:,1));
%         figure 
%         plot(check_phase);
%         title('Phase check');
%

        [maxval,maxindex] = max(indeces_to_plot)
        phase_to_plot = atan(Q_range_profile_matrix(:,maxindex) ./ I_range_profile_matrix(:,maxindex));
        sizep = size(phase_to_plot)

%remove DC
        avg = mean(phase_to_plot);
        for i=1:chirp_num
            phase_to_plot(i) = phase_to_plot(i) - avg;
        end
    
%plotting result
        heart_rate = zeros(chirp_num,1);
        breath_rate = zeros(chirp_num,1);
        time = (0:length(phase_to_plot)-1) / result_smapling_rate;
        phase_to_plot = unwrap(phase_to_plot);

% IIR highpass filter (HR)

        [N,Wc] = cheb2ord(0.8/(result_smapling_rate/2),0.4/(result_smapling_rate/2),3,60)
        [b,a] = cheby2(N,60,Wc,'high');
        heart_rate = filtfilt(b,a,phase_to_plot);

        [N,Wc] = cheb2ord(1.6/(result_smapling_rate/2),3.6/(result_smapling_rate/2),3,60)
        [b,a] = cheby2(N,60,Wc,'low');
        heart_rate = filtfilt(b,a,heart_rate);

%         heart_rate = movmean(phase_to_plot,3);

% IIR lowpass filter (BR)
        
%         [N,Wc] = cheb2ord(0.7/(result_smapling_rate/2),0.5/(result_smapling_rate/2),3,60)
%         [b,a] = cheby2(N,60,Wc,'low');
%         breath_rate = filtfilt(b,a,phase_to_plot);
        

% Moving average lowpass filter (BR)

        breath_rate = movmean(phase_to_plot,23);
     
        %phase_to_plot = range_coeff * phase_to_plot;
        figure
        plot(time,phase_to_plot);
        title('Phase plot in time domain');
        xlabel('time (s)');
        ylabel('phase');
        writematrix(phase_to_plot,'data.csv');

        figure
        plot(time,heart_rate);
        title('Phase plot in time domain (HR)');
        xlabel('time (s)');
        ylabel('phase');
        %writematrix(heart_rate,'data.csv');

        figure
        plot(time,breath_rate);
        title('Phase plot in time domain (BR)');
        xlabel('time (s)');
        ylabel('phase');
        %writematrix(breath_rate,'data.csv');

%plotting FFT result

        figure
        frequency = (0:length(phase_to_plot)-1)*result_smapling_rate/length(phase_to_plot);   %used to frequency domain (FFT result)
        phase_to_plot = abs(fft(heart_rate));
        frequency = frequency(1:length(frequency)/2);
        phase_to_plot = phase_to_plot(1:length(phase_to_plot)/2);
        plot(frequency,phase_to_plot);
        title('FFT result');
        xlabel('frequency');
        ylabel('magnitude');

%music spectrum evaluation
%        pmusic(heart_rate,100,length(heart_rate),20);
        rp = I_range_profile_matrix;



end

