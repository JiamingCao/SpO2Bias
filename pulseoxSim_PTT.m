clear

for subject=1:22
    disp(subject)
    tmp = readmatrix(['PulseTransitTimeData/s',num2str(subject),'_sit.csv'],'NumHeaderLines',1);
    
    % Use red and ir at distal location
    data = tmp(:,[5,4])';
    data = resample(data',100,500)'; % original Fs=500Hz; downsample to 50
    data = data(:, 20:end-20);
    Fs = 100;
    
    %% Filter the data
    [b,a] = butter(2, [0.66,5]/(Fs/2));
    filtData = filtfilt(b,a,data')';
    % Extract the slow drifting
    [b2,a2] = butter(2, 0.66/(Fs/2));
    slow = filtfilt(b2,a2,data')';
    
    I_ac1 = mean(findpeaks(filtData(1,:),'MinPeakDistance',0.8*Fs))+mean(findpeaks(-filtData(1,:),'MinPeakDistance',0.8*Fs));
    filtData = filtData/I_ac1;
    slow = slow/I_ac1;
    
    %% Peak-to-Peak Amplitude Estimation and Calculating r-values
    amps=[];
    win=2;  % Window size in seconds
    zero_noise_amps = rvalEstimation(filtfilt(b,a,filtData'+slow')', data, Fs, win);
    
    %% Noisy Data and Bias Estimation 
    
    % vars_730=linspace(10^0,10^4,1001);
    vars_730=0:0.001:0.1;
    % range=linspace(4,1,13);
    snr730_div_snr830 = [4,3,2,1,0.5,1/3,1/4];
     
    len_var = length(vars_730);
    len_div = length(snr730_div_snr830);
    total_repeat = 50;
    rvals = zeros(len_div, len_var, total_repeat);
    rvals_var = zeros(len_div, len_var, total_repeat);
    for repeat = 1:total_repeat
        % fprintf('Repetition %d\n', repeat);
        amps = zeros(length(vars_730), length(snr730_div_snr830), size(zero_noise_amps,1), size(zero_noise_amps,2));
        parfor i=1:len_var
            v = vars_730(i);
            for j=1:len_div
                snr_ratio = snr730_div_snr830(j);
                std_unif_730 = (rand(1,size(data,2)) - 0.5) .* sqrt(12);
                std_unif_830 = (rand(1,size(data,2)) - 0.5) .* sqrt(12);
        
                noise_730 = std_unif_730 .* sqrt(v);
        
                snr_730 = (mean(findpeaks(filtData(1,:),'MinPeakDistance',0.8*Fs))+mean(findpeaks(-filtData(1,:),'MinPeakDistance',0.8*Fs)))^2 ./ mean(noise_730.^2,2);
                snr_830 = snr_730 ./ snr_ratio;
        
                v_830 = (mean(findpeaks(filtData(2,:),'MinPeakDistance',0.8*Fs))+mean(findpeaks(-filtData(2,:),'MinPeakDistance',0.8*Fs)))^2 ./ snr_830;
                noise_830 = std_unif_830 .* sqrt(v_830);
        
                noiseData = filtData + [noise_730; noise_830] + slow;
        
                noise_amps=rvalEstimation(filtfilt(b,a,noiseData')', noiseData, Fs, win);
                amps(i,j,:,:) = noise_amps;
            end
        end
        
        rvals(:,:,repeat) = transpose(mean(squeeze(amps(:,:,1,:) ./ amps(:,:,2,:)),3));
        rvals_var(:,:,repeat) = transpose(var(squeeze(amps(:,:,1,:) ./ amps(:,:,2,:)),0,3));
    end
    zero_noise_rvals = mean(zero_noise_amps(1,:) ./ zero_noise_amps(2,:),2);
    rvals_bias = rvals - zero_noise_rvals;
    
    all_bias(:,:,subject) = mean(rvals_bias,3);
    all_var(:,:,subject) = mean(rvals_var,3);
    all_zeror(subject) = zero_noise_rvals;
    all_zerovar(subject) = var(zero_noise_amps(1,:) ./ zero_noise_amps(2,:),0,2);
end
idx = all_zerovar<0.01;
% save('r_results', 'rvals', 'rvals_bias', 'rvals_var', 'zero_noise_rvals');

%% Plot the bias
pos = [1070,320,850,700];
figure
set(gcf, 'Position', pos)
colorwheel = parula;
newcolors = colorwheel(round(linspace(1,256,7)),:);
hold on
for i=1:length(snr730_div_snr830)
    shadedErrorBar(sqrt(vars_730), mean(all_bias(i,:,idx),3), std(all_bias(i,:,idx),0,3)/sqrt(sum(idx)),'lineprops',{ 'LineWidth', 2, 'Color', newcolors(i,:)})
end

xlabel('Noise StDev of 730nm Wave ($$\Delta_1$$)','Interpreter','Latex')
ylabel('$$E[\hat{R}] - R$$','Interpreter','Latex')
hcb=colorbar('Ticks',[0,0.5,1,1.5,2],...
         'TickLabels',{'4','5/2','1','2/5','1/4'});
c = get(hcb,'Title');
set(c,'String','$\frac{SNR_1}{SNR_2}$','Interpreter','Latex');
set(gca,"FontSize",20)
xlim([0, 0.3])

% saveas(gcf,'bias.pdf');
%% Plot the variance
figure
set(gcf, 'Position', pos)
newcolors = colorwheel(round(linspace(1,256,7)),:);
hold on
for i=1:length(snr730_div_snr830)
    shadedErrorBar(sqrt(vars_730), mean(all_var(i,:,idx),3), std(all_var(i,:,idx),0,3)/sqrt(sum(idx)),'lineprops',{'LineWidth', 2, 'Color', newcolors(i,:)})
end

xlabel('Noise StDev of 730nm Wave ($$\Delta_1$$)','Interpreter','Latex')
ylabel('$$Var[\hat{R}]$$','Interpreter','Latex')
hcb=colorbar('Ticks',[0,0.5,1,1.5,2],...
         'TickLabels',{'4','5/2','1','2/5','1/4'});
c = get(hcb,'Title');
set(c,'String','$\frac{SNR_1}{SNR_2}$','Interpreter','Latex');
set(gca,"FontSize",20)
xlim([0, 0.3])

