clear
load("Finger_SpO2sim_0120.mat")

%% Plot Raw Data
% The full data uses 8 wavelengths. We choose two from them (730 nm and 830 nm)
data = data(:,[1,7])';

%% Filter the data
[b,a] = butter(2, [0.66,15]/(Fs/2));
filtData = filtfilt(b,a,data')';
% Extract the slow drifting
[b2,a2] = butter(2, 0.66/(Fs/2));
slow = filtfilt(b2,a2,data')';

%% Peak-to-Peak Amplitude Estimation and Calculating r-values
amps=[];
win=2;  % Window size in seconds
zero_noise_amps = rvalEstimation(filtData, data, Fs, win);

%% Noisy Data and Bias Estimation 

vars_730=linspace(10^0,10^4,1001);
% vars_730=0:1:50;
% range=linspace(4,1,13);
snr730_div_snr830 = [4,3,2,1,0.5,1/3,1/4];
 
len_var = length(vars_730);
len_div = length(snr730_div_snr830);
total_repeat = 100;
rvals = zeros(len_div, len_var, total_repeat);
rvals_var = zeros(len_div, len_var, total_repeat);
for repeat = 1:total_repeat
    fprintf('Repetition %d\n', repeat);
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

save('r_results', 'rvals', 'rvals_bias', 'rvals_var', 'zero_noise_rvals');

%% Plot the bias
pos = [1070,320,850,700];
figure
set(gcf, 'Position', pos)
colorwheel = parula;
newcolors = colorwheel(round(linspace(1,256,7)),:);
hold on
for i=1:length(snr730_div_snr830)
    plot(sqrt(vars_730), mean(rvals_bias(i,:,:),3), 'LineWidth', 2, 'Color', newcolors(i,:))
end

xlabel('Noise StDev of 730nm Wave ($$\Delta_1$$)','Interpreter','Latex')
ylabel('$$E[\hat{R}] - R$$','Interpreter','Latex')
hcb=colorbar('Ticks',[0,0.25,0.5,0.75,1],...
         'TickLabels',{'4','5/2','1','2/5','1/4'});
c = get(hcb,'Title');
set(c,'String','$\frac{SNR_1}{SNR_2}$','Interpreter','Latex');
set(gca,"FontSize",20)

% saveas(gcf,'bias.pdf');
%% Plot the variance
figure
set(gcf, 'Position', pos)
newcolors = colorwheel(round(linspace(1,256,7)),:);
hold on
for i=1:length(snr730_div_snr830)
    plot(sqrt(vars_730), mean(rvals_var(i,:,:),3), 'LineWidth', 2, 'Color', newcolors(i,:))
end

xlabel('Noise StDev of 730nm Wave ($$\Delta_1$$)','Interpreter','Latex')
ylabel('$$Var[\hat{R}]$$','Interpreter','Latex')
hcb=colorbar('Ticks',[0,0.25,0.5,0.75,1],...
         'TickLabels',{'4','5/2','1','2/5','1/4'});
c = get(hcb,'Title');
set(c,'String','$\frac{SNR_1}{SNR_2}$','Interpreter','Latex');
set(gca,"FontSize",20)

%% Plot bias-over-stdev
figure, hold on
set(gcf, 'Position', pos)
newcolors = colorwheel(round(linspace(1,256,7)),:);
for i=1:length(snr730_div_snr830)
    plot(sqrt(vars_730), mean(rvals_bias(i,:,:),3)./sqrt(mean(rvals_var(i,:,:),3)), 'LineWidth', 2, 'Color', newcolors(i,:))
end

xlabel('Noise StDev of 730nm Wave ($$\Delta_1$$)','Interpreter','Latex')
ylabel('$$Bias[\hat{R}]/StDev(\hat{R})$$','Interpreter','Latex')
hcb=colorbar('Ticks',[0,0.25,0.5,0.75,1],...
         'TickLabels',{'4','5/2','1','2/5','1/4'});
c = get(hcb,'Title');
set(c,'String','$\frac{SNR_1}{SNR_2}$','Interpreter','Latex');
set(gca,"FontSize",20)
