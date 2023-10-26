function amps = rvalEstimation(data, rawdata, Fs, win)

win_size = round(win*Fs);
n_windows = floor((size(data,2) / win_size));
amps = zeros(2, n_windows);
for wl = [1 2]
    for n = 1:n_windows
        win_idx = (n-1)*win_size+1:n*win_size;
        dc_comp = mean(rawdata(wl, win_idx));
        amps(wl,n) = (max(data(wl,win_idx)) - min(data(wl,win_idx)))/dc_comp;
    end
end
