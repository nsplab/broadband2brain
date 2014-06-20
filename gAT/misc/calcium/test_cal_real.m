% Load real data
load('Imaging-SNR-Data');
dat = ExpDat20080716a24;

F = dat.Fluorescence;
V.T = length(F);
V.dt = mean(diff(dat.FluorescenceTime));