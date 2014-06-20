function [ plo, phi, pdlo, pdhi, distrlo, distrhi, datlo, dathi, cutoff, pd_glm ] = nlx_spikephase_2(elec, t, name, cutoff)
% function nlx_spikephase
%   Compute spike phase locking from Neuralynx recordings.

% elec is channel
% t is spike times
% name has fields data_set, method_name, T, title (of plot)
% cutoff is threshold for high/low amplitude, optional (default is median)

% NOTES
% distrhi is currently normalized to have mean spk_rate
% should it be normalized to have mean 1
% used by fig4a.m and fig4.m

% parameters
%sdr = [root_dir() '../Data/matthew/2008-07-21_11-09-57/']; % prepared session directory
fdr = [root_dir() '../Data/']; % raw data directory
monkey = 'matthew'; % monkey
session = 20080721; % recording session
%elec = 16; % electrode
%clus = 1; % cluster
rpc = 5000; % records to look at per cycle (512 samples/record) to prevent out-of-memory errors
blfp = [20 30]; % lfp band (Hz)
rems = 1; % remove spikes? (0 or 1)
bspk = [100 6000]; % spike band (Hz)
wspk = [-0.5 2]; % spike window (ms)

% load spike times
cspk = (round((t*1e6-3)/33)*33+3)/1e6;  % make cspk e6 congruent to 3 (mod 33)
PHASE = NaN*ones(size(cspk));
AMP = NaN*ones(size(cspk));
fs = 30303; % assumed sampling rate (S/s)

% convert spike times to vector (needs to be saved for glmfit for fig 4)
% spike times are already congruent to 3 (mod 33)
% make them congruent to 3 (mod 990)
% 990 = 33*30
% this is ~ 1 ms bins
% also only keep lfp values at times congruent to 3 (mod 990)
%spk_rounded = (round((t*1e6-3)/990)*990+3);  % these have a factor of 10^6 (they are integers)
time_vec = [];
lfp_vec = [];

spk_rate = length(t)/(t(end)-t(1));  % average firing rate

% compute average spike waveform
if rems
    cd(fdr);
    swin = round(wspk(1)/1000*fs):round(wspk(2)/1000*fs);
    twin = swin/fs*1000; Nwin = length(swin);
    TSWP = NaN*ones(length(cspk),Nwin);
    [b,a] = butter(2,bspk/(fs/2));
    scsess = num2str(session);
    cnlx = dir([monkey '/' scsess(1:4) '-' scsess(5:6) '-' scsess(7:8) '_*']);
    if isempty(cnlx), disp('no neuralynx files for this session were found'); return; end
    endTS = 0; fprintf('computing average spike waveform...');
    for ifile = 1:length(cnlx) % in case more than one set of CSC files was used during the recording session
        cd([monkey '/' cnlx(ifile).name]);
        idx = [1 rpc];
        [allTStamps,NlxHeader] = Nlx2MatCSC(['CSC' num2str(elec) '.ncs'],[1 0 0 0 0],1,1);
        nrec = length(allTStamps); iiBV = strmatch('-ADBitVolts',NlxHeader);
        ad2uv = str2num(NlxHeader{iiBV}(strfind(NlxHeader{iiBV},' '):end))*1e6; % AD units to micro-volts
        if ifile ~= 1, endTS = endTS - allTStamps(1) + 1e6; end % 1 sec gap in timestamps between files (this is what was done for spike timestamps)
        while idx(1) <= nrec - rpc  % big fix
            if idx(1) == 1
                fprintf('%03d%%',round(idx(1)/nrec*100));
            else
                fprintf('\b\b\b\b%03d%%',round(idx(1)/nrec*100));
            end
            [TStamps,ChannelNumbers,SampleFrequencies,NumberValidSamples,Samples] = Nlx2MatCSC(['CSC' num2str(elec) '.ncs'],ones(1,5),0,2,idx);
            tlfp = zeros(prod(size(Samples)),1); lastsamp = 0;
            for isamp = 1:length(TStamps)
                tlfp(lastsamp+1:lastsamp+NumberValidSamples(isamp)) = TStamps(isamp)+round(1e6/SampleFrequencies(isamp))*[0:NumberValidSamples(isamp)-1]';
                lastsamp = lastsamp + NumberValidSamples(isamp);
            end
            tlfp = tlfp + endTS; % correct timestamp for multiple files
            dlfp = reshape(Samples,sum(NumberValidSamples),1)*ad2uv;
            if session < 20080507, dlfp = -1*dlfp; end % Neuralynx recordings had flipped polarity prior to 20080507 (pre Cheetah v5.1)
            idx = idx + rpc;
            dlfp = filtfilt(b,a,dlfp);
            [tf,loc] = ismember(round(cspk*1e6),tlfp);
            ispk = find(tf); sspk = loc(ispk);
            if ~isempty(ispk)
                mspk = sspk*ones(1,Nwin) + ones(length(sspk),1)*swin;
                [iskp,jskp] = find(mspk<1 | mspk>length(tlfp)); iskp = unique(iskp);
                ispk(iskp) = []; mspk(iskp,:) = [];
                TSWP(ispk,:) = dlfp(mspk);
            end
        end
        endTS = endTS + allTStamps(end);
        cd(fdr);
    end
    fprintf('\n');
    inan = find(isnan(sum(TSWP,2)));
    TSWP(inan,:) = [];
    msw = mean(TSWP);
end

% compute spike phase locking
cd(fdr);
scsess = num2str(session);
cnlx = dir([monkey '/' scsess(1:4) '-' scsess(5:6) '-' scsess(7:8) '_*']);
if isempty(cnlx), disp('no neuralynx files for this session were found'); return; end
endTS = 0; fprintf('computing spike phase distributions...');
for ifile = 1:length(cnlx) % in case more than one set of CSC files was used during the recording session
    cd([monkey '/' cnlx(ifile).name]);
    idx = [1 rpc];
    [allTStamps,NlxHeader] = Nlx2MatCSC(['CSC' num2str(elec) '.ncs'],[1 0 0 0 0],1,1);
    nrec = length(allTStamps); iiBV = strmatch('-ADBitVolts',NlxHeader);
    ad2uv = str2num(NlxHeader{iiBV}(strfind(NlxHeader{iiBV},' '):end))*1e6; % AD units to micro-volts
    if ifile ~= 1, endTS = endTS - allTStamps(1) + 1e6; end % 1 sec gap in timestamps between files (this is what was done for spike timestamps)
    while idx(1) <= nrec - rpc  % big fix
        if idx(1) == 1
            fprintf('%03d%%',round(idx(1)/nrec*100));
        else
            fprintf('\b\b\b\b%03d%%',round(idx(1)/nrec*100));
        end
        [TStamps,ChannelNumbers,SampleFrequencies,NumberValidSamples,Samples] = Nlx2MatCSC(['CSC' num2str(elec) '.ncs'],ones(1,5),0,2,idx);
        tlfp = zeros(prod(size(Samples)),1); lastsamp = 0;
        for isamp = 1:length(TStamps)
            tlfp(lastsamp+1:lastsamp+NumberValidSamples(isamp)) = TStamps(isamp)+round(1e6/SampleFrequencies(isamp))*[0:NumberValidSamples(isamp)-1]';
            lastsamp = lastsamp + NumberValidSamples(isamp);
        end
        tlfp = tlfp + endTS; % correct timestamp for multiple files
        dlfp = reshape(Samples,sum(NumberValidSamples),1)*ad2uv;
        if session < 20080507, dlfp = -1*dlfp; end % Neuralynx recordings had flipped polarity prior to 20080507 (pre Cheetah v5.1)
        idx = idx + rpc;
        fs = mean(SampleFrequencies);
        
        % remove spike waveforms (similar to Zanos et al 2012)
        if rems
            [tf,loc] = ismember(round(cspk*1e6),tlfp);
            ispk = find(tf); sspk = loc(ispk);
            if ~isempty(ispk)
                mspk = sspk*ones(1,Nwin) + ones(length(sspk),1)*swin;
                [iskp,jskp] = find(mspk<1 | mspk>length(tlfp)); iskp = unique(iskp);
                ispk(iskp) = []; mspk(iskp,:) = [];
                c = zeros(length(ispk),1);
                for i = 1:length(ispk)
                    c(i) = msw'\dlfp(mspk(i,:)); % scale factor that minimized sum squared differences between the two (Tolias et al 2007)
                end
                dlfp(mspk) = dlfp(mspk) - (c*ones(1,Nwin)).*(ones(size(mspk,1),1)*msw);
            end
        end
        
        % compute analytic signal for specified frequency band
        [b,a] = butter(2,blfp/(fs/2));
        dlfp = filtfilt(b,a,dlfp);
        alfp = hilbert(dlfp);
        
        % phase/amplitude of LFP oscillation at spike times
        [tf,loc] = ismember(round(cspk*1e6),tlfp);
        ispk = find(tf); sspk = loc(ispk);
        if ~isempty(ispk)
            PHASE(ispk) = angle(alfp(sspk));
            AMP(ispk) = abs(alfp(sspk));
        end
        
        
        % TESTING
        %{
        figure;
        hold on;
        %plot(tlfp, dlfp, 'b');  % I think this plots LFP
        plot(tlfp, angle(alfp), 'b');
        stem(cspk*1e6, median(angle(alfp))*ones(size(cspk)), 'r');  % I think this plots spikes
        title('spk-LFP plot');
        %}
        
        % for glmfit (for fig 4)
        % take only lfp values angle(alfp) with time tlfp that is congruent
        % to 3 (mod 990)
        first_ind = 31 - (mod(tlfp(1), 990) - 3)/33;
        %test1 = mod(tlfp(first_ind), 990)
        time_ind = first_ind : 30 : length(tlfp);
        time_vec = [time_vec ; tlfp(time_ind)];
        lfp_vec = [lfp_vec ; angle(alfp(time_ind))];
        %test1 = mod(time_vec, 990)
        %test3 = 0
        
        % correct for LFP waveform asymmetry (Siapas et al 2005)
        ECDFt = sort(angle(alfp));
        ECDFx = (1:length(ECDFt))/length(ECDFt);
        
        % Bug fix - remove duplicates
        dup_ind = find(diff(ECDFt) == 0);
        ECDFt(dup_ind) = [];
        ECDFx(dup_ind) = [];
        
        cx = interp1(ECDFt,ECDFx,PHASE(i),'linear');
        PHASE(i) = 2*pi*cx-pi;
    end
    endTS = endTS + allTStamps(end);
    cd(fdr);
end
fprintf('\n');
in = isnan(AMP);
AMP(in) = [];
PHASE(in) = [];

% TESTING: plot AMP
%figure;
%hist(AMP,100)

% plot
bins = -pi:pi/8:pi;
if nargin < 4
    cutoff = median(AMP);
end
ilo = AMP < cutoff;

% TESTING
frac_lo = (nnz(ilo))/length(AMP)

nlo = histc(PHASE(ilo),bins); nlo(end) = [];
plo = rayleightest(PHASE(ilo)); % valid only if correct for waveform asymmetry
pdlo = atan2(sum(sin(PHASE(ilo))),sum(cos(PHASE(ilo))))*180/pi;
ihi = AMP > cutoff;
nhi = histc(PHASE(ihi),bins); nhi(end) = [];
phi = rayleightest(PHASE(ihi));
pdhi = atan2(sum(sin(PHASE(ihi))),sum(cos(PHASE(ihi))))*180/pi;
deg = (bins(1:end-1)+mean(diff(bins))/2)*180/pi;
fig = figure('Units','normalize','Position',[1/4 1/4 1/2 1/2]);
distrlo = nlo/(sum(nlo)/(length(bins)-1));
distrhi = nhi/(sum(nhi)/(length(bins)-1));
subplot(1,2,1), plot(deg,distrlo,'k','LineWidth',2); hold on; axis square; ylim1 = get(gca,'YLim');
subplot(1,2,2), plot(deg,distrhi,'k','LineWidth',2); hold on; axis square; ylim2 = get(gca,'YLim');
subplot(1,2,1), set(gca,'Box','off','XLim',[-180 180],'YLim',[min(ylim1(1),ylim2(1)), max(ylim1(2),ylim2(2))]);
line(get(gca,'XLim'),ones(1,2),'Color','k','LineStyle',':');
xlabel('Phase (degrees)', 'FontSize', 14);
%xlabel(['Phase of ' num2str(blfp(1)) '-' num2str(blfp(2)) 'Hz Oscillations (deg)'], 'FontSize', 14);
ylabel('Normalized P(spike | phase)', 'FontSize', 14); title(['Low Amplitude (<' num2str(cutoff,'%5.1f') '\muV)'], 'FontSize', 14);
%text(min(get(gca,'XLim')),max(get(gca,'YLim')),{[' N = ' num2str(sum(nlo))],[' p = ' num2str(plo)],[' pd = ' num2str(pdlo,'%5.1f') 'deg']},'HorizontalAlignment','left','VerticalAlignment','top');
set(gca, 'FontSize', 12);
subplot(1,2,2), set(gca,'Box','off','XLim',[-180 180],'YLim',[min(ylim1(1),ylim2(1)), max(ylim1(2),ylim2(2))]);
line(get(gca,'XLim'),ones(1,2),'Color','k','LineStyle',':');
xlabel('Phase (degrees)', 'FontSize', 14);
%xlabel(['Phase of ' num2str(blfp(1)) '-' num2str(blfp(2)) 'Hz Oscillations (deg)'], 'FontSize', 14);
ylabel('Normalized P(spike | phase)', 'FontSize', 14); title(['High Amplitude (>' num2str(cutoff,'%5.1f') '\muV)'], 'FontSize', 14);
%text(min(get(gca,'XLim')),max(get(gca,'YLim')),{[' N = ' num2str(sum(nhi))],[' p = ' num2str(phi)],[' pd = ' num2str(pdhi,'%5.1f') 'deg']},'HorizontalAlignment','left','VerticalAlignment','top');
set(gca, 'FontSize', 12);
if nargin < 3
    suptitle([monkey ' ' num2str(session) ', ' num2str(elec) ', ' num2str(clus)]);
else
    suptitle(name.title);
end

datlo = PHASE(ilo);
dathi = PHASE(ihi);

direc = [root_dir() '../paper2plots/final/fig4/'];
saveas(fig, [direc 'LFP'], 'epsc');

%print(fig, [root_dir() '../paper2plots/LFP/phase_' int2str(elec)], '-dpdf', '-r0');

% Preferred direction using glmfit
spk_vec = [0 bin_spikes(time_vec, t*1e6)];
vx = cos(lfp_vec); % not actually velocities but analogous to velocities in tuning curve calculation
vy = sin(lfp_vec);
[b, dev, stats] = glmfit([vx vy], spk_vec, 'poisson', 'const', 'on');
pd_glm = atan2(b(3),b(2))*180/pi;

% normalize distrhi to have mean spk_rate
%distrhi = distrhi * spk_rate;

% Save data to MAT file (for figure 4a)
save([root_dir() '../LFP/distr/distr_' name.method_name '_' int2str(name.data_set) '_' int2str(elec) '_' int2str(1000*name.T)], 'deg', 'distrlo', 'distrhi', 'cutoff', 'pdhi', 'pdlo', 'time_vec', 'lfp_vec', 'spk_vec', 'spk_rate');

% TESTING: plot
%{
figure;
hold on;
plot(time_vec, lfp_vec, 'r');
stem(time_vec, spk_vec, 'b');
%}
