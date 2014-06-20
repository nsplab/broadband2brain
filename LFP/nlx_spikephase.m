function nlx_spikephase(elec, clus, neuron_num)
% function nlx_spikephase
%   Compute spike phase locking from Neuralynx recordings.

% parameters
sdr = [root_dir() '../Data/matthew/2008-07-21_11-09-57/']; % prepared session directory
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
cd(sdr);
if ~exist([monkey num2str(session) '.mat'],'file'), disp('not a valid session'); return; end
load([monkey num2str(session) '.mat'],'spk');
thiselec = []; thisclus = [];
for i = 1:length(spk), thiselec(i) = (spk(i).elec == elec); end;
for i = 1:length(spk), thisclus(i) = (spk(i).clus == clus); end
indspk = find(thiselec & thisclus);
if isempty(indspk) | length(indspk)>1
    disp('electrode/cluster not found');
    return;
else
    cspk = spk(indspk).data;
    celec = spk(indspk).elec;
    cclus = spk(indspk).clus;
end
clear spk
PHASE = NaN*ones(size(cspk));
AMP = NaN*ones(size(cspk));
fs = 30303; % assumed sampling rate (S/s)

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
        
        % correct for LFP waveform asymmetry (Siapas et al 2005)
        ECDFt = sort(angle(alfp));
        ECDFx = (1:length(ECDFt))/length(ECDFt);
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

% plot
bins = -pi:pi/8:pi;
ilo = AMP < median(AMP);
nlo = histc(PHASE(ilo),bins); nlo(end) = [];
plo = rayleightest(PHASE(ilo)); % valid only if correct for waveform asymmetry
pdlo = atan2(sum(sin(PHASE(ilo))),sum(cos(PHASE(ilo))))*180/pi;
ihi = AMP > median(AMP);
nhi = histc(PHASE(ihi),bins); nhi(end) = [];
phi = rayleightest(PHASE(ihi));
pdhi = atan2(sum(sin(PHASE(ihi))),sum(cos(PHASE(ihi))))*180/pi;
deg = (bins(1:end-1)+mean(diff(bins))/2)*180/pi;
fig = figure('Units','normalize','Position',[1/4 1/4 1/2 1/2]);
subplot(1,2,1), plot(deg,nlo/(sum(nlo)/(length(bins)-1)),'k','LineWidth',2); hold on; axis square; ylim1 = get(gca,'YLim');
subplot(1,2,2), plot(deg,nhi/(sum(nhi)/(length(bins)-1)),'k','LineWidth',2); hold on; axis square; ylim2 = get(gca,'YLim');
subplot(1,2,1), set(gca,'Box','off','XLim',[-180 180],'YLim',[min(ylim1(1),ylim2(1)), max(ylim1(2),ylim2(2))]);
line(get(gca,'XLim'),ones(1,2),'Color','k','LineStyle',':');
xlabel(['phase of ' num2str(blfp(1)) '-' num2str(blfp(2)) 'Hz oscillations (deg)']); ylabel('P(spike|phase)'); title(['amp < ' num2str(median(AMP),'%5.1f') '\muV']);
text(min(get(gca,'XLim')),max(get(gca,'YLim')),{[' N = ' num2str(sum(nlo))],[' p = ' num2str(plo)],[' pd = ' num2str(pdlo,'%5.1f') 'deg']},'HorizontalAlignment','left','VerticalAlignment','top');
subplot(1,2,2), set(gca,'Box','off','XLim',[-180 180],'YLim',[min(ylim1(1),ylim2(1)), max(ylim1(2),ylim2(2))]);
line(get(gca,'XLim'),ones(1,2),'Color','k','LineStyle',':');
xlabel(['phase of ' num2str(blfp(1)) '-' num2str(blfp(2)) 'Hz oscillations (deg)']); ylabel('P(spike|phase)'); title(['amp > ' num2str(median(AMP),'%5.1f') '\muV']);
text(min(get(gca,'XLim')),max(get(gca,'YLim')),{[' N = ' num2str(sum(nhi))],[' p = ' num2str(phi)],[' pd = ' num2str(pdhi,'%5.1f') 'deg']},'HorizontalAlignment','left','VerticalAlignment','top');
suptitle([monkey ' ' num2str(session) ', ' num2str(elec) ', ' num2str(clus)]);

print(fig, [root_dir() '../paper2plots/LFP/phase' int2str(neuron_num)], '-dpdf', '-r0');
