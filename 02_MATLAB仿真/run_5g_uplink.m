function outDir = run_5g_uplink(mode)
% RUN_5G_UPLINK Reproducible CP-OFDM / localized DFT-s-OFDM study.
% Usage: run_5g_uplink('quick'), run_5g_uplink('full'), run_5g_uplink('test').
% Base MATLAB only. Intended for MATLAB R2016b or later.
% Simplified SISO 5G-uplink-inspired link, not a complete NR implementation.
% Eb/N0 includes CP energy and is measured BEFORE the PA. No FEC or pilots.
% All transforms are unitary. A memoryless hard limiter is a PA surrogate.
% See README_CN.txt and experiment_design.txt before interpreting results.

if nargin < 1, mode = 'quick'; end
assert(any(strcmp(mode, {'quick','full','test'})), 'Use quick, full, or test.');
cfg.Nfft = 256; cfg.K = 72; cfg.Ncp = 18; cfg.L = 4;
cfg.deltaF = 15e3; cfg.seed = 20261027; cfg.batch = 128;
cfg.delays = [0 1 3 5]; cfg.pdpDB = [0 -3 -6 -9];
cfg.modOrders = [4 16]; cfg.awgnDB = 0:2:16; cfg.fadingDB = 0:3:24;
cfg.iboDB = 0:1:10; cfg.paEbN0DB = 14;
cfg.paprBlocks = 3000; cfg.maxBits = 150000; cfg.minBlocks = 256;
cfg.targetErrors = 200; cfg.paBlocks = 1000; cfg.spectrumBlocks = 100;
if strcmp(mode, 'full')
    cfg.paprBlocks = 30000; cfg.maxBits = 2000000; cfg.minBlocks = 2048;
    cfg.targetErrors = 1000; cfg.paBlocks = 8000; cfg.spectrumBlocks = 600;
end
cfg.N = cfg.Nfft*cfg.L; cfg.cp = cfg.Ncp*cfg.L;
cfg.idx = (cfg.N/2-cfg.K/2+1):(cfg.N/2+cfg.K/2);
cfg.bin = (-cfg.K/2:cfg.K/2-1).';
cfg.overhead = (cfg.Nfft+cfg.Ncp)/cfg.Nfft;
cfg.mode = mode; cfg.matlabVersion = version;
rng(cfg.seed, 'twister');
selftest(cfg);
if strcmp(mode, 'test'), outDir = ''; return; end

root = fileparts(mfilename('fullpath'));
outDir = fullfile(root, 'results', [mode '_' datestr(now,'yyyymmdd_HHMMSS')]);
if exist(outDir,'dir'), outDir = [outDir '_new']; end
mkdir(outDir);
diary(fullfile(outDir, 'run_log.txt')); cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('Mode %s, seed %d, MATLAB %s\n',mode,cfg.seed,version);
fprintf('Output %s\n', outDir);
set(groot,'defaultAxesFontName','Times New Roman','defaultAxesFontSize',11);
set(groot,'defaultTextFontName','Times New Roman');
colors = [0.00 0.32 0.58; 0.85 0.30 0.10; 0.12 0.55 0.36; 0.55 0.28 0.67];
labels = {'CP-OFDM QPSK','DFT-s-OFDM QPSK','CP-OFDM 16-QAM','DFT-s-OFDM 16-QAM'};

%% 1. Oversampled PAPR statistics. CP excluded from measurement interval.
fprintf('1/6 PAPR statistics\n');
papr = zeros(cfg.paprBlocks,4);
for m = 1:2
    for first = 1:cfg.batch:cfg.paprBlocks
        n = min(cfg.batch,cfg.paprBlocks-first+1);
        [~,d] = source(cfg.modOrders(m),cfg.K,n);
        for w = 1:2
            x = transmit(d,w,cfg);
            v = 10*log10(max(abs(x).^2,[],1)./mean(abs(x).^2,1));
            papr(first:first+n-1,2*(m-1)+w) = v.';
        end
    end
end
threshold = (0:0.1:15).'; ccdf = zeros(length(threshold),4);
fig = newfig; hold on;
for c = 1:4
    ccdf(:,c) = mean(bsxfun(@gt,papr(:,c).',threshold),2);
    semilogy(threshold,positive(ccdf(:,c)),'LineWidth',1.6,'Color',colors(c,:));
end
set(gca,'YScale','log'); grid(gca,'on'); ylim([1/cfg.paprBlocks 1]);
xlabel('PAPR threshold (dB)'); ylabel('Pr(PAPR > threshold)');
title(sprintf('PAPR CCDF, %d blocks, %dx oversampling',cfg.paprBlocks,cfg.L));
legend(labels,'Location','southwest');
savefigs(fig,outDir,'fig01_papr_ccdf');
writecsv(fullfile(outDir,'papr_ccdf.csv'),{'threshold_dB','CP_QPSK','DFTs_QPSK','CP_16QAM','DFTs_16QAM'},[threshold ccdf]);
pq = zeros(4,2);
for c=1:4
    v=sort(papr(:,c)); pq(c,:)=[v(ceil(0.99*numel(v))) v(ceil(0.999*numel(v)))];
end
writecsv(fullfile(outDir,'papr_quantiles.csv'),{'curve_id','q99_dB','q999_dB'},[(1:4).' pq]);

%% 2 and 3. Linear AWGN and block Rayleigh channels.
% FFT-domain simulation is exactly equivalent to CP removal / FFT under
% a CP longer than the channel and perfect timing. selftest verifies this.
berResults = struct;
for channel = {'AWGN','Rayleigh'}
    cname = channel{1};
    if strcmp(cname,'AWGN'), axisDB=cfg.awgnDB; num=2; else, axisDB=cfg.fadingDB; num=3; end
    fprintf('%d/6 BER in %s\n',num,cname);
    matrix=zeros(numel(axisDB),4); rows=[];
    for m=1:2
        M=cfg.modOrders(m);
        for s=1:numel(axisDB)
            st=linearPoint(M,axisDB(s),cname,cfg);
            matrix(s,2*m-1:2*m)=st.ber;
            for w=1:2
                rows(end+1,:)=[M w axisDB(s) st.errors(w) st.bits st.blocks st.ber(w) st.blockSE(w)]; %#ok<AGROW>
            end
            fprintf('  M=%d, Eb/N0=%g dB, BER=[%.3g %.3g], bits=%d\n',M,axisDB(s),st.ber,st.bits);
        end
    end
    fig=newfig; hold on;
    for c=1:4
        semilogy(axisDB,positive(matrix(:,c)),'-o','LineWidth',1.4,'Color',colors(c,:),'MarkerSize',4);
    end
    leg=labels;
    if strcmp(cname,'AWGN')
        g=10.^(axisDB/10)/cfg.overhead;
        qpsk=0.5*erfc(sqrt(g));
        a=sqrt(0.8*g); Q=@(x) 0.5*erfc(x/sqrt(2));
        qam16=0.75*Q(a)+0.5*Q(3*a)-0.25*Q(5*a);
        semilogy(axisDB,positive(qpsk),'k--','LineWidth',1.2);
        semilogy(axisDB,positive(qam16),'k:','LineWidth',1.2);
        leg=[labels {'QPSK theory with CP cost','16-QAM theory with CP cost'}];
    end
    set(gca,'YScale','log'); grid(gca,'on'); ylim([1/cfg.maxBits 0.5]);
    xlabel('Pre-PA E_b/N_0 including CP energy (dB)'); ylabel('Uncoded bit error rate');
    title([cname ' channel, perfect CSI, ZF equalization']);
    legend(leg,'Location','southwest');
    savefigs(fig,outDir,sprintf('fig%02d_ber_%s',num,lower(cname)));
    writecsv(fullfile(outDir,['ber_' lower(cname) '.csv']), ...
        {'M','waveform_id','EbN0_dB','error_bits','total_bits','blocks','BER','block_standard_error'},rows);
    berResults.(cname).curves=matrix; berResults.(cname).raw=rows;
end

%% 4 and 5. PA limiter: BER and noiseless in-band EVM, 16-QAM.
% Both waveforms have the same ensemble mean input power and threshold.
% There is NO post-PA power renormalization. Gain estimated on a separate
% independent calibration sequence is compensated at the receiver.
fprintf('4/6 and 5/6 PA back-off experiment\n');
paRows=[]; paBER=zeros(numel(cfg.iboDB),2); paEVM=paBER;
for s=1:numel(cfg.iboDB)
    ibo=cfg.iboDB(s); A=sqrt(cfg.K/cfg.N)*10^(ibo/20);
    gain=ones(1,2); pin=zeros(1,2); pout=pin;
    [~,dcal]=source(16,cfg.K,512);
    for w=1:2
        x=transmit(dcal,w,cfg); z=limiter(x,A);
        gain(w)=real(sum(conj(x(:)).*z(:))/sum(abs(x(:)).^2));
    end
    errors=zeros(1,2); evmNumer=errors; evmDenom=0;
    for first=1:cfg.batch:cfg.paBlocks
        n=min(cfg.batch,cfg.paBlocks-first+1);
        [b,d]=source(16,cfg.K,n);
        nv=cfg.overhead/(4*10^(cfg.paEbN0DB/10));
        noise=sqrt(nv/2)*(randn(cfg.K,n)+1i*randn(cfg.K,n));
        evmDenom=evmDenom+sum(abs(d(:)).^2);
        for w=1:2
            x=transmit(d,w,cfg); z=limiter(x,A);
            rxGrid=fftshift(fft(z,[],1)/sqrt(cfg.N),1); v=rxGrid(cfg.idx,:);
            dhat=despread((v+noise)/gain(w),w);
            dec=qamDemap(dhat,16); errors(w)=errors(w)+sum(dec(:)~=b(:));
            clean=despread(v/gain(w),w);
            evmNumer(w)=evmNumer(w)+sum(abs(clean(:)-d(:)).^2);
            pin(w)=pin(w)+sum(abs(x(:)).^2); pout(w)=pout(w)+sum(abs(z(:)).^2);
        end
    end
    bits=cfg.paBlocks*cfg.K*4;
    paBER(s,:)=errors/bits; paEVM(s,:)=100*sqrt(evmNumer/evmDenom);
    for w=1:2
        paRows(end+1,:)=[ibo w cfg.paEbN0DB errors(w) bits paBER(s,w) paEVM(s,w) gain(w) 10*log10(pout(w)/pin(w))]; %#ok<AGROW>
    end
    fprintf('  IBO=%g dB, BER=[%.3g %.3g], EVM=[%.2f %.2f] %%\n',ibo,paBER(s,:),paEVM(s,:));
end
fig=newfig;
semilogy(cfg.iboDB,positive(paBER(:,1)),'-o',cfg.iboDB,positive(paBER(:,2)),'-s','LineWidth',1.5);
grid(gca,'on'); xlabel('PA input back-off (dB)'); ylabel('Uncoded bit error rate');
title(sprintf('16-QAM with hard limiting, pre-PA E_b/N_0 = %g dB',cfg.paEbN0DB));
legend({'CP-OFDM','DFT-s-OFDM'},'Location','northeast');
savefigs(fig,outDir,'fig04_ber_pa_backoff');
fig=newfig;
plot(cfg.iboDB,paEVM(:,1),'-o',cfg.iboDB,paEVM(:,2),'-s','LineWidth',1.5);
grid(gca,'on'); xlabel('PA input back-off (dB)'); ylabel('RMS in-band EVM after gain compensation (%)');
title('16-QAM distortion without receiver noise');
legend({'CP-OFDM','DFT-s-OFDM'},'Location','northeast');
savefigs(fig,outDir,'fig05_evm_pa_backoff');
writecsv(fullfile(outDir,'pa_backoff.csv'), ...
    {'IBO_dB','waveform_id','EbN0_dB','error_bits','total_bits','BER','EVM_percent','gain','output_input_power_dB'},paRows);

%% 6. Welch PSD implemented here, including CP and symbol boundaries.
fprintf('6/6 Power spectral density\n');
[~,d]=source(16,cfg.K,cfg.spectrumBlocks);
psdAll=[]; fs=cfg.N*cfg.deltaF;
for w=1:2
    x=transmit(d,w,cfg); frame=[x(end-cfg.cp+1:end,:);x]; stream=frame(:);
    for state=1:2
        sig=stream;
        if state==2, sig=limiter(sig,sqrt(cfg.K/cfg.N)*10^(2/20)); end
        [f,p]=welchLocal(sig,fs,4*cfg.N);
        psdAll(:,2*(w-1)+state)=p; %#ok<AGROW>
    end
end
referenceLevel=max(psdAll(:,1)); psdDB=10*log10(max(psdAll/referenceLevel,1e-14));
fig=newfig; plot(f/1e6,psdDB,'LineWidth',1.3); grid(gca,'on');
xlim([-3 3]); ylim([-80 5]); xlabel('Baseband frequency (MHz)');
ylabel('PSD relative to linear CP-OFDM peak (dB)');
title('16-QAM spectrum, rectangular pulses, common normalization');
legend({'CP-OFDM linear','CP-OFDM IBO 2 dB','DFT-s-OFDM linear','DFT-s-OFDM IBO 2 dB'},'Location','southwest');
savefigs(fig,outDir,'fig06_psd');
writecsv(fullfile(outDir,'psd.csv'),{'frequency_Hz','CP_linear_dB','CP_clipped_dB','DFTs_linear_dB','DFTs_clipped_dB'},[f psdDB]);

save(fullfile(outDir,'results.mat'),'cfg','papr','threshold','ccdf','pq','berResults','paRows','paBER','paEVM','f','psdAll');
fid=fopen(fullfile(outDir,'summary.txt'),'w');
fprintf(fid,'Mode: %s\nSeed: %d\nMATLAB: %s\n',mode,cfg.seed,version);
fprintf(fid,'Curve order: CP QPSK, DFTs QPSK, CP 16QAM, DFTs 16QAM.\n');
fprintf(fid,'PAPR q99 / q999 in dB:\n');
for c=1:4, fprintf(fid,'%s: %.4f / %.4f\n',labels{c},pq(c,:)); end
fprintf(fid,'No zero-error point has been replaced by a fabricated positive BER.\n');
fprintf(fid,'Zero-error BER/CCDF points are omitted on logarithmic plots; raw counts are in CSV.\n');
fprintf(fid,'Rayleigh uncertainty: block-level standard error, not an iid-bit interval.\n');
fprintf(fid,'PA results use pre-PA energy. Output power loss is saved; no battery-efficiency claim.\n');
fprintf(fid,'Full mode is recommended for the report. Quick mode is for debugging.\n');
fclose(fid);
fprintf('Complete. Return this entire results folder for report and slides.\n');
fprintf('Keep raw CSV and MAT files, not just screenshots.\n');
end

function [b,d]=source(M,K,n)
b=randi([0 1],log2(M),K*n); d=reshape(qamMap(b,M),K,n);
end

function d=qamMap(b,M)
p=log2(M)/2;
if M==4, levels=[-1 1]; elseif M==16, levels=[-3 -1 3 1]; else, error('Only QPSK and 16-QAM'); end
weights=2.^(p-1:-1:0);
ii=weights*b(1:p,:); qq=weights*b(p+1:end,:);
d=(levels(ii+1)+1i*levels(qq+1))/sqrt(2*(M-1)/3);
end

function b=qamDemap(d,M)
p=log2(M)/2;
if M==4, levels=[-1 1]; else, levels=[-3 -1 3 1]; end
v=d(:)*sqrt(2*(M-1)/3);
[~,ii]=min(abs(bsxfun(@minus,real(v),levels)),[],2);
[~,qq]=min(abs(bsxfun(@minus,imag(v),levels)),[],2);
b=zeros(2*p,numel(v));
for j=1:p
    b(j,:)=bitget(uint32(ii-1),p-j+1).';
    b(p+j,:)=bitget(uint32(qq-1),p-j+1).';
end
end

function z=spread(d,w)
if w==2, z=fft(d,[],1)/sqrt(size(d,1)); else, z=d; end
end

function d=despread(z,w)
if w==2, d=ifft(z,[],1)*sqrt(size(z,1)); else, d=z; end
end

function x=transmit(d,w,cfg)
txGrid=zeros(cfg.N,size(d,2)); txGrid(cfg.idx,:)=spread(d,w);
x=ifft(ifftshift(txGrid,1),[],1)*sqrt(cfg.N);
end

function z=limiter(x,A)
z=x.*min(1,A./max(abs(x),eps));
end

function st=linearPoint(M,ebDB,cname,cfg)
bitsPerBlock=cfg.K*log2(M);
maxBlocks=max(cfg.minBlocks,floor(cfg.maxBits/bitsPerBlock));
e=zeros(1,2); sums2=e; blocks=0;
pdp=10.^(cfg.pdpDB/10); pdp=pdp/sum(pdp);
response=exp(-1i*2*pi*cfg.bin*cfg.delays/cfg.Nfft);
nv=cfg.overhead/(log2(M)*10^(ebDB/10));
while blocks<maxBlocks && (blocks<cfg.minBlocks || any(e<cfg.targetErrors))
    n=min(cfg.batch,maxBlocks-blocks); [b,d]=source(M,cfg.K,n);
    if strcmp(cname,'Rayleigh')
        taps=bsxfun(@times,sqrt(pdp(:)/2),randn(numel(pdp),n)+1i*randn(numel(pdp),n));
        H=response*taps; % E[|H|^2]=1; do NOT normalize each realization.
    else
        H=ones(cfg.K,n);
    end
    noise=sqrt(nv/2)*(randn(cfg.K,n)+1i*randn(cfg.K,n));
    for w=1:2
        z=H.*spread(d,w)+noise;
        dhat=despread(z./H,w); % exact ZF; deep fades legitimately amplify noise.
        dec=qamDemap(dhat,M); err=reshape(dec~=b,bitsPerBlock,n);
        eb=sum(err,1); e(w)=e(w)+sum(eb);
        sums2(w)=sums2(w)+sum((eb/bitsPerBlock).^2);
    end
    blocks=blocks+n;
end
st.bits=blocks*bitsPerBlock; st.blocks=blocks; st.errors=e; st.ber=e/st.bits;
samplevar=max(0,(sums2-blocks*st.ber.^2)/(blocks-1));
st.blockSE=sqrt(samplevar/blocks);
end

function selftest(cfg)
fprintf('Running deterministic modem and CP tests...\n');
for M=[4 16]
    b=zeros(log2(M),M);
    for j=1:log2(M), b(j,:)=bitget(uint32(0:M-1),log2(M)-j+1); end
    assert(isequal(b,qamDemap(qamMap(b,M),M)),'QAM round trip failed.');
    assert(abs(mean(abs(qamMap(b,M)).^2)-1)<1e-12,'Constellation energy failed.');
end
[b,d]=source(16,cfg.K,3);
h=zeros(1,max(cfg.delays)*cfg.L+1); h(cfg.delays*cfg.L+1)=[1 0.3+0.2i -0.1i 0.1];
H=fftshift(fft(h,cfg.N)).';
for w=1:2
    x=transmit(d,w,cfg);
    assert(max(abs(sum(abs(x).^2,1)-sum(abs(d).^2,1)))<1e-9,'Energy preservation failed.');
    xc=[x(end-cfg.cp+1:end,:);x]; y=filter(h,1,xc,[],1);
    g=fftshift(fft(y(cfg.cp+(1:cfg.N),:),[],1)/sqrt(cfg.N),1);
    dhat=despread(bsxfun(@rdivide,g(cfg.idx,:),H(cfg.idx)),w);
    assert(max(abs(dhat(:)-d(:)))<1e-9,'CP convolution / FFT equivalence failed.');
    assert(isequal(b,qamDemap(dhat,16)),'Noiseless link failed.');
    c1=cfg; c1.N=cfg.Nfft; c1.idx=(c1.N/2-c1.K/2+1):(c1.N/2+c1.K/2);
    xb=transmit(d,w,c1); xs=x(1:cfg.L:end,:)*sqrt(cfg.L);
    assert(max(abs(xs(:)-xb(:)))<1e-9,'Oversampling consistency failed.');
end
one=ones(cfg.K,1); a=transmit(one,1,cfg); s=transmit(one,2,cfg);
assert(abs(max(abs(a).^2)/mean(abs(a).^2)-cfg.K)<1e-9,'OFDM peak test failed.');
assert(abs(max(abs(s).^2)/mean(abs(s).^2)-1)<1e-9,'DFT tone test failed.');
assert(cfg.Ncp>max(cfg.delays),'CP does not cover channel.');
fprintf('PASS: Gray mapping, unit energy, CP channel, zero BER, oversampling, known PAPR.\n');
end

function [f,p]=welchLocal(x,fs,nfft)
len=nfft; hop=floor(len/2); window=0.5-0.5*cos(2*pi*(0:len-1)'/(len-1));
p=zeros(nfft,1); count=0;
for a=1:hop:(numel(x)-len+1)
    s=x(a:a+len-1).*window;
    p=p+abs(fft(s,nfft)).^2/(fs*sum(window.^2)); count=count+1;
end
assert(count>0,'Not enough samples for PSD.');
p=fftshift(p/count); f=(-nfft/2:nfft/2-1)'*fs/nfft;
end

function y=positive(x)
y=x; y(y<=0)=NaN;
end

function fig=newfig
fig=figure('Color','w','Position',[100 100 900 600]);
end

function savefigs(fig,path,name)
savefig(fig,fullfile(path,[name '.fig']));
print(fig,fullfile(path,[name '.png']),'-dpng','-r300');
set(fig,'PaperPositionMode','auto');
print(fig,fullfile(path,[name '.pdf']),'-dpdf','-bestfit');
end

function writecsv(path,header,data)
fid=fopen(path,'w'); assert(fid>=0,'Cannot create CSV.');
fprintf(fid,'%s',header{1});
for k=2:numel(header), fprintf(fid,',%s',header{k}); end
fprintf(fid,'\n'); fclose(fid);
dlmwrite(path,data,'-append','delimiter',',','precision','%.12g');
end
