% This script loads and analyzes data relevant to the paper entitled
% "Terahertz electro-optic modulation of single-photon level pulses"
% published in Physical Review Letters (2025)
% Details on data acquisition are available in the published paper and
% commented throughout this script

% Author: Nicolas Couture, 2025

%% clean-up
clearvars,close all,clc
set(groot, 'DefaultFigureWindowStyle','normal',...
    'DefaultLineLineWidth', 1.5, ...
    'DefaultLineMarkerSize',10,...
    'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultAxesFontSize', 12, ...
    'DefaultAxesBox', 'on', ...
    'DefaultAxesTickLength',[0.02, 0.002],...
    'DefaultAxesLineWidth',1,...
    'DefaultAxesFontweight','bold');
%% Characterization of THz electric field - electro-optic sampling
root = dir();
EOS = load(fullfile('.','Data','EOS_BNA_gen_ZnTe_det.mat')); % t = time in ps, A = amplitude in V

[~, p] = max(EOS.A);
offset = 113; % [ps]

t = -(EOS.t - offset);
amp_norm = EOS.A./max(EOS.A) - mean(EOS.A./max(EOS.A)) - 0.012; %remove mean and residual offset for fft

% Super-Gaussian window for time-domain data
t_min_set = 1.15;
t_max_set = 5;

W = 0.8;
G_min = exp(-(t - (t_min_set )).^6./W.^6);
[v_ p_] = max(G_min);
G_min(1:p_) = 1;

G_max = exp(- (t - t_max_set).^6./W.^6);
[v_ p_] = max(G_max);
G_max(p_:length(G_max)) = 1;

amp_windowed = amp_norm.*G_min.*G_max ;

[Ef, fr, rho, phi] = THz_fft(amp_windowed,flipud(t),1000); %flipud so that time is increasing

figure(1),clf
tiledlayout(1,2,'TileSpacing','compact','padding','compact');
nexttile
plot(t, amp_norm, 'k')
hold on
plot(t,amp_windowed, '--r')
xlabel('Time delay [ps]'),ylabel('Amplitude [norm.]'),grid on
axis([0 10 -0.4 1.05])
legend('Raw','Windowed')

nexttile
plot(fr, rho/max(rho), 'k')
xlabel('Frequency [THz]'),ylabel('Amplitude [norm.]'),grid on
axis([0 8 0 1.05])
%% Single-photon spectral shearing as function of delay
file_list = dir(fullfile('.','Data','*DelayScan.mat'));
offset = 134; % [ps]
colo = [0.4940 0.1840 0.5560; 0 0.4470 0.7410; 0.4660 0.6740 0.1880;...
    0.8500 0.3250 0.0980;0.6350 0.0780 0.1840]; % colormap

for ii = 1:length(file_list)
    R(ii) = load(fullfile('.','Data',file_list(ii).name));
    mean_coinc(:,ii) = mean(R(ii).A_Coinc); % mean over 21 iterations
    std_coinc(:,ii) =  std(R(ii).A_Coinc);
    norm_coinc(ii) = mean(mean_coinc(end-20:end,ii));
end

figure(2),clf
tiledlayout(5,1,'Padding','compact','TileSpacing','compact');
for ii = 1:length(file_list)
    nexttile
    
    if any(ii == [1 5])
        plot(-(R(ii).t-offset),mean_coinc(:,ii),'color',[colo(ii,:),0.3])
        hold on
        plot(-(R(ii).t-offset),smooth(mean_coinc(:,ii),0.04),'color',colo(ii,:))
        if ii == 1
            title('Single photon')
        end
    else
        plot(-(R(ii).t-offset),mean_coinc(:,ii),'color',colo(ii,:))
    end
    if ii ~= length(file_list)
        set(gca,'xticklabel','');
    end
    if ii == 3
        ylabel('Counts [Hz]');
    end
    grid on
end
xlabel('Time delay [ps]'),set(gcf,'Position',[199 162 560 556])
%% Classical spectral shearing as function of delay and power

file_list = dir(fullfile('.','Data','*purged*'));
% RS = load('BNADelay_90mW_purged.mat');
offset = 117.45; % [ps]
pow_norm = 0.1:0.1:0.9;

w_fit = linspace(750,850,2E3); % wavelength axis for fitting
range = 1100:1500; % limit range of spectrometer data for fitting

fprintf('Fitting spectra...\n')
% fitting of spectra
for i_pow = 1:length(file_list)
    
    RS = load(fullfile('.','Data',file_list(i_pow).name));
    s = RS.A; % spectral amplitude, normalized
    t = -(RS.t - offset); % time delay, in ps
    w = RS.w(:,1);  % wavelength axis of spectrometer, always the same, just need one
    
    for ix=1:length(t)
        %     fit spectrum to single/triple Gaussian
        [fitted_curve,gof] = fit(w(range,1),smooth(s(range,ix))./max(RS.A(:,ix)),'gauss1');
        [fitted_curve3,gof3] = fit(w(range,1),smooth(s(range,ix))./max(RS.A(:,ix)),'gauss3');
        
        coeffvals1(:,ix) = coeffvalues(fitted_curve);
        
        %    fit to longer vector
        s_fit(:,ix) = fitted_curve3(w_fit);
        %     find indices in the most intense 10% of spectrum
        ii = s_fit(:,ix)/max(s_fit(:,ix)) > 0.9;
        w_ = w_fit(ii);
        %     find mean wavelength within most intense portion of spectrum
        wc(ix,i_pow) = mean(w_);
        fwhm(ix,i_pow) = 2*coeffvals1(3,ix); % fwhm from single gaussian fit
        clear ii;
        clear w_;
    end
    fprintf('Progress: %i/%i\n',i_pow,length(pow_norm))
end
fprintf('DONE\n')

wc_delta = wc - wc(end,:); % find delta central wavelength with/without THz
%%
figure(3),clf
subplot(2,4,[1 2 5 6])
surf(t, w, s./max(s))
colormap(jet), view(2), shading flat
axis([t(end) t(1) 793 807])
cc=colorbar(); cc.Label.String = 'Intensity (norm.)';
set(gca,'tickdir','out')
xlabel('Time delay (ps)');ylabel('Wavelength (nm)')

c = jet(length(pow_norm))*0.8;

figure(3),subplot(2,4,3)
for ii = 1:length(pow_norm)
    plot(t, wc_delta(:,ii),'color',c(ii,:))
    hold on
end
xlabel('Time delay [ps]'),ylabel('\lambda_{shift} [nm]')
leg = legend(num2str(pow_norm'));
leg.Title.String = 'E_{THz} [norm.]';
leg.Location = 'northwest';
xlim([-6 2]),grid on

figure(3),subplot(2,4,7)
for ii = 1:length(pow_norm)
    plot(t, fwhm(:,ii),'color',c(ii,:))
    hold on
end
xlabel('Time delay [ps]'),ylabel('FWHM [nm]')
leg = legend(num2str(pow_norm'));
leg.Title.String = 'E_{THz} [norm.]';
leg.Location = 'northwest';
xlim([-6 2]),grid on

% indices of max blue/red shifts
p_redshift = 68;
p_blueshift = 61;
% indices of max compression/expansion
p_compression = 74;
p_expansion = 65;

figure(3),subplot(2,4,4)
plot(pow_norm,wc_delta(p_redshift,:),'r')
hold on
plot(pow_norm,wc_delta(p_blueshift,:),'b')
xlabel('E_{THz} [norm.]')
grid on
legend('Redshift','Blueshift','Location','east')
set(gca,'yticklabel','')

figure(3),subplot(2,4,8)
plot(pow_norm,fwhm(p_compression,:))
hold on
plot(pow_norm,fwhm(p_expansion,:))
xlabel('E_{THz} [norm.]'),grid on
legend('Compression','Expansion','Location','east')
set(gca,'yticklabel','')
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

% indices for spectral shearing
idx = [1278 1287 1303 1318 1328];

% note that there is a horizontal shift relative to the single-photon
% experiment
%%
figure(4),clf
tiledlayout(5,1,'TileSpacing','compact','Padding','compact');
for ii = 1:length(idx)
    nexttile
    plot(t, sum(s(idx(ii)-5:idx(ii)+5,:)/100)+30,'color',colo(ii,:))
    grid on
    if ii ~= length(idx)
        set(gca,'xticklabel','');
    end
    if ii == 3
        ylabel('Intensity [arb. u.]');
    end
    if ii == 1
        title('Classical')
    end
end
xlabel('Time delay [ps]')
set(gcf,'Position',[759 162 560 556])
%% single-photon shear and noise count power dependence
file_list = dir(fullfile('.','Data','*Delay133p4*'));
width = 10; % bins, bin_width = 50 ps
rep_rate = 100; % Hz
pow_norm = [0.2 0.4 0.5 0.6 0.8 1];

for ii = 1:length(file_list)
    counts_mean(ii,:) = mean(load(fullfile(file_list(ii).folder,file_list(ii).name)).finaldata);
end

[~,p] = max(counts_mean(1,:)); % index at which to center coincidence window

% sum coincidences within 1 ns window
for ii=1:length(file_list)
    coinc(ii) = sum(counts_mean(ii,p-width:p+width)); % Hz
end

figure(5),clf
plot(pow_norm,1E3*coinc/rep_rate,'.-k') % divide by rep rate for counts/pulse
xlabel('E_{THz} [norm.]'),ylabel('Counts/pulse (\times10^{-3})')
file_list = dir(fullfile('.','Data','*PumpOnly*'));
pow_noise = [0.2 0.3 0.4 0.5 0.6 0.7 0.9 1];

for ii = 1:length(file_list)
    noise_mean(ii,:) = mean(load(fullfile(file_list(ii).folder,file_list(ii).name)).finaldata);
    noise_coinc(ii) = sum(noise_mean(ii,p-width:p+width));
end

figure(5)
yyaxis right
plot(pow_noise,1E3*noise_coinc/rep_rate,'.-')
ylabel('Counts/pulse (\times10^{-3})'),ylim([0 1.2])