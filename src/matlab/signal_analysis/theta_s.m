%% ------------------ LOAD LFP DATA ------------------
% Loads allfp_ds for this session if it is not already in the workspace.
% This is the raw/downsampled LFP matrix used later for channel selection.
% Expected shape for this file is channels x samples.

scriptDir = fileparts(mfilename('fullpath'));
addpath(genpath(scriptDir));
compatDir = fullfile(scriptDir, 'compat');
if exist(compatDir, 'dir')
    addpath(compatDir, '-begin');
end

if ~exist('allfp_ds', 'var')
    repoRoot = 'C:\Users\tadse\OneDrive\Documenti\GitHub\CA1-cellclass';
    lfpPath = fullfile(repoRoot, 'data', 'raw', 'lfp', ...
        'VS103_2024-11-10_14-58-13_lfp1250Hz.mat');

    assert(exist(lfpPath, 'file') == 2, 'LFP file not found: %s', lfpPath);
    load(lfpPath, 'allfp_ds', 'channels');
    fprintf('Loaded LFP data from: %s\n', lfpPath);
end

%% ------------------ LOAD TRAJECTORY DATA ------------------
% Loads Traj and behavior variables for the same session.
% theta_s.m needs trial metadata (Traj), speed (XSpeed), position (X_ds_n),
% and the behavior time vector (t_ds). If the MAT file only contains Traj,
% behavior_vectors_from_traj reconstructs those vectors from trial fields.

if ~exist('Traj', 'var') || ~exist('XSpeed', 'var') || ...
        ~exist('X_ds_n', 'var') || ~exist('t_ds', 'var')
    trajDataPath = ['\\Epsztein-nas02\TEAM\Tadeo\data\trajdata\' ...
        'VS103_2024-11-10_14-58-13_TrajData.mat'];

    assert(exist(trajDataPath, 'file') == 2, ...
        'Trajectory data file not found: %s', trajDataPath);

    trajData = load(trajDataPath);

    if isfield(trajData, 'Traj')
        Traj = trajData.Traj;
    end
    if isfield(trajData, 'XSpeed')
        XSpeed = trajData.XSpeed;
    end
    if isfield(trajData, 'X_ds_n')
        X_ds_n = trajData.X_ds_n;
    end
    if isfield(trajData, 't_ds')
        t_ds = trajData.t_ds;
    end

    assert(exist('Traj', 'var') == 1, ...
        'The trajectory file did not contain a Traj variable.');

    if ~exist('XSpeed', 'var') || ~exist('X_ds_n', 'var') || ...
            ~exist('t_ds', 'var')
        [t_ds, X_ds_n, XSpeed] = behavior_vectors_from_traj(Traj);
        fprintf(['Reconstructed t_ds, X_ds_n, and XSpeed from Traj ' ...
            'trial fields.\n']);
    end

    fprintf('Loaded trajectory data from: %s\n', trajDataPath);
end

%% ------------------ PARAMETRES ------------------
% Define analysis settings and copy loaded variables into the names used by
% the legacy theta code. This cell creates prm and the theta behavior fields:
% theta.veld = speed, theta.posd = position, theta.timed = behavior time.
% Chargement des donn�es
lfp = allfp_ds;       % [nChannels x nSamples] LFP downsampled
prm.fq_lfp = 1250;    % fr�quence LFP (Hz)
prm.fd = 1000;        % fr�quence d'analyse finale (Hz)
prm.theta_band = [6 9]; % Bande de fr�quence th�ta (Hz)
prm.delta_band = [0.5 4]; % Bande de fr�quence delta (Hz)
prm.vel_th = 2;       % Seuil de vitesse (cm/s)
prm.frange = 0:0.5:20; % Plage de fr�quences pour l'analyse (Hz)
verbose = true;      % Affichage des informations d�taill�es

% Fen�tres temporelles (en secondes)
prm.win1 = [Traj(find([Traj.Cond]==1,1,'first')).tstart Traj(find([Traj.Cond]==1,1,'last')).tstop];
prm.win2 = [Traj(find([Traj.Cond]==2,1,'first')).tstart Traj(find([Traj.Cond]==2,1,'last')).tstop];
% V�rification des fen�tres temporelles
disp('V�rification des fen�tres temporelles (en secondes) :');
disp(['Fen�tre 1 : ' num2str(prm.win1(1)) ' � ' num2str(prm.win1(2))]);
disp(['Fen�tre 2 : ' num2str(prm.win2(1)) ' � ' num2str(prm.win2(2))]);

% Comportement
[theta.veld] = XSpeed; % Vitesse
[theta.posd] = X_ds_n; % Position
theta.timed = t_ds;    % Temps (en secondes)

% V�rification des plages temporelles
disp('V�rification des plages temporelles :');
disp(['Plage de theta.timed : ' num2str(min(theta.timed)) ' � ' num2str(max(theta.timed)) ' secondes']);

% Taille des donn�es
[S, N] = size(lfp);  % S = nombre de canaux, N = nombre d'�chantillons
Nd = length(fct_dsampling(lfp(1,:), prm.fq_lfp, prm.fd));

%% ------------------ ALIGNEMENT LFP / COMPORTEMENT ------------------
% Put LFP and behavior on one shared sample grid so later indexing uses the
% same time axis. This legacy version interpolates allfp_ds to the behavior
% length. After interpolation the aligned LFP is treated as prm.fd Hz.
% V�rification de la diff�rence de taille entre LFP et comportement
disp(['Diff�rence de points entre LFP et comportement : ' num2str(size(allfp_ds,2) - size(X_ds_n,1))]);

% Nombre de points
nLFP = size(allfp_ds,2);
nBHV = size(X_ds_n,1);

% Vecteurs temps normalis�s
tLFP = linspace(0,1,nLFP);
tBHV = linspace(0,1,nBHV);

% Alignement du comportement sur le LFP
lfp = interp1(tLFP, allfp_ds', tBHV, 'linear')';
prm.fq_lfp = prm.fd;
fprintf(['LFP was interpolated to behavior length (%d samples); ' ...
    'using prm.fq_lfp = prm.fd = %g Hz from here.\n'], nBHV, prm.fd);

% S�lection de la fen�tre d'int�r�t
idx = ceil(prm.win1(1)):ceil(prm.win1(2));

%% ------------------ SELECTION DU MEILLEUR CANAL ------------------
% Score every LFP channel by theta/delta amplitude during condition 1.
% For each channel the script filters in theta and delta, computes Hilbert
% envelopes, then selects the channel with the largest median theta/delta
% envelope ratio. The selected channel is stored in prm.igoodch and v.
disp('S�lection du meilleur canal (ratio th�ta/delta)');
mtd = 1; % M�thode Hilbert (plus rapide que les ondelettes)
thetapower = nan(S,1);
deltapower = nan(S,1);

% M�thode Hilbert pour le calcul du ratio th�ta/delta
switch mtd
    case 0 % M�thode ondelettes
        for k = S:-1:1
            lfpd_tmp = fct_dsampling(lfp(k, idx), prm.fq_lfp, prm.fd);
            [cwt_tmp, ffreq] = fct_cwt(lfpd_tmp, 'fs', prm.fd, 'frange', prm.frange);
            cwt_tmp = abs(cwt_tmp);

            idx_theta = ffreq >= prm.theta_band(1) & ffreq <= prm.theta_band(2);
            thetapower(k) = nanmedian(nanmedian(cwt_tmp(idx_theta, :), 1));

            idx_delta = ffreq >= prm.delta_band(1) & ffreq <= prm.delta_band(2);
            deltapower(k) = nanmedian(nanmedian(cwt_tmp(idx_delta, :), 1));
        end
    case 1 % M�thode Butterworth + Hilbert
        [B, A] = butter(2, [prm.theta_band(1)/(prm.fq_lfp/2) prm.theta_band(2)/(prm.fq_lfp/2)], 'bandpass');
        [B2, A2] = butter(2, [prm.delta_band(1)/(prm.fq_lfp/2) prm.delta_band(2)/(prm.fq_lfp/2)], 'bandpass');
        for k = S:-1:1
            lfp_theta = filtfilt(B, A, lfp(k, idx));
            lfp_delta = filtfilt(B2, A2, lfp(k, idx));
            [~, ~, theta_mod] = fct_hilbert(lfp_theta, prm.fq_lfp);
            thetapower(k) = nanmedian(theta_mod);
            [~, ~, delta_mod] = fct_hilbert(lfp_delta, prm.fq_lfp);
            deltapower(k) = nanmedian(delta_mod);
        end
end

% Calcul du ratio th�ta/delta
thetadeltaratio = thetapower ./ deltapower;
[~, prm.igoodch] = max(thetadeltaratio);
disp(['Le canal th�ta sera analys� dans le canal : ' num2str(prm.igoodch)]);

% S�lection du meilleur canal
v = lfp(prm.igoodch, :);

%% ------------------ DETECTION DES ARTEFACTS ------------------
% Placeholder for manual artifact removal. The code is commented out, so no
% artifact correction is currently applied and v is passed forward unchanged.
% Optionnel : D�tection et suppression des artefacts
% [~, iseq_evt] = ui_event_selection(v, [], []);
% vm_rem = v;
% for k = 1:size(iseq_evt, 1)
%     xind = iseq_evt(k,1):iseq_evt(k,2);
%     vm_rem(xind) = linspace(v(iseq_evt(k,1)), v(iseq_evt(k,2)), length(xind));
% end

%% ------------------ ANALYSE TH�TA ET DELTA ------------------
disp('Analyse des amplitudes th�ta et delta');

% Filtrage et enveloppe th�ta
% Cell purpose: analyze only the selected best channel v. This cell
% bandpass-filters v in theta and delta, computes Hilbert envelopes, stores
% filtered waveforms and amplitude envelopes, then compares win1 versus win2
% after masking low-speed samples.
[B_theta, A_theta] = butter(2, [prm.theta_band(1)/(prm.fq_lfp/2) prm.theta_band(2)/(prm.fq_lfp/2)], 'bandpass');
lfp_theta = filtfilt(B_theta, A_theta, v);
[~, ~, lfp_theta_mod] = fct_thetahilb(lfp_theta, prm.fq_lfp);
theta.lfp_thetad = fct_dsampling(lfp_theta, prm.fq_lfp, prm.fd);
theta.lfp_thetad_mod = fct_dsampling(lfp_theta_mod, prm.fq_lfp, prm.fd);

% Filtrage et enveloppe delta
[B_delta, A_delta] = butter(2, [prm.delta_band(1)/(prm.fq_lfp/2) prm.delta_band(2)/(prm.fq_lfp/2)], 'bandpass');
lfp_delta = filtfilt(B_delta, A_delta, v);
[~, ~, lfp_delta_mod] = fct_thetahilb(lfp_delta, prm.fq_lfp);
theta.lfp_deltad = fct_dsampling(lfp_delta, prm.fq_lfp, prm.fd);
theta.lfp_deltad_mod = fct_dsampling(lfp_delta_mod, prm.fq_lfp, prm.fd);

% Comparaison avant/apr�s muscimol
tmp_theta = theta.lfp_thetad_mod;
tmp_theta(theta.veld < prm.vel_th) = NaN;
tmp_delta = theta.lfp_deltad_mod;
tmp_delta(theta.veld < prm.vel_th) = NaN;

idx1 = ceil(prm.win1(1)):ceil(prm.win1(2));
idx2 = ceil(prm.win2(1)):ceil(prm.win2(2));

theta.mean_thetadelta = [nanmedian(tmp_theta(idx1) - tmp_delta(idx1)) nanmedian(tmp_theta(idx2) - tmp_delta(idx2))];
theta.mean_theta = [nanmedian(tmp_theta(idx1)) nanmedian(tmp_theta(idx2))];

%% ------------------ TRANSFORM�E EN ONDELETTES ------------------
disp('Transform�e en ondelettes');
wltnorm = 'freq'; % Normalisation par fr�quence

% LFP downsampled
% Cell purpose: compute a time-frequency representation of the selected LFP
% channel. fct_cwt returns frequency x time amplitudes over prm.frange. The
% chosen normalization ('freq') multiplies amplitudes by sqrt(frequency),
% making higher frequencies more visible relative to low-frequency LFP power.
theta.lfpd = fct_dsampling(v, prm.fq_lfp, prm.fd);
[cwtd, ffreq] = fct_cwt(theta.lfpd, 'fs', prm.fd, 'frange', prm.frange);
cwtd = abs(cwtd);

% Normalisation
switch wltnorm
    case 'zscore'
        tmp_m = bsxfun(@minus, cwtd, nanmean(cwtd, 1));
        tmp_v = 1 ./ nanstd(cwtd, [], 1);
        cwtdz = bsxfun(@times, tmp_m, tmp_v);
    case 'freq'
        cwtdz = bsxfun(@rdivide, abs(cwtd), 1./sqrt(ffreq)');
    case '01'
        cwtdz = fct_normnegmat(cwtd');
        cwtdz = cwtdz';
    otherwise
        cwtdz = cwtd;
end

%% ------------------ AMPLITUDE TH�TA PENDANT LES MOUVEMENTS ------------------
prm.thresh_speed = 2; % Seuil de vitesse (cm/s)
prm.duration = 2*prm.fd; % Dur�e minimale de mouvement/immobilit� (en �chantillons)
prm.intermov_dur = 0.5*prm.fd; % Dur�e maximale entre mouvements pour fusionner (en �chantillons)

% D�tection des p�riodes de mouvement
% Cell purpose: classify time samples into movement/immobility from speed.
% Movement is speed > prm.thresh_speed, adjacent runs can be merged, and
% short runs are removed. The resulting movement intervals index theta
% amplitude and wavelet theta-band amplitude.
MI = fct_find_MI(theta.veld, prm);
theta.idxM = MI.ind_M;
nmov = MI.nb_M;

if nmov > 0
    maxMovementIndex = max(theta.idxM(:, 2));
    assert(maxMovementIndex <= numel(theta.lfp_thetad_mod), ...
        ['Movement indices reach %d, but theta.lfp_thetad_mod has only %d ' ...
        'samples. Behavior and theta traces are not on the same time base.'], ...
        maxMovementIndex, numel(theta.lfp_thetad_mod));
    assert(maxMovementIndex <= size(cwtdz, 2), ...
        ['Movement indices reach %d, but cwtdz has only %d time bins. ' ...
        'Behavior and wavelet traces are not on the same time base.'], ...
        maxMovementIndex, size(cwtdz, 2));
end

% Bande de fr�quence th�ta
idx_theta = ffreq >= prm.theta_band(1) & ffreq <= prm.theta_band(2);

% Calcul de l'amplitude th�ta pendant les mouvements
for k = nmov:-1:1
    idx = theta.idxM(k, 1):theta.idxM(k, 2);
    theta.theta_mov(k) = nanmedian(theta.lfp_thetad_mod(idx));
    theta.theta_mov2(k) = nanmedian(nanmedian(cwtdz(idx_theta, idx), 1));
end

%% ------------------ VISUALISATION ------------------
if verbose
    % Cell purpose: make diagnostic figures.
    % Figure 1 stacks wavelet amplitude, theta-filtered LFP, and position;
    % green segments mark movement periods. Figure 2 compares the LFP
    % spectrum in win1 versus win2 after excluding low-speed samples.
    visualsmooth = 10; % Lissage visuel (en secondes)

    % Figure 1 : Ondelettes, th�ta et position
    h(1) = figure;
    ax(1) = subplot(3, 2, [1 2]);
    imagesc(theta.timed, ffreq, fct_smoothgauss(cwtdz, floor(visualsmooth*prm.fd)));
    ylim([0 12]);
    colormap jet;
    ylabel('Fr�quences (Hz)', 'FontSize', 14);

    ax(2) = subplot(3, 2, [3 4]);
    plot(theta.timed, theta.lfp_thetad);
    hold on;
    for k = 1:nmov
        idx = theta.idxM(k, 1):theta.idxM(k, 2);
        plot(theta.timed(idx), theta.lfp_thetad(idx), 'g');
    end
    ylabel('Th�ta (5-10 Hz)', 'FontSize', 14);

    ax(3) = subplot(3, 2, [5 6]);
    plot(theta.timed, theta.posd);
    hold on;
    for k = 1:nmov
        idx = theta.idxM(k, 1):theta.idxM(k, 2);
        plot(theta.timed(idx), theta.posd(idx), 'g');
    end
    ylabel('Position (cm)', 'FontSize', 14);
    linkaxes(ax, 'x');
    fct_fullscreen(gcf);

    % Figure 2 : Spectre de puissance
    h(2) = figure;
    color_set = [65 97 120; 120 175 218] / 255;
    color1 = color_set(1, :);
    color2 = color_set(2, :);

    idx1 = ceil(prm.win1(1)):ceil(prm.win1(2));
    idx2 = ceil(prm.win2(1)):ceil(prm.win2(2));

    lfpd = theta.lfpd;
    lfpd(theta.veld < prm.vel_th) = NaN;
    lfpd1 = lfpd(idx1);
    lfpd2 = lfpd(idx2);
    lfpd1(isnan(lfpd1)) = [];
    lfpd2(isnan(lfpd2)) = [];

    [yf, xf] = fct_spectrum(lfpd1, 250, 'freq_range', [0 20], 'output', 'amplitude');
    yf = fct_smoothgauss(yf, 200);
    plot(xf, yf, 'Color', color1, 'LineWidth', 2);
    bfmusci = max(yf(xf > prm.theta_band(1) & xf < prm.theta_band(2)));

    hold on;
    [yf, xf] = fct_spectrum(lfpd2, 250, 'freq_range', [0 20], 'output', 'amplitude');
    yf = fct_smoothgauss(yf, 200);
    plot(xf, yf, 'Color', color2, 'LineWidth', 2);
    afmusci = max(yf(xf > prm.theta_band(1) & xf < prm.theta_band(2)));

    theta.theta_attenuation_spectrum = 100 - (afmusci*100) / bfmusci;
    title(['Att�nuation th�ta : ' num2str(theta.theta_attenuation_spectrum) ' %']);
    xlim([0 20]);
    ylabel('Amplitude', 'FontSize', 17);
    xlabel('Fr�quences (Hz)', 'FontSize', 17);
end

%% ------------------ PR�CESSION TH�TA ------------------
% Conversion des temps des spikes en secondes

spk_time = double(allcel.itime_spk) / 25000; % Si en �chantillons � 25 kHz
% spk_time = double(allcel.itime_spk) / 1000; % Si en millisecondes

% V�rification des plages temporelles
disp('V�rification des plages temporelles :');
disp(['Plage de spk_time : ' num2str(min(spk_time)) ' � ' num2str(max(spk_time))]);
disp(['Plage de theta.timed : ' num2str(min(theta.timed)) ' � ' num2str(max(theta.timed))]);

% Ajout d'un offset si n�cessaire (exemple : 2677.344 secondes)
offset = 2677.344; % � ajuster selon vos donn�es
spk_time = spk_time + offset;

% Restriction � la fen�tre win1
% Cell purpose: estimate theta phase precession for spikes. Spike times are
% converted to seconds, restricted to win1, assigned theta phase from the
% filtered LFP, and assigned position by interpolation. The final scatter is
% spike position versus theta phase.
valid_idx = (spk_time >= prm.win1(1)) & (spk_time <= prm.win1(2));
spk_time = spk_time(valid_idx);

% V�rification du nombre de spikes retenus
disp(['Nombre de spikes retenus dans win1 : ' num2str(numel(spk_time))]);

% Filtrage et phase th�ta
[B_theta, A_theta] = butter(2, [prm.theta_band(1)/(prm.fq_lfp/2) prm.theta_band(2)/(prm.fq_lfp/2)], 'bandpass');
lfp_theta = filtfilt(B_theta, A_theta, v);
analytic_signal = hilbert(lfp_theta);
theta_phase = mod(angle(analytic_signal), 2*pi);

% Fen�tre temporelle correspondante
idx_win1 = ceil(prm.win1(1)):ceil(prm.win1(2));
time_win1 = theta.timed(idx_win1);
pos_win1 = theta.posd(idx_win1);
phase_win1 = theta_phase(idx_win1);

% Interpolation des spikes -> phase et position
spk_phase = interp1(time_win1, phase_win1, spk_time, 'nearest', NaN);
spk_pos = interp1(theta.timed, theta.posd, spk_time, 'linear', NaN);

% V�rification des NaN dans spk_phase et spk_pos
disp(['Nombre de NaN dans spk_phase : ' num2str(sum(isnan(spk_phase)))]);
disp(['Nombre de NaN dans spk_pos : ' num2str(sum(isnan(spk_pos)))]);

% Visualisation de la pr�cession th�ta
figure;
scatter(spk_pos, spk_phase, 15, 'k', 'filled');
hold on;
xlabel('Position (cm)');
ylabel('Phase th�ta (rad)');
title('Pr�cession de la phase th�ta - condition 1');
ylim([0 2*pi]);

% Ajout d'une ligne de tendance circulaire-lin�aire
p = polyfit(spk_pos(~isnan(spk_pos)), spk_phase(~isnan(spk_phase)), 1);
xfit = linspace(min(spk_pos(~isnan(spk_pos))), max(spk_pos(~isnan(spk_pos))), 200);
yfit = polyval(p, xfit);
plot(xfit, mod(yfit,2*pi), 'r','LineWidth',2);

% Mesure de la force de pr�cession
[r, pval] = circ_corrcl(spk_phase(~isnan(spk_phase)), spk_pos(~isnan(spk_pos)));
disp(['Corr�lation circulaire-lin�aire : r=' num2str(r) ', p=' num2str(pval)]);

%% ------------------ SAUVEGARDE ------------------
% Cell purpose: optionally save figures and the theta struct. This block
% only runs if issave is true and expects datalc.save_path plus
% datalc.session_name to exist in the workspace.
if issave
    if ~exist(datalc.save_path, 'dir')
        warndlg('Le chemin de sauvegarde est manquant, les donn�es trait�es ne seront pas sauvegard�es.', 'Avertissement');
    else
        disp('Sauvegarde en cours');
        path_save = [datalc.save_path filesep datalc.session_name];
        if ~(exist(path_save, 'dir') == 7)
            mkdir(path_save);
            warning('Le dossier d''enregistrement n''existe pas et sera cr�� dans le chemin de sauvegarde');
        end
        if verbose
            for k = 1:length(h)
                fct_save_figure(h(k), [path_save filesep datalc.session_name '_theta' num2str(k)],'jpg');
            end
        end
        save([path_save filesep datalc.session_name '_Theta'], 'theta');
    end
end

varargout{1} = theta;
disp('----- Analyse Th�ta : TERMIN� -----');
