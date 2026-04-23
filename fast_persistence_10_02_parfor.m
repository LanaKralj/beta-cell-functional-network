%% edge_persistence_FAST_60s_fixed_CLEAN.m
% FAST edge persistence (temporal network) with short windows (60s)

clear; clc; close all;

% --- Parallel pool (open once) ---
if isempty(gcp('nocreate'))
    parpool(4);   % processes, 4 workerji
end

%% ===== Nastavitve =====
dataFile  = 'data_jki_2_ds_5.txt';   % [T x N] downsampled
coordFile = 'koordinate_jki_2.txt';  % [N x 2]
dt0       = load('sampling.txt');            % original dt (npr 0.1)
ds        = 5;                               % downsampling factor

dt  = dt0 * ds;
Fs  = 1/dt;
nyq = Fs/2;

cells = 2:169;   % analizirane celice
nC = numel(cells);

% globalni interval (sekunde)
tStart =1000;
tEnd   = 2800;

% FAST band (Hz) — mora biti < Nyquist!
fastBand = [0.02 0.15];% popravili spodnjo mejo  na 0.02-0.05

% okna
winSec  = 60;              % s
stepSec = 30;              % s

% top % povezav v vsakem oknu
keepTopPct = 3;            % %

% prag persistence za končno mrežo (0..1)
persistThr = 0.20;

% Alternativa: diskreten prag (rob mora biti prisoten vsaj v minWin oknih)
useMinWin = false;
minWin    = 2;

% wcoherence frequency limits (tale verzija tega ne podpira)
useFreqLimits = false;
%% ======================

% --- band clamp na Nyquist ---
fastBandUse = [fastBand(1), min(fastBand(2), 0.95*nyq)];
if fastBandUse(1) >= nyq
    warning('FAST band start (%.4f Hz) >= Nyquist (%.4f Hz). FAST bo prazen. Zmanjšaj ds ali popravi fastBand.', fastBandUse(1), nyq);
end
if fastBandUse(2) <= fastBandUse(1)
    warning('FAST band je degeneriran po clamping-u: [%.4f %.4f] Hz', fastBandUse(1), fastBandUse(2));
end

% FrequencyLimits je smiselno dati blizu banda
freqLimits = [max(0, fastBandUse(1)), max(fastBandUse(2), fastBandUse(1)+eps)];

fprintf('Fs=%.6f Hz, dt=%.6f s, Nyquist=%.6f Hz\n', Fs, dt, nyq);
fprintf('FAST band (orig): [%.4f %.4f] Hz | used: [%.4f %.4f] Hz\n', fastBand(1), fastBand(2), fastBandUse(1), fastBandUse(2));

%% 1) Load
tracesAll = load(dataFile);
coordsAll = load(coordFile);

coordsSel = coordsAll(cells,:);
x = coordsSel(:,1); y = coordsSel(:,2);

tAll = (0:size(tracesAll,1)-1) * dt;

idxGlobal = (tAll >= tStart) & (tAll <= tEnd);
traces = tracesAll(idxGlobal, cells);
t = tAll(idxGlobal);

fprintf('Global interval %.1f–%.1f s => %d samples\n', t(1), t(end), numel(t));

%% 2) Okna
edges = tStart:stepSec:(tEnd - winSec);
nW = numel(edges);
fprintf('FAST windows: %d (win=%.0fs, step=%.0fs)\n', nW, winSec, stepSec);

if nW < 3
    warning('Premalo oken (%d). Podaljšaj interval (tEnd) ali zmanjšaj winSec/stepSec.', nW);
end

U = triu(true(nC),1);
Pcount = zeros(nC,nC,'single');

nValid = 0;  % šteje samo okna, ki jih dejansko procesiramo

% --- RISANJE MREŽ PO OKNIH ---
doPlotEachWindow = true;     % vklop/izklop
plotEvery = 2;               % npr. 1 = vsako okno, 2 = vsako drugo ...
outDirWin = 'FAST_windows_networks_okno_60';
if doPlotEachWindow && ~exist(outDirWin,'dir')
    mkdir(outDirWin);
end


Csum = zeros(nC,nC,'single');   % akumulator za povprecno C matriko
%% 3) Sliding-window loop
for w = 1:nW
    tw = tic;   % START timer

    w0 = edges(w);
    w1 = w0 + winSec;

    idxW = (t >= w0) & (t < w1);
    X = traces(idxW, :);     % [Tw x nC]
    Tw = size(X,1);

    if Tw < 10
        warning('Okno %d prekratko (%d vzorcev), preskakujem.', w, Tw);
        continue;
    end

    C = zeros(nC, nC, 'single');

% ===== PARALLEL po vrsticah =====
parfor a = 1:nC
    xa = X(:,a);
    row = zeros(1, nC, 'single');  % parfor-safe lokalna vrstica

    for b = a+1:nC
        xb = X(:,b);

        % === wcoherence ===
        if useFreqLimits
            try
                [wcoh, ~, f] = wcoherence(xa, xb, Fs, 'FrequencyLimits', freqLimits);
            catch
                % če FrequencyLimits ni podprt, preklopi lokalno (v parfor ne spreminjaj globalne spremenljivke!)
                [wcoh, ~, f] = wcoherence(xa, xb, Fs);
            end
        else
            [wcoh, ~, f] = wcoherence(xa, xb, Fs);
        end

        fHz = f(:);

        % --- alignment fHz <-> wcoh ---
        nF = size(wcoh,1);
        if numel(fHz) ~= nF
            n = min(numel(fHz), nF);
            fHz  = fHz(1:n);
            wcoh = wcoh(1:n,:);
        end

        idxF = (fHz >= fastBandUse(1)) & (fHz <= fastBandUse(2));
        if any(idxF)
            m = mean(wcoh(idxF,:), 'all', 'omitnan');
            if isnan(m), m = 0; end
            row(b) = single(m);  % samo upper triangle
        else
            row(b) = 0;
        end

        % (Debug v parfor: ne printaj pogosto; če res rabiš, shrani zastavico)
    end

    C(a,:) = row;  % vsak worker napiše svojo vrstico
end

% Simetrizacija + diagonala
C = C + C.';
C(1:nC+1:end) = 0;
Csum = Csum + C;
    % obdrži top keepTopPct% povezav v tem oknu
    A = keepTopPctSym(C, keepTopPct, U);

    % prištej persistence count
    Pcount = Pcount + single(A > 0);
    nValid = nValid + 1;

    % --- Plot network for this window ---
if doPlotEachWindow && mod(w, plotEvery) == 0
    % node strength in this window (utežena mreža A)
    strength_w = sum(A, 2);

    Gw = graph(A, 'upper');

    fig = figure('Visible','off','Color','w');
    pw = plot(Gw, 'XData', x, 'YData', y);
    pw.NodeLabel = {};
    pw.MarkerSize = 10;
    pw.NodeCData = strength_w;
    colorbar;

    if ~isempty(Gw.Edges)
        pw.LineWidth = 0.5 + 5*normalize(Gw.Edges.Weight,'range');
        pw.EdgeAlpha = 0.6;
    end

    axis equal; grid on;
    title(sprintf('FAST window %d/%d: %.1f–%.1f s | top %.1f%%', ...
        w, nW, w0, w1, keepTopPct));
    xlabel('x'); ylabel('y');

    fn = fullfile(outDirWin, sprintf('FAST_win_%03d_t%04d_%04d.png', w, round(w0), round(w1)));
    exportgraphics(fig, fn, 'Resolution', 200);
    close(fig);
end

    % Debug: koliko robov + porazdelitev
    valsPos = C(U);
    valsPos = valsPos(~isnan(valsPos) & valsPos>0);
    if isempty(valsPos)
        q95 = NaN; mx = 0;
    else
        q95 = quantile(valsPos, 0.95);
        mx  = max(valsPos);
    end
    nEdges = nnz(triu(A,1));
    fprintf('  window %d/%d (%.1f–%.1f s): kept edges=%d | maxC=%.3f | q95=%.3f\n|time=%.2fs\n', ...
        w, nW, w0, w1, nEdges, mx, q95, toc (tw));
end

fprintf('Valid windows used: %d / %d\n', nValid, nW);
if nValid == 0
    error('Ni veljavnih oken (nValid=0). Preveri tStart/tEnd/winSec/stepSec.');
end
%povprecna coupling matr
Cavg = Csum / nValid;

figure('Color','w');
imagesc(Cavg); axis image;
colormap(parula); colorbar;
caxis([0 1]);
title(sprintf('Average FAST wavelet coupling | band [%.2f %.2f] Hz', ...
    fastBandUse(1), fastBandUse(2)));
xlabel('cell index');
ylabel('cell index');

exportgraphics(gcf, 'FAST_Cavg.png', 'Resolution', 200);

%% 4) Persistence
P = Pcount / nValid;
P = triu(P,1) + triu(P,1)';
P(1:nC+1:end) = 0;

% končna mreža
A_persist = P;
if useMinWin
    A_persist(Pcount < minWin) = 0;
else
    A_persist(A_persist < persistThr) = 0;
end
A_persist(1:nC+1:end) = 0;

strength = sum(A_persist,2);


%% 4b) Izbor reprezentativnih FAST parov

% --- upper triangle ---
Uonly = triu(true(nC),1);

% ============================================================
% A) DOBRO POVEZAN PAR = maksimum v A_persist
% ============================================================
Aup = A_persist;
Aup(~Uonly) = 0;

if nnz(Aup) == 0
    error('A_persist nima nobene nenicelne povezave.');
end

[~, idxMax] = max(Aup(:));
[iConn, jConn] = ind2sub(size(Aup), idxMax);

% razdalja tega para
dConn = hypot(x(iConn)-x(jConn), y(iConn)-y(jConn));

fprintf('\n=== DOBRO POVEZAN FAST PAR ===\n');
fprintf('local idx: (%d, %d)\n', iConn, jConn);
fprintf('cell IDs : (%d, %d)\n', cells(iConn), cells(jConn));
fprintf('A_persist: %.4f\n', A_persist(iConn,jConn));
fprintf('Cavg     : %.4f\n', Cavg(iConn,jConn));
fprintf('distance : %.2f\n', dConn);

% ============================================================
% B) NEPOVEZAN PAR
% ============================================================
bestScore = inf;
iDisc = NaN;
jDisc = NaN;

for i = 1:nC
    for j = i+1:nC
        if A_persist(i,j) == 0
            dij = hypot(x(i)-x(j), y(i)-y(j));

            % kombiniran score:
            score = abs(dij - dConn) + 100*Cavg(i,j);

            if score < bestScore
                bestScore = score;
                iDisc = i;
                jDisc = j;
            end
        end
    end
end  % <-- pomemben end za zanko!

% izpis (ZUNAJ zanke!)
if ~isnan(iDisc)
    dDisc = hypot(x(iDisc)-x(jDisc), y(iDisc)-y(jDisc));

    fprintf('\n=== NEPOVEZAN FAST PAR ===\n');
    fprintf('local idx: (%d, %d)\n', iDisc, jDisc);
    fprintf('cell IDs : (%d, %d)\n', cells(iDisc), cells(jDisc));
    fprintf('A_persist: %.4f\n', A_persist(iDisc,jDisc));
    fprintf('Cavg     : %.4f\n', Cavg(iDisc,jDisc));
    fprintf('distance : %.2f\n', dDisc);
else
    warning('Nisem nasel nepovezanega FAST para.');
end

%% 4c) Scalogram dobro povezanega FAST para

% signala iz izbranega para
xa_conn = traces(:, iConn);
xb_conn = traces(:, jConn);

% wavelet coherence
if useFreqLimits
    try
        [wcoh_conn, ~, f_conn, coi_conn] = wcoherence(xa_conn, xb_conn, Fs, 'FrequencyLimits', freqLimits);
    catch
        [wcoh_conn, ~, f_conn, coi_conn] = wcoherence(xa_conn, xb_conn, Fs);
    end
else
    [wcoh_conn, ~, f_conn, coi_conn] = wcoherence(xa_conn, xb_conn, Fs);
end

% poravnava dimenzij, ce je treba
f_conn = f_conn(:);
coi_conn = coi_conn(:);

nF = size(wcoh_conn, 1);
if numel(f_conn) ~= nF
    n = min(numel(f_conn), nF);
    f_conn = f_conn(1:n);
    wcoh_conn = wcoh_conn(1:n, :);
end

nT = size(wcoh_conn, 2);
if numel(coi_conn) ~= nT
    n = min(numel(coi_conn), nT);
    coi_conn = coi_conn(1:n);
    wcoh_conn = wcoh_conn(:, 1:n);
end

% izris
figScal = figure('Color','w','Name','FAST connected pair scalogram');
surface(t, f_conn, wcoh_conn, 'EdgeColor', 'none');
view(0,90);
axis tight;
ylim([fastBandUse(1) fastBandUse(2)]);
colormap(parula);
colorbar;
hold on;

% FAST band
yline(fastBandUse(1), 'r--', 'LineWidth', 1.2);
yline(fastBandUse(2), 'r--', 'LineWidth', 1.2);

% cone of influence
plot(t(1:numel(coi_conn)), coi_conn, 'w--', 'LineWidth', 1.5);

hold off;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(sprintf('FAST connected pair: cells %d-%d | P = %.3f | Cavg = %.3f', ...
    cells(iConn), cells(jConn), A_persist(iConn,jConn), Cavg(iConn,jConn)));

% po zelji omejis y-os, da bolje vidis FAST pas
ylim([0, min(0.5, max(f_conn))]);

exportgraphics(figScal, ...
    sprintf('FAST_scalogram_connected_cells_%d_%d.png', cells(iConn), cells(jConn)), ...
    'Resolution', 250);

%% 4d) Scalogram NEPOVEZANEGA FAST para

xa_disc = traces(:, iDisc);
xb_disc = traces(:, jDisc);

if useFreqLimits
    try
        [wcoh_disc, ~, f_disc, coi_disc] = wcoherence(xa_disc, xb_disc, Fs, 'FrequencyLimits', freqLimits);
    catch
        [wcoh_disc, ~, f_disc, coi_disc] = wcoherence(xa_disc, xb_disc, Fs);
    end
else
    [wcoh_disc, ~, f_disc, coi_disc] = wcoherence(xa_disc, xb_disc, Fs);
end

f_disc = f_disc(:);
coi_disc = coi_disc(:);

% alignment
nF = size(wcoh_disc,1);
if numel(f_disc) ~= nF
    n = min(numel(f_disc), nF);
    f_disc = f_disc(1:n);
    wcoh_disc = wcoh_disc(1:n,:);
end

nT = size(wcoh_disc,2);
if numel(coi_disc) ~= nT
    n = min(numel(coi_disc), nT);
    coi_disc = coi_disc(1:n);
    wcoh_disc = wcoh_disc(:,1:n);
end

% izris
figScalDisc = figure('Color','w','Name','FAST disconnected pair scalogram');
surface(t, f_disc, wcoh_disc, 'EdgeColor', 'none');
view(0,90);
ylim([fastBandUse(1) fastBandUse(2)]);
axis tight;
colormap(parula);
colorbar;
hold on;

% FAST band
yline(fastBandUse(1), 'r--', 'LineWidth', 1.2);
yline(fastBandUse(2), 'r--', 'LineWidth', 1.2);

% cone of influence
plot(t(1:numel(coi_disc)), coi_disc, 'w--', 'LineWidth', 1.5);

hold off;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(sprintf('FAST disconnected pair: cells %d-%d | P = %.3f | Cavg = %.3f', ...
    cells(iDisc), cells(jDisc), A_persist(iDisc,jDisc), Cavg(iDisc,jDisc)));

ylim([0, min(0.5, max(f_disc))]);

exportgraphics(figScalDisc, ...
    sprintf('FAST_scalogram_disconnected_cells_%d_%d.png', cells(iDisc), cells(jDisc)), ...
    'Resolution', 250);


%% 4e) SIDE-BY-SIDE: connected vs disconnected

fig = figure('Color','w','Position',[100 100 1200 450]);

% =========================
% CONNECTED (levo)
% =========================
subplot(1,2,1)

surface(t, f_conn, wcoh_conn, 'EdgeColor','none');
view(0,90);
axis tight;


%hold on;
%yline(fastBandUse(1), 'w--', 'LineWidth', 1.2);
%yline(fastBandUse(2), 'w--', 'LineWidth', 1.2);
%plot(t(1:numel(coi_conn)), coi_conn, 'k--', 'LineWidth', 1.2);
%hold off;

title(sprintf('Connected pair'));

xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([0.02 0.15]);

% =========================
% DISCONNECTED (desno)
% =========================
subplot(1,2,2)

surface(t, f_disc, wcoh_disc, 'EdgeColor','none');
view(0,90);
axis tight;

colormap(turbo);
caxis([0 1]);
colorbar;

%hold on;
%yline(fastBandUse(1), 'w--', 'LineWidth', 1.2);
%yline(fastBandUse(2), 'w--', 'LineWidth', 1.2);
%plot(t(1:numel(coi_disc)), coi_disc, 'k--', 'LineWidth', 1.2);
%hold off;

title(sprintf('Disconnected pair'));

xlabel('Time (s)');
ylim([0.02 0.15]);

% =========================
% skupni naslov
% =========================
sgtitle('Fast wavelet coherence');
set(findall(fig,'-property','FontSize'),'FontSize',16)
% shrani
exportgraphics(fig, 'FAST_connected_vs_disconnected.png', 'Resolution', 300);

%% Hub cells (top K by strength)
K = min(5, nC);
[~, ix] = sort(strength, 'descend');
hubIdx_local = ix(1:K);
hubCells = cells(hubIdx_local);

disp('Hub celice (ID):');
disp(hubCells);



%% 5) Plot persistence network
G = graph(A_persist, 'upper');

figure; clf; set(gcf,'Color','w','Name','FAST persistence network (60s)');
p = plot(G, 'XData', x, 'YData', y);
p.NodeLabel = {};
p.MarkerSize = 10;
p.NodeCData = strength;
colorbar;

if ~isempty(G.Edges)
    p.LineWidth = 0.5 + 5*normalize(G.Edges.Weight,'range');
    p.EdgeAlpha = 0.6;
end

axis equal; grid on;
if useMinWin
    title(sprintf('FAST edge persistence (%.0f–%.0f s, win=%.0fs step=%.0fs, top %.1f%%, minWin=%d)', ...
        tStart, tEnd, winSec, stepSec, keepTopPct, minWin));
else
    title(sprintf('FAST edge persistence (%.0f–%.0f s, win=%.0fs step=%.0fs, top %.1f%%, thr %.0f%%)', ...
        tStart, tEnd, winSec, stepSec, keepTopPct, 100*persistThr));
end
xlabel('x'); ylabel('y');

fprintf('DEBUG Fs=%.6f Hz | Nyquist=%.6f Hz | fastBandUse=[%.4f %.4f] Hz\n', Fs, Fs/2, fastBandUse(1), fastBandUse(2));
fprintf('Nonzero P count: %d | Max persistence: %.3f | Edges in A_persist: %d\n', ...
    nnz(triu(P,1)), max(P(:)), nnz(triu(A_persist,1)));

%% 6) Shrani
outName = sprintf('FAST_persistence_%ds_%ds_win%.0f_step%.0f_top%.1f_thr%.0f.mat', ...
    round(tStart), round(tEnd), winSec, stepSec, keepTopPct, 100*persistThr);
save(outName, 'P','Pcount','A_persist','cells','hubCells','strength','x','y','coordsSel','Cavg',...
    'tStart','tEnd','winSec','stepSec','keepTopPct','persistThr','useMinWin','minWin','Fs','dt','fastBand','fastBandUse','-v7.3');
fprintf('Saved: %s\n', outName);



%% ===== local function =====
function A = keepTopPctSym(C, pct, U)
% Obdrži top pct% povezav v upper-tri in vrni simetrično uteženo matriko.
    A = C;
    vals = A(U);
    vals = vals(~isnan(vals));
    vals = vals(vals > 0);

    if isempty(vals)
        A(:) = 0;
        return;
    end

    thr = quantile(vals, 1 - pct/100);

    A(~U) = 0;
    A(A < thr) = 0;

    A = triu(A,1) + triu(A,1)';
    A(1:size(A,1)+1:end) = 0;
end

out = [(0:nC-1)', strength(:)];
writematrix(out, 'FAST_strength_byIndex0_jkl_5_3%posto_okno_60_0_15.txt', 'Delimiter',' ');