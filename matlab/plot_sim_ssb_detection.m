function plot_sim_ssb_detection(R, file)
%PLOT_SIM_SSB_DETECTION  Figure of sim_ssb_detection results (simulation, not measured).
%   plot_sim_ssb_detection(R, 'results/sim_ssb_detection.png')
%   plot_sim_ssb_detection('results/sim_ssb_detection.mat')     % redraw from the saved run
if ischar(R) || isstring(R)
    matFile = char(R);
    S = load(matFile, 'R');
    R = S.R;
    if nargin < 2, file = strrep(matFile, '.mat', '.png'); end
end
% categorical slots 1-3 of the reference palette (validated: CVD and normal-vision
% separation pass; aqua is below 3:1 contrast, hence direct labels, markers and the table)
cols = double([0x2a 0x78 0xd6; 0xeb 0x68 0x34; 0x1b 0xaf 0x7a])/255;
ink = [0.20 0.20 0.20];
mk = {'o', 's', '^'};
labels = {'PCI correct (per SSB)', 'MIB, one SSB', 'MIB in burst (combined)'};
labelY = [0.70, 0.30, 0.50];          % where each direct label sits on its curve
labelSide = [1, 1, -1];               % right of the curve, or left of it
f = figure('Visible', 'off', 'Color', [0.988 0.988 0.984], 'Position', [100 100 1100 430]);
tl = tiledlayout(f, 1, numel(R.mode), 'TileSpacing', 'compact', 'Padding', 'compact');
for im = 1:numel(R.mode)
    M = R.mode(im);
    ax = nexttile(tl);
    hold(ax, 'on');
    Y = [M.pPci, M.pMib, M.pBurstMibComb];
    h = gobjects(1, 3);
    for k = 1:3
        h(k) = plot(ax, M.snrDb, Y(:, k), '-', 'Color', cols(k, :), 'LineWidth', 2, 'Marker', mk{k}, ...
            'MarkerSize', 7, 'MarkerFaceColor', cols(k, :), 'MarkerEdgeColor', [1 1 1]);
        x = crossing(M.snrDb(:), Y(:, k), labelY(k));
        if ~isnan(x)
            if labelSide(k) > 0
                text(ax, x + 0.3, labelY(k), labels{k}, 'Color', ink, 'FontSize', 9, 'HorizontalAlignment', 'left');
            else
                text(ax, x - 0.3, labelY(k), labels{k}, 'Color', ink, 'FontSize', 9, 'HorizontalAlignment', 'right');
            end
        end
    end
    ax.YLim = [0 1.02];
    ax.XLim = [min(M.snrDb) max(M.snrDb)];
    ax.YGrid = 'on'; ax.XGrid = 'off';
    ax.GridColor = [0.88 0.88 0.88]; ax.GridAlpha = 1;
    ax.XColor = [0.35 0.35 0.35]; ax.YColor = [0.35 0.35 0.35];
    ax.Box = 'off'; ax.Color = 'none';
    xlabel(ax, 'SNR per resource element (dB)');
    if im == 1, ylabel(ax, 'Probability'); end
    title(ax, M.name, 'FontWeight', 'normal');
    if im == 1, legend(ax, h, labels, 'Location', 'northwest', 'Box', 'off', 'TextColor', ink); end
end
title(tl, sprintf('Simulated SSB detection, ESP32-C5 capture model (AWGN, %d bursts per point, not measured)', ...
    R.opt.Trials), 'FontSize', 11);
exportgraphics(f, file, 'Resolution', 150);
close(f);
end

function x = crossing(xs, ys, y0)
% first x where the curve rises through y0 (linear interpolation)
x = NaN;
i = find(ys(1:end-1) < y0 & ys(2:end) >= y0, 1);
if ~isempty(i)
    x = xs(i) + (y0 - ys(i))*(xs(i+1) - xs(i))/(ys(i+1) - ys(i));
end
end
