clear;
clc;
close all;

fprintf(['Running the pure-script version of revised Figure 7. ' ...
    'No local or external functions are used.\n']);

%% =========================================================
% Figure 7: SLA-constrained idle backfilling
%
% This is a pure MATLAB script. It contains no local functions,
% subfunctions, or external function calls.
%
% q^E(v,eta) = max {q in [0,eta*S_G]:
%                   S_G-q-D_h(q;v) >= b_h}
%
% Panel (a): q^E versus q^B = eta*S_G
% Panel (b): correction of 0/P/F market-state migration
% Panel (c): Delta SW = SW^E - SW^B
%% =========================================================

%% 1. Unified calibration

nTotal = 100;
alpha = 0.30;

nH = alpha * nTotal;
nL = (1 - alpha) * nTotal;

tH = 1;
tL = 1;

pH = 2;
pL = 1.5;

gammaH = 20;
gammaL = 12;

SG = 50;
SC = 100;

% SLA requirement:
% W_h <= WbarH  <=>  S_G-q-D_h(q;v) >= bH.

WbarH = 0.10;
bH = 1 / WbarH;

etaMax = 0.60;

vMin = 0;
vMax = 6;
vStep = 0.02;
vGrid = vMin:vStep:vMax;

etaStep = 0.01;
etaGrid = 0:etaStep:etaMax;

etaLineList = [0.20,0.30,0.40];

nV = length(vGrid);
nEta = length(etaGrid);
nEtaLine = length(etaLineList);

tol = 1.0e-10;
checkTol = 1.0e-7;
bisectionIterations = 80;

AH = nH / tH;
AL = 2 * nL / tL;

%% 2. Parameter checks

assert(nH > 0 && nL > 0, ...
    'Both market sizes must be positive.');

assert(tH > 0 && tL > 0, ...
    'Matching-cost parameters must be positive.');

assert(SG > 0 && SC > 0, ...
    'Capacity parameters must be positive.');

assert(WbarH > 0 && bH > 0, ...
    'The SLA limit and reserve must be positive.');

assert(etaMax > 0 && etaMax < 1, ...
    'etaMax must lie strictly between zero and one.');

%% 3. Storage

qFixedAll = nan(nEta,nV);
qIdleAll = nan(nEta,nV);
qGapAll = nan(nEta,nV);

SWFixedAll = nan(nEta,nV);
SWIdleAll = nan(nEta,nV);
DeltaSW = nan(nEta,nV);

regHFixedAll = nan(nEta,nV);
regLFixedAll = nan(nEta,nV);
regHIdleAll = nan(nEta,nV);
regLIdleAll = nan(nEta,nV);

baselineFeasible = false(nEta,nV);
slaBinding = false(nEta,nV);

% Correction categories:
% 1 = no 0/P/F state correction
% 2 = high-tier state correction only
% 3 = low-tier state correction only
% 4 = both tiers corrected

correctionType = nan(nEta,nV);

%% 4. Main grid calculation

for e = 1:nEta

    eta = etaGrid(e);

    for j = 1:nV

        v = vGrid(j);
        Vh = v - pH;
        Vl = v - pL;

        qFixed = eta * SG;
        qFixedAll(e,j) = qFixed;

        %% -------------------------------------------------
        % 4.1 Evaluate the SLA reserve at q = 0 and q = qFixed
        %% -------------------------------------------------

        qProbe = [0,qFixed];
        reserveProbe = nan(1,2);

        for probe = 1:2

            qTest = qProbe(probe);
            sHTest = SG - qTest;

            if sHTest <= 0

                dHTest = NaN;

            elseif Vh - gammaH / sHTest <= tol

                dHTest = 0;

            elseif sHTest - nH > tol ...
                    && Vh - tH ...
                    - gammaH / (sHTest - nH) >= -tol

                dHTest = nH;

            else

                radicandH = ...
                    (sHTest - AH * Vh)^2 ...
                    + 4 * AH * gammaH;

                dHTest = ...
                    (sHTest ...
                    + AH * Vh ...
                    - sqrt(max(radicandH,0))) / 2;

                dHUpper = ...
                    min(nH,sHTest - max(1.0e-10,1.0e-9 * sHTest));

                dHTest = min(max(dHTest,0),dHUpper);

            end

            reserveProbe(probe) = sHTest - dHTest;

        end

        reserveAtZero = reserveProbe(1);
        reserveAtCap = reserveProbe(2);

        if reserveAtZero < bH - checkTol

            % The no-backfilling baseline itself violates the SLA.
            % This point is excluded rather than silently assigned q = 0.
            continue;

        end

        baselineFeasible(e,j) = true;

        %% -------------------------------------------------
        % 4.2 Largest SLA-feasible q
        %% -------------------------------------------------

        if qFixed <= tol

            qIdle = 0;
            isBinding = false;

        elseif reserveAtCap >= bH - checkTol

            qIdle = qFixed;
            isBinding = false;

        else

            % The SLA boundary is monotone in q, so bisection
            % identifies the largest feasible q without any function call.

            qLower = 0;
            qUpper = qFixed;

            for bisectIter = 1:bisectionIterations

                qTest = 0.5 * (qLower + qUpper);
                sHTest = SG - qTest;

                if Vh - gammaH / sHTest <= tol

                    dHTest = 0;

                elseif sHTest - nH > tol ...
                        && Vh - tH ...
                        - gammaH / (sHTest - nH) >= -tol

                    dHTest = nH;

                else

                    radicandH = ...
                        (sHTest - AH * Vh)^2 ...
                        + 4 * AH * gammaH;

                    dHTest = ...
                        (sHTest ...
                        + AH * Vh ...
                        - sqrt(max(radicandH,0))) / 2;

                    dHUpper = ...
                        min(nH,sHTest ...
                        - max(1.0e-10,1.0e-9 * sHTest));

                    dHTest = min(max(dHTest,0),dHUpper);

                end

                reserveTest = sHTest - dHTest;

                if reserveTest >= bH
                    qLower = qTest;
                else
                    qUpper = qTest;
                end

            end

            qIdle = qLower;
            isBinding = true;

        end

        qIdleAll(e,j) = qIdle;
        qGapAll(e,j) = qFixed - qIdle;
        slaBinding(e,j) = isBinding;

        if qIdle > qFixed + checkTol

            error(['The SLA-constrained q exceeds the fixed cap at ' ...
                'eta = %.4f and v = %.4f.'],eta,v);

        end

        %% -------------------------------------------------
        % 4.3 Calculate outcomes under fixed and idle backfilling
        %% -------------------------------------------------

        qMechanism = [qFixed,qIdle];

        welfareTemp = nan(1,2);
        regHTemp = nan(1,2);
        regLTemp = nan(1,2);

        for mechanism = 1:2

            qUse = qMechanism(mechanism);

            sH = SG - qUse;
            sL = SC + qUse;

            %% High-tier market

            if Vh - gammaH / sH <= tol

                dH = 0;
                csH = 0;
                regH = 0;

            elseif sH - nH > tol ...
                    && Vh - tH - gammaH / (sH - nH) >= -tol

                dH = nH;

                csH = ...
                    nH ...
                    * (Vh ...
                    - gammaH / (sH - nH) ...
                    - tH / 2);

                csH = max(csH,0);
                regH = 2;

            else

                radicandH = ...
                    (sH - AH * Vh)^2 ...
                    + 4 * AH * gammaH;

                dH = ...
                    (sH ...
                    + AH * Vh ...
                    - sqrt(max(radicandH,0))) / 2;

                dHUpper = ...
                    min(nH,sH - max(1.0e-10,1.0e-9 * sH));

                if dH < -tol || dH > dHUpper + checkTol

                    error(['Infeasible high-tier demand at eta = %.4f ' ...
                        'and v = %.4f.'],eta,v);

                end

                dH = min(max(dH,0),dHUpper);

                residualH = ...
                    Vh ...
                    - dH / AH ...
                    - gammaH / (sH - dH);

                if abs(residualH) ...
                        > checkTol * max(1,abs(Vh))

                    error(['High-tier equilibrium residual exceeds ' ...
                        'the tolerance at eta = %.4f and v = %.4f.'], ...
                        eta,v);

                end

                csH = dH^2 / (2 * AH);
                regH = 1;

            end

            %% Low-tier market

            if Vl - gammaL / sL <= tol

                dL = 0;
                csL = 0;
                regL = 0;

            elseif sL - nL > tol ...
                    && Vl - tL / 2 ...
                    - gammaL / (sL - nL) >= -tol

                dL = nL;

                csL = ...
                    nL ...
                    * (Vl ...
                    - gammaL / (sL - nL) ...
                    - tL / 4);

                csL = max(csL,0);
                regL = 2;

            else

                radicandL = ...
                    (sL - AL * Vl)^2 ...
                    + 4 * AL * gammaL;

                dL = ...
                    (sL ...
                    + AL * Vl ...
                    - sqrt(max(radicandL,0))) / 2;

                dLUpper = ...
                    min(nL,sL - max(1.0e-10,1.0e-9 * sL));

                if dL < -tol || dL > dLUpper + checkTol

                    error(['Infeasible low-tier demand at eta = %.4f ' ...
                        'and v = %.4f.'],eta,v);

                end

                dL = min(max(dL,0),dLUpper);

                residualL = ...
                    Vl ...
                    - dL / AL ...
                    - gammaL / (sL - dL);

                if abs(residualL) ...
                        > checkTol * max(1,abs(Vl))

                    error(['Low-tier equilibrium residual exceeds ' ...
                        'the tolerance at eta = %.4f and v = %.4f.'], ...
                        eta,v);

                end

                csL = dL^2 / (2 * AL);
                regL = 1;

            end

            profit = pH * dH + pL * dL;
            welfareTemp(mechanism) = profit + csH + csL;
            regHTemp(mechanism) = regH;
            regLTemp(mechanism) = regL;

        end

        SWFixedAll(e,j) = welfareTemp(1);
        SWIdleAll(e,j) = welfareTemp(2);
        DeltaSW(e,j) = welfareTemp(2) - welfareTemp(1);

        regHFixedAll(e,j) = regHTemp(1);
        regLFixedAll(e,j) = regLTemp(1);
        regHIdleAll(e,j) = regHTemp(2);
        regLIdleAll(e,j) = regLTemp(2);

        %% -------------------------------------------------
        % 4.4 Market-state correction category
        %% -------------------------------------------------

        highCorrected = regHTemp(2) > regHTemp(1);
        lowCorrected = regLTemp(2) < regLTemp(1);

        if regHTemp(2) < regHTemp(1) ...
                || regLTemp(2) > regLTemp(1)

            error(['A reverse state correction was detected at ' ...
                'eta = %.4f and v = %.4f.'],eta,v);

        end

        correctionType(e,j) = ...
            1 ...
            + double(highCorrected) ...
            + 2 * double(lowCorrected);

    end

end

%% 5. Extract line data for Panel (a)

qIdleLines = nan(nEtaLine,nV);
qFixedLines = nan(nEtaLine,nV);

for kLine = 1:nEtaLine

    etaIndex = round(etaLineList(kLine) / etaStep) + 1;

    qIdleLines(kLine,:) = ...
        qIdleAll(etaIndex,:) / SG;

    qFixedLines(kLine,:) = ...
        qFixedAll(etaIndex,:) / SG;

end

%% 6. Numerical diagnostics

feasibleMask = baselineFeasible;

assert(any(feasibleMask(:)), ...
    'No baseline-feasible grid points were found.');

assert(all(isfinite(qIdleAll(feasibleMask))), ...
    'Nonfinite q values were detected.');

assert(all(isfinite(DeltaSW(feasibleMask))), ...
    'Nonfinite welfare differences were detected.');

assert(all(qGapAll(feasibleMask) >= -checkTol), ...
    'The numerical result violates q^E <= q^B.');

bindingCount = nnz(slaBinding & feasibleMask);
feasibleCount = nnz(feasibleMask);

fprintf('\nPure-script Figure 7 checks passed.\n');
fprintf('SLA reserve b_h: %.4f\n',bH);
fprintf('Waiting-time limit Wbar_h: %.4f\n',WbarH);
fprintf('Feasible grid points: %d\n',feasibleCount);
fprintf('SLA-binding grid points: %d (%.2f%%)\n', ...
    bindingCount,100 * bindingCount / feasibleCount);
fprintf('Maximum q^B-q^E: %.6f\n', ...
    max(qGapAll(feasibleMask)));
fprintf('Delta SW range: [%.6f, %.6f]\n', ...
    min(DeltaSW(feasibleMask)), ...
    max(DeltaSW(feasibleMask)));

for category = 1:4

    fprintf('Correction category %d: %d points\n', ...
        category,nnz(correctionType == category));

end

%% 7. Plot settings

fontName = 'Times New Roman';

lineColors = [ ...
    0.0000,0.4470,0.7410;
    0.8500,0.3250,0.0980;
    0.4940,0.1840,0.5560 ...
    ];

[X,Y] = meshgrid(vGrid,etaGrid);

fig7 = figure( ...
    'Color','w', ...
    'Position',[30,100,1650,520], ...
    'Renderer','painters');

%% 8. Panel (a): q^E versus q^B

ax1 = subplot(1,3,1,'Parent',fig7);
hold(ax1,'on');

legendLabels = cell(1,2 * nEtaLine);

for kLine = 1:nEtaLine

    currentColor = lineColors(kLine,:);
    eta = etaLineList(kLine);

    plot(ax1,vGrid,qIdleLines(kLine,:), ...
        '-', ...
        'Color',currentColor, ...
        'LineWidth',2.2);

    plot(ax1,vGrid,qFixedLines(kLine,:), ...
        '--', ...
        'Color',currentColor, ...
        'LineWidth',1.6);

    legendLabels{2 * kLine - 1} = ...
        sprintf('$q^E/S_G, \\eta=%.2f$',eta);

    legendLabels{2 * kLine} = ...
        sprintf('$q^B/S_G, \\eta=%.2f$',eta);

end

xlabel(ax1,'Service value  $v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',14);

ylabel(ax1,'Normalized backfilling volume', ...
    'FontName',fontName, ...
    'FontSize',14);

title(ax1,'(a) SLA-constrained backfilling', ...
    'FontName',fontName, ...
    'FontSize',13, ...
    'FontWeight','bold');

set(ax1, ...
    'FontName',fontName, ...
    'FontSize',11, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','on', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'GridAlpha',0.15);

legend(ax1,legendLabels, ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',9, ...
    'Location','southwest', ...
    'Box','on');

xlim(ax1,[vMin,vMax]);
ylim(ax1,[0,max(etaLineList) + 0.05]);

%% 9. Panel (b): market-state correction

ax2 = subplot(1,3,2,'Parent',fig7);

imagesc(ax2,vGrid,etaGrid,correctionType);
axis(ax2,'xy');
hold(ax2,'on');

correctionColors = [ ...
    0.88,0.88,0.88;
    0.38,0.64,0.82;
    0.95,0.72,0.30;
    0.48,0.74,0.52 ...
    ];

colormap(ax2,correctionColors);
caxis(ax2,[0.5,4.5]);

contour(ax2,X,Y,double(slaBinding),[0.5,0.5], ...
    'k-', ...
    'LineWidth',1.5);

xlabel(ax2,'Service value  $v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',14);

ylabel(ax2,'Institutional cap  $\eta$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',14);

title(ax2,'(b) Correction of market-state migration', ...
    'FontName',fontName, ...
    'FontSize',13, ...
    'FontWeight','bold');

set(ax2, ...
    'FontName',fontName, ...
    'FontSize',11, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','on', ...
    'Layer','top');

xlim(ax2,[vMin,vMax]);
ylim(ax2,[0,etaMax]);

cb2 = colorbar(ax2);
cb2.Ticks = 1:4;
cb2.TickLabels = { ...
    'No state correction', ...
    'High-tier correction', ...
    'Low-tier correction', ...
    'Both tiers corrected' ...
    };
cb2.FontName = fontName;
cb2.FontSize = 9;

%% 10. Panel (c): welfare correction

ax3 = subplot(1,3,3,'Parent',fig7);

imagesc(ax3,vGrid,etaGrid,DeltaSW);
axis(ax3,'xy');
hold(ax3,'on');

nColor = 256;
nHalf = nColor / 2;

blueSide = [ ...
    linspace(0.20,1.00,nHalf)', ...
    linspace(0.42,1.00,nHalf)', ...
    linspace(0.72,1.00,nHalf)' ...
    ];

redSide = [ ...
    linspace(1.00,0.75,nHalf)', ...
    linspace(1.00,0.18,nHalf)', ...
    linspace(1.00,0.18,nHalf)' ...
    ];

divergingMap = [blueSide;redSide];
colormap(ax3,divergingMap);

finiteDelta = DeltaSW(isfinite(DeltaSW));
colorLimit = max(abs(finiteDelta));

if colorLimit <= tol
    colorLimit = 1;
end

caxis(ax3,[-colorLimit,colorLimit]);

contour(ax3,X,Y,double(slaBinding),[0.5,0.5], ...
    'k--', ...
    'LineWidth',1.5);

xlabel(ax3,'Service value  $v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',14);

ylabel(ax3,'Institutional cap  $\eta$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',14);

title(ax3,'(c) Welfare correction', ...
    'FontName',fontName, ...
    'FontSize',13, ...
    'FontWeight','bold');

set(ax3, ...
    'FontName',fontName, ...
    'FontSize',11, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','on', ...
    'Layer','top');

xlim(ax3,[vMin,vMax]);
ylim(ax3,[0,etaMax]);

cb3 = colorbar(ax3);
cb3.Label.String = '$\Delta SW=SW^E-SW^B$';
cb3.Label.Interpreter = 'latex';
cb3.Label.FontName = fontName;
cb3.Label.FontSize = 12;
cb3.FontName = fontName;
cb3.FontSize = 9;

%% 11. Export using version-compatible print commands

drawnow;

print(fig7, ...
    'Figure7_idle_backfilling_pure_script.pdf', ...
    '-dpdf', ...
    '-painters');

print(fig7, ...
    'Figure7_idle_backfilling_pure_script.png', ...
    '-dpng', ...
    '-r600');

fprintf(['Figure 7 has been exported as PDF and PNG in the ' ...
    'current MATLAB folder.\n']);
