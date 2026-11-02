clear;
clc;
close all;

%% =========================================================
% Figure 5: Incremental effects of fixed backfilling
%
% Baseline: post-merger strict segregation, eta = 0
% Comparison: post-merger fixed backfilling, eta = etaBar
%
% Delta Y(v,etaBar) = Y(v,etaBar) - Y(v,0)
%
% The script generates four independent panels:
%   Figure6a: Delta high-tier consumer surplus
%   Figure6b: Delta low-tier consumer surplus
%   Figure6c: Delta platform profit
%   Figure6d: Delta social welfare
%% =========================================================

%% 1. Unified calibration

nTotal = 100;
alpha = 0.30;

nH = alpha * nTotal;
nL = (1 - alpha) * nTotal;

tH = 1;
tL = 1;

% Pre-merger capacities are retained only to preserve the
% unified calibration used in the other numerical figures.
sHPre = 50;
sLPre = 50;

SG = 50;
SC = 2 * sLPre;              % 100

pH = 2;
pL = 1.5;

vMin = 0;
vMax = 6;
vStep = 0.02;
vGrid = vMin:vStep:vMax;
nV = length(vGrid);

ghBench = 20;
glBench = 12;

deltaH = 0;
deltaL = 0;

etaBarList = [0.05, 0.10, 0.20, 0.30];
nEta = length(etaBarList);

% The first row is the strict-segregation benchmark eta = 0.
etaAll = [0, etaBarList];
nEtaAll = length(etaAll);

% Numerical tolerances
tol = 1.0e-10;
capacityTolFactor = 1.0e-8;

%% 2. Parameter validation

assert(nTotal > 0, ...
    'nTotal must be positive.');

assert(alpha > 0 && alpha < 1, ...
    'alpha must lie strictly between zero and one.');

assert(nH > 0 && nL > 0, ...
    'Both market sizes must be positive.');

assert(tH > 0 && tL > 0, ...
    'Differentiation parameters must be positive.');

assert(SG > 0 && SC > 0, ...
    'Post-merger capacities must be positive.');

assert(pH >= 0 && pL >= 0, ...
    'Prices must be nonnegative.');

assert(ghBench >= 0 && glBench >= 0, ...
    'Congestion-sensitivity parameters must be nonnegative.');

assert(all(etaBarList > 0 & etaBarList < 1), ...
    'Each etaBar value must lie strictly between zero and one.');

assert(all(diff(etaBarList) > 0), ...
    'etaBarList must be strictly increasing.');

%% 3. Storage

ProfitAll = zeros(nEtaAll,nV);
CSHAll = zeros(nEtaAll,nV);
CSLAll = zeros(nEtaAll,nV);
WelfareAll = zeros(nEtaAll,nV);

DHAll = zeros(nEtaAll,nV);
DLAll = zeros(nEtaAll,nV);

HighRegimeAll = strings(nEtaAll,nV);
LowRegimeAll = strings(nEtaAll,nV);

%% 4. Simulation over eta and v

for e = 1:nEtaAll

    eta = etaAll(e);

    %% Effective capacities under fixed backfilling

    sH = (1 - eta) * SG;
    sL = SC + eta * SG;

    if sH <= 0
        error('High-tier capacity is nonpositive at eta = %.4f.',eta);
    end

    for idx = 1:nV

        v = vGrid(idx);

        %% Net service values

        Vh = v - pH + deltaH;
        Vl = v - pL + deltaL;

        %% =================================================
        % 4.1 Post-merger high-tier equilibrium
        %
        % Marginal-consumer equation:
        %
        % Vh - tH*dH/nH - ghBench/(sH-dH) = 0
        %% =================================================

        highUtilityAtZero = ...
            Vh - ghBench / sH;

        if highUtilityAtZero <= tol

            highReg = "0h";
            dH = 0;

        else

            highFullFeasible = ...
                sH - nH > tol;

            if highFullFeasible

                highUtilityAtFull = ...
                    Vh ...
                    - tH ...
                    - ghBench / (sH - nH);

            else

                highUtilityAtFull = -Inf;

            end

            if highFullFeasible ...
                    && highUtilityAtFull >= -tol

                highReg = "Fh";
                dH = nH;

            else

                capacityBufferH = ...
                    max(1.0e-8,capacityTolFactor * sH);

                upperBoundH = ...
                    min(nH,sH - capacityBufferH);

                if upperBoundH <= 0

                    error(['Invalid high-tier demand interval at ' ...
                           'eta = %.4f and v = %.4f.'],eta,v);

                end

                highDemandEquation = @(d) ...
                    Vh ...
                    - tH * d / nH ...
                    - ghBench / (sH - d);

                highUtilityAtUpper = ...
                    highDemandEquation(upperBoundH);

                if highUtilityAtUpper >= -tol ...
                        && abs(upperBoundH - nH) <= tol

                    highReg = "Fh";
                    dH = nH;

                elseif highUtilityAtUpper >= 0

                    error(['Unable to bracket the post-merger ' ...
                           'high-tier equilibrium at eta = %.4f ' ...
                           'and v = %.4f.'],eta,v);

                else

                    dH = fzero( ...
                        highDemandEquation, ...
                        [0,upperBoundH]);

                    if dH <= tol

                        highReg = "0h";
                        dH = 0;

                    elseif abs(dH - nH) <= tol

                        highReg = "Fh";
                        dH = nH;

                    else

                        highReg = "Ph";

                    end

                end

            end

        end

        %% =================================================
        % 4.2 Post-merger low-tier equilibrium
        %
        % Marginal-consumer equation:
        %
        % Vl - tL*dL/(2*nL) - glBench/(sL-dL) = 0
        %% =================================================

        lowUtilityAtZero = ...
            Vl - glBench / sL;

        if lowUtilityAtZero <= tol

            lowReg = "0l";
            dL = 0;

        else

            lowFullFeasible = ...
                sL - nL > tol;

            if lowFullFeasible

                lowUtilityAtFull = ...
                    Vl ...
                    - tL / 2 ...
                    - glBench / (sL - nL);

            else

                lowUtilityAtFull = -Inf;

            end

            if lowFullFeasible ...
                    && lowUtilityAtFull >= -tol

                lowReg = "Fl";
                dL = nL;

            else

                capacityBufferL = ...
                    max(1.0e-8,capacityTolFactor * sL);

                upperBoundL = ...
                    min(nL,sL - capacityBufferL);

                if upperBoundL <= 0

                    error(['Invalid low-tier demand interval at ' ...
                           'eta = %.4f and v = %.4f.'],eta,v);

                end

                lowDemandEquation = @(d) ...
                    Vl ...
                    - tL * d / (2 * nL) ...
                    - glBench / (sL - d);

                lowUtilityAtUpper = ...
                    lowDemandEquation(upperBoundL);

                if lowUtilityAtUpper >= -tol ...
                        && abs(upperBoundL - nL) <= tol

                    lowReg = "Fl";
                    dL = nL;

                elseif lowUtilityAtUpper >= 0

                    error(['Unable to bracket the post-merger ' ...
                           'low-tier equilibrium at eta = %.4f ' ...
                           'and v = %.4f.'],eta,v);

                else

                    dL = fzero( ...
                        lowDemandEquation, ...
                        [0,upperBoundL]);

                    if dL <= tol

                        lowReg = "0l";
                        dL = 0;

                    elseif abs(dL - nL) <= tol

                        lowReg = "Fl";
                        dL = nL;

                    else

                        lowReg = "Pl";

                    end

                end

            end

        end

        %% =================================================
        % 4.3 Platform profit
        %% =================================================

        profit = ...
            pH * dH ...
            + pL * dL;

        %% =================================================
        % 4.4 High-tier consumer surplus
        %% =================================================

        if highReg == "0h"

            csh = 0;

        elseif highReg == "Fh"

            csh = ...
                nH ...
                * (Vh ...
                - ghBench / (sH - nH) ...
                - tH / 2);

        else

            csh = ...
                tH * dH^2 ...
                / (2 * nH);

        end

        %% =================================================
        % 4.5 Low-tier consumer surplus
        %% =================================================

        if lowReg == "0l"

            csl = 0;

        elseif lowReg == "Fl"

            csl = ...
                nL ...
                * (Vl ...
                - glBench / (sL - nL) ...
                - tL / 4);

        else

            csl = ...
                tL * dL^2 ...
                / (4 * nL);

        end

        %% Numerical safeguards

        if csh < -tol

            warning(['Negative high-tier consumer surplus at ' ...
                     'eta = %.4f and v = %.4f.'],eta,v);

        end

        if csl < -tol

            warning(['Negative low-tier consumer surplus at ' ...
                     'eta = %.4f and v = %.4f.'],eta,v);

        end

        csh = max(csh,0);
        csl = max(csl,0);

        %% =================================================
        % 4.6 Social welfare
        %% =================================================

        welfare = ...
            profit ...
            + csh ...
            + csl;

        %% =================================================
        % 4.7 Store results
        %% =================================================

        DHAll(e,idx) = dH;
        DLAll(e,idx) = dL;

        ProfitAll(e,idx) = profit;
        CSHAll(e,idx) = csh;
        CSLAll(e,idx) = csl;
        WelfareAll(e,idx) = welfare;

        HighRegimeAll(e,idx) = highReg;
        LowRegimeAll(e,idx) = lowReg;

    end

end

%% 5. Incremental effects relative to eta = 0

DeltaProfit = ...
    ProfitAll(2:end,:) - ProfitAll(1,:);

DeltaCSH = ...
    CSHAll(2:end,:) - CSHAll(1,:);

DeltaCSL = ...
    CSLAll(2:end,:) - CSLAll(1,:);

DeltaWelfare = ...
    WelfareAll(2:end,:) - WelfareAll(1,:);

%% 6. Numerical consistency checks

assert(all(isfinite(ProfitAll(:))), ...
    'ProfitAll contains nonfinite values.');

assert(all(isfinite(CSHAll(:))), ...
    'CSHAll contains nonfinite values.');

assert(all(isfinite(CSLAll(:))), ...
    'CSLAll contains nonfinite values.');

assert(all(isfinite(WelfareAll(:))), ...
    'WelfareAll contains nonfinite values.');

assert(all(DHAll(:) >= -tol), ...
    'A negative high-tier demand was detected.');

assert(all(DLAll(:) >= -tol), ...
    'A negative low-tier demand was detected.');

assert(all(DHAll(:) <= nH + tol), ...
    'High-tier demand exceeds the market size.');

assert(all(DLAll(:) <= nL + tol), ...
    'Low-tier demand exceeds the market size.');

% Proposition 3 predicts weakly decreasing high-tier CS and
% weakly increasing low-tier CS under fixed backfilling.
if any(DeltaCSH(:) > 1.0e-7)

    warning(['Positive high-tier CS increments were detected. ' ...
             'Inspect market-boundary transitions.']);

end

if any(DeltaCSL(:) < -1.0e-7)

    warning(['Negative low-tier CS increments were detected. ' ...
             'Inspect market-boundary transitions.']);

end

%% 7. Display regime transitions

disp(' ');
disp('Market-regime transitions over v:');

for e = 1:nEtaAll

    highChange = ...
        [true, ...
        HighRegimeAll(e,2:end) ~= ...
        HighRegimeAll(e,1:end-1)];

    lowChange = ...
        [true, ...
        LowRegimeAll(e,2:end) ~= ...
        LowRegimeAll(e,1:end-1)];

    fprintf('\neta = %.2f\n',etaAll(e));

    fprintf('  High-tier: ');

    highIndex = find(highChange);

    for k = 1:length(highIndex)

        currentIndex = highIndex(k);

        fprintf('%s at v = %.2f', ...
            char(HighRegimeAll(e,currentIndex)), ...
            vGrid(currentIndex));

        if k < length(highIndex)
            fprintf(' -> ');
        end

    end

    fprintf('\n');

    fprintf('  Low-tier:  ');

    lowIndex = find(lowChange);

    for k = 1:length(lowIndex)

        currentIndex = lowIndex(k);

        fprintf('%s at v = %.2f', ...
            char(LowRegimeAll(e,currentIndex)), ...
            vGrid(currentIndex));

        if k < length(lowIndex)
            fprintf(' -> ');
        end

    end

    fprintf('\n');

end


%% 8. Display sign-transition intervals

disp(' ');
disp('Incremental-effect sign transitions:');

effectNames = { ...
    'Delta Profit', ...
    'Delta Welfare' ...
    };

effectArrays = { ...
    DeltaProfit, ...
    DeltaWelfare ...
    };

signTol = 1.0e-8;

for m = 1:length(effectArrays)

    currentEffect = effectArrays{m};

    fprintf('\n%s\n',effectNames{m});

    for i = 1:nEta

        currentSeries = currentEffect(i,:);

        transitionFound = false;
        previousNonzeroSign = 0;
        previousNonzeroIndex = NaN;

        fprintf('  etaBar = %.2f: ',etaBarList(i));

        for idx = 1:nV

            currentValue = currentSeries(idx);

            if abs(currentValue) <= signTol
                continue;
            end

            currentSign = sign(currentValue);

            if previousNonzeroSign ~= 0 ...
                    && currentSign ~= previousNonzeroSign

                transitionFound = true;

                fprintf('[%.2f, %.2f] ', ...
                    vGrid(previousNonzeroIndex), ...
                    vGrid(idx));

            end

            previousNonzeroSign = currentSign;
            previousNonzeroIndex = idx;

        end

        if ~transitionFound
            fprintf('no strict sign reversal');
        end

        fprintf('\n');

    end

end


%% 9. Plot settings

fontName = 'Times New Roman';

colors = [ ...
    0.0000,0.4470,0.7410;
    0.8500,0.3250,0.0980;
    0.4660,0.6740,0.1880;
    0.4940,0.1840,0.5560 ...
    ];

lineStyles = {'-','--','-.',':'};

legendLabels = cell(1,nEta);

for i = 1:nEta

    legendLabels{i} = ...
        sprintf('$\\bar{\\eta}=%.2f$',etaBarList(i));

end


%% =========================================================
% 10. Figure 6(a): Incremental high-tier consumer surplus
%% =========================================================

fig6a = figure( ...
    'Color','w', ...
    'Position',[100,100,850,620], ...
    'Renderer','painters');

hold on;

for i = 1:nEta

    plot( ...
        vGrid, ...
        DeltaCSH(i,:), ...
        'LineStyle',lineStyles{i}, ...
        'Color',colors(i,:), ...
        'LineWidth',2.2);

end

yline(0,'k:','LineWidth',1.2);

xlabel('$v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',16);

ylabel('$\Delta CS_h$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',16);

set(gca, ...
    'FontName',fontName, ...
    'FontSize',13, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','on', ...
    'Layer','top', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'GridAlpha',0.15);

legend(legendLabels, ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',12, ...
    'Location','southwest', ...
    'Box','on');

text(0.025,0.965,'(a)', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');

xlim([vMin,vMax]);

drawnow;

exportgraphics(fig6a,'Figure6a.pdf','ContentType','vector');
exportgraphics(fig6a,'Figure6a.png','Resolution',600);


%% =========================================================
% 11. Figure 6(b): Incremental low-tier consumer surplus
%% =========================================================

fig6b = figure( ...
    'Color','w', ...
    'Position',[120,120,850,620], ...
    'Renderer','painters');

hold on;

for i = 1:nEta

    plot( ...
        vGrid, ...
        DeltaCSL(i,:), ...
        'LineStyle',lineStyles{i}, ...
        'Color',colors(i,:), ...
        'LineWidth',2.2);

end

yline(0,'k:','LineWidth',1.2);

xlabel('$v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',16);

ylabel('$\Delta CS_\ell$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',16);

set(gca, ...
    'FontName',fontName, ...
    'FontSize',13, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','on', ...
    'Layer','top', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'GridAlpha',0.15);

legend(legendLabels, ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',12, ...
    'Location','southeast', ...
    'Box','on');

text(0.025,0.965,'(b)', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');

xlim([vMin,vMax]);

drawnow;

exportgraphics(fig6b,'Figure6b.pdf','ContentType','vector');
exportgraphics(fig6b,'Figure6b.png','Resolution',600);


%% =========================================================
% 12. Figure 6(c): Incremental platform profit
%% =========================================================

fig6c = figure( ...
    'Color','w', ...
    'Position',[140,140,850,620], ...
    'Renderer','painters');

hold on;

for i = 1:nEta

    plot( ...
        vGrid, ...
        DeltaProfit(i,:), ...
        'LineStyle',lineStyles{i}, ...
        'Color',colors(i,:), ...
        'LineWidth',2.2);

end

yline(0,'k:','LineWidth',1.2);

xlabel('$v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',16);

ylabel('$\Delta \Pi$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',16);

set(gca, ...
    'FontName',fontName, ...
    'FontSize',13, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','on', ...
    'Layer','top', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'GridAlpha',0.15);

legend(legendLabels, ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',12, ...
    'Location','southwest', ...
    'Box','on');

text(0.025,0.965,'(c)', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');

xlim([vMin,vMax]);

drawnow;

exportgraphics(fig6c,'Figure6c.pdf','ContentType','vector');
exportgraphics(fig6c,'Figure6c.png','Resolution',600);


%% =========================================================
% 13. Figure 6(d): Incremental social welfare
%% =========================================================

fig6d = figure( ...
    'Color','w', ...
    'Position',[160,160,850,620], ...
    'Renderer','painters');

hold on;

for i = 1:nEta

    plot( ...
        vGrid, ...
        DeltaWelfare(i,:), ...
        'LineStyle',lineStyles{i}, ...
        'Color',colors(i,:), ...
        'LineWidth',2.2);

end

yline(0,'k:','LineWidth',1.2);

xlabel('$v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',16);

ylabel('$\Delta SW$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',16);

set(gca, ...
    'FontName',fontName, ...
    'FontSize',13, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','on', ...
    'Layer','top', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'GridAlpha',0.15);

legend(legendLabels, ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',12, ...
    'Location','southwest', ...
    'Box','on');

text(0.025,0.965,'(d)', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');

xlim([vMin,vMax]);

drawnow;

exportgraphics(fig6d,'Figure6d.pdf','ContentType','vector');
exportgraphics(fig6d,'Figure6d.png','Resolution',600);
