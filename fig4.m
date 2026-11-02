clear;
clc;
close all;

%% =========================================================
% Figure 4
%
% Comparative statics with respect to the common original
% service value v.
%
% Institutional scenarios:
%   1. Pre-merger
%   2. Post-merger with eta = 0
%   3. Post-merger with eta = etaBar
%
% This version uses tolerance-based boundary identification
% and fzero to avoid false "no feasible equilibrium" errors
% at regime-transition points such as v = 1.6200.
%% =========================================================

%% =========================================================
% 1. Unified calibration
%% =========================================================

nTotal = 100;
alpha = 0.30;

nH = alpha * nTotal;          % High-tier market size: 30
nL = (1 - alpha) * nTotal;    % Low-tier market size: 70

tH = 1;
tL = 1;

%% Capacities before the merger

sHPre = 50;                   % High-tier capacity
sLPre = 50;                   % Low-tier capacity per platform

%% Capacities after the merger

SG = 50;                      % Merged high-tier/GPU pool
SC = 2 * sLPre;               % Merged low-tier/CPU pool: 100

%% Fixed backfilling ratio

etaBar = 0.30;

%% Prices

pH = 2;
pL = 1.5;

%% Original service-value range

vMin = 0;
vMax = 6;
vStep = 0.02;

vGrid = vMin:vStep:vMax;
nV = length(vGrid);

%% Congestion sensitivity

ghBench = 20;
glBench = 12;

%% Task-type-specific service-value adjustments

deltaH = 0;
deltaL = 0;

%% Numerical tolerance

tol = 1.0e-10;

%% =========================================================
% 2. Parameter validation
%% =========================================================

assert(nTotal > 0, ...
    'nTotal must be positive.');

assert(alpha > 0 && alpha < 1, ...
    'alpha must lie strictly between zero and one.');

assert(nH > 0 && nL > 0, ...
    'Both market sizes must be positive.');

assert(tH > 0 && tL > 0, ...
    'Differentiation parameters must be positive.');

assert(sHPre > 0 && sLPre > 0, ...
    'Pre-merger capacities must be positive.');

assert(SG > 0 && SC > 0, ...
    'Post-merger capacities must be positive.');

assert(etaBar >= 0 && etaBar < 1, ...
    'etaBar must lie in [0,1).');

assert(ghBench >= 0 && glBench >= 0, ...
    'Congestion parameters must be nonnegative.');

%% =========================================================
% 3. Preallocate result arrays
%% =========================================================

ProfitPre = zeros(1,nV);
ProfitPost0 = zeros(1,nV);
ProfitPostEta = zeros(1,nV);

CSHPre = zeros(1,nV);
CSHPost0 = zeros(1,nV);
CSHPostEta = zeros(1,nV);

CSLPre = zeros(1,nV);
CSLPost0 = zeros(1,nV);
CSLPostEta = zeros(1,nV);

WelfarePre = zeros(1,nV);
WelfarePost0 = zeros(1,nV);
WelfarePostEta = zeros(1,nV);

%% Store demands

DHPre = zeros(1,nV);
DHPost0 = zeros(1,nV);
DHPostEta = zeros(1,nV);

DLPre = zeros(1,nV);
DLPost0 = zeros(1,nV);
DLPostEta = zeros(1,nV);

%% Store market regimes

HPreRegime = strings(1,nV);
HPost0Regime = strings(1,nV);
HPostEtaRegime = strings(1,nV);

LPreRegime = strings(1,nV);
LPost0Regime = strings(1,nV);
LPostEtaRegime = strings(1,nV);

%% =========================================================
% 4. Main simulation loop
%% =========================================================

for idx = 1:nV

    v = vGrid(idx);

    %% -----------------------------------------------------
    % Gross value v is converted into tier-specific net value
    %% -----------------------------------------------------

    Vh = v - pH + deltaH;
    Vl = v - pL + deltaL;

    %% =====================================================
    % 4.1 Pre-merger high-tier market
    %
    % Marginal-consumer equation:
    %
    % Vh - tH*dH/nH - ghBench/(sHPre-dH) = 0
    %% =====================================================

    sH = sHPre;

    utilityAtZero = ...
        Vh - ghBench / sH;

    if utilityAtZero <= tol

        hPreReg = "0h";
        dHPre = 0;

    else

        fullCoverageFeasible = ...
            sH - nH > tol;

        if fullCoverageFeasible

            utilityAtFull = ...
                Vh ...
                - tH ...
                - ghBench / (sH - nH);

        else

            utilityAtFull = -Inf;

        end

        if fullCoverageFeasible ...
                && utilityAtFull >= -tol

            hPreReg = "Fh";
            dHPre = nH;

        else

            capacityBuffer = ...
                max(1.0e-8,1.0e-8 * sH);

            upperBound = ...
                min(nH,sH - capacityBuffer);

            if upperBound <= 0

                error(['Invalid high-tier demand interval ' ...
                       'before the merger at v = %.4f.'],v);

            end

            highDemandEquation = @(d) ...
                Vh ...
                - tH * d / nH ...
                - ghBench / (sH - d);

            utilityAtUpper = ...
                highDemandEquation(upperBound);

            if utilityAtUpper >= -tol ...
                    && abs(upperBound - nH) <= tol

                hPreReg = "Fh";
                dHPre = nH;

            elseif utilityAtUpper >= 0

                error(['Unable to bracket the pre-merger ' ...
                       'high-tier equilibrium at v = %.4f.'],v);

            else

                dHPre = fzero( ...
                    highDemandEquation, ...
                    [0,upperBound]);

                if dHPre <= tol

                    hPreReg = "0h";
                    dHPre = 0;

                elseif abs(dHPre - nH) <= tol

                    hPreReg = "Fh";
                    dHPre = nH;

                else

                    hPreReg = "Ph";

                end

            end

        end

    end

    %% =====================================================
    % 4.2 Pre-merger low-tier market
    %
    % Total low-tier demand is equally divided between
    % the two platforms. The per-platform load is dL/2.
    %
    % Marginal-consumer equation:
    %
    % Vl - tL*dL/(2*nL)
    %    - glBench/(sLPre-dL/2) = 0
    %% =====================================================

    sL = sLPre;

    utilityAtZero = ...
        Vl - glBench / sL;

    if utilityAtZero <= tol

        lPreReg = "0l";
        dLPre = 0;

    else

        fullCoverageFeasible = ...
            sL - nL / 2 > tol;

        if fullCoverageFeasible

            utilityAtFull = ...
                Vl ...
                - tL / 2 ...
                - glBench / (sL - nL / 2);

        else

            utilityAtFull = -Inf;

        end

        if fullCoverageFeasible ...
                && utilityAtFull >= -tol

            lPreReg = "Fl";
            dLPre = nL;

        else

            capacityBuffer = ...
                max(1.0e-8,1.0e-8 * sL);

            maximumStableDemand = ...
                2 * (sL - capacityBuffer);

            upperBound = ...
                min(nL,maximumStableDemand);

            if upperBound <= 0

                error(['Invalid low-tier demand interval ' ...
                       'before the merger at v = %.4f.'],v);

            end

            lowDemandEquation = @(d) ...
                Vl ...
                - tL * d / (2 * nL) ...
                - glBench / (sL - d / 2);

            utilityAtUpper = ...
                lowDemandEquation(upperBound);

            if utilityAtUpper >= -tol ...
                    && abs(upperBound - nL) <= tol

                lPreReg = "Fl";
                dLPre = nL;

            elseif utilityAtUpper >= 0

                error(['Unable to bracket the pre-merger ' ...
                       'low-tier equilibrium at v = %.4f.'],v);

            else

                dLPre = fzero( ...
                    lowDemandEquation, ...
                    [0,upperBound]);

                if dLPre <= tol

                    lPreReg = "0l";
                    dLPre = 0;

                elseif abs(dLPre - nL) <= tol

                    lPreReg = "Fl";
                    dLPre = nL;

                else

                    lPreReg = "Pl";

                end

            end

        end

    end

    %% =====================================================
    % 4.3 Post-merger high-tier market: eta = 0
    %% =====================================================

    eta = 0;
    sH = (1 - eta) * SG;

    utilityAtZero = ...
        Vh - ghBench / sH;

    if utilityAtZero <= tol

        hPost0Reg = "0h";
        dHPost0 = 0;

    else

        fullCoverageFeasible = ...
            sH - nH > tol;

        if fullCoverageFeasible

            utilityAtFull = ...
                Vh ...
                - tH ...
                - ghBench / (sH - nH);

        else

            utilityAtFull = -Inf;

        end

        if fullCoverageFeasible ...
                && utilityAtFull >= -tol

            hPost0Reg = "Fh";
            dHPost0 = nH;

        else

            capacityBuffer = ...
                max(1.0e-8,1.0e-8 * sH);

            upperBound = ...
                min(nH,sH - capacityBuffer);

            if upperBound <= 0

                error(['Invalid post-merger high-tier demand ' ...
                       'interval at eta = 0 and v = %.4f.'],v);

            end

            highDemandEquation = @(d) ...
                Vh ...
                - tH * d / nH ...
                - ghBench / (sH - d);

            utilityAtUpper = ...
                highDemandEquation(upperBound);

            if utilityAtUpper >= -tol ...
                    && abs(upperBound - nH) <= tol

                hPost0Reg = "Fh";
                dHPost0 = nH;

            elseif utilityAtUpper >= 0

                error(['Unable to bracket the post-merger ' ...
                       'high-tier equilibrium at eta = 0 ' ...
                       'and v = %.4f.'],v);

            else

                dHPost0 = fzero( ...
                    highDemandEquation, ...
                    [0,upperBound]);

                if dHPost0 <= tol

                    hPost0Reg = "0h";
                    dHPost0 = 0;

                elseif abs(dHPost0 - nH) <= tol

                    hPost0Reg = "Fh";
                    dHPost0 = nH;

                else

                    hPost0Reg = "Ph";

                end

            end

        end

    end

    %% =====================================================
    % 4.4 Post-merger low-tier market: eta = 0
    %
    % This block contains the correction for the error at
    % v = 1.6200. At the zero-coverage boundary, demand is
    % set equal to zero instead of being treated as an
    % infeasible partial-coverage equilibrium.
    %% =====================================================

    sL = SC + eta * SG;

    utilityAtZero = ...
        Vl - glBench / sL;

    if utilityAtZero <= tol

        lPost0Reg = "0l";
        dLPost0 = 0;

    else

        fullCoverageFeasible = ...
            sL - nL > tol;

        if fullCoverageFeasible

            utilityAtFull = ...
                Vl ...
                - tL / 2 ...
                - glBench / (sL - nL);

        else

            utilityAtFull = -Inf;

        end

        if fullCoverageFeasible ...
                && utilityAtFull >= -tol

            lPost0Reg = "Fl";
            dLPost0 = nL;

        else

            capacityBuffer = ...
                max(1.0e-8,1.0e-8 * sL);

            upperBound = ...
                min(nL,sL - capacityBuffer);

            if upperBound <= 0

                error(['Invalid post-merger low-tier demand ' ...
                       'interval at eta = 0 and v = %.4f.'],v);

            end

            lowDemandEquation = @(d) ...
                Vl ...
                - tL * d / (2 * nL) ...
                - glBench / (sL - d);

            utilityAtUpper = ...
                lowDemandEquation(upperBound);

            if utilityAtUpper >= -tol ...
                    && abs(upperBound - nL) <= tol

                lPost0Reg = "Fl";
                dLPost0 = nL;

            elseif utilityAtUpper >= 0

                error(['Unable to bracket the post-merger ' ...
                       'low-tier equilibrium at eta = 0 ' ...
                       'and v = %.4f.'],v);

            else

                dLPost0 = fzero( ...
                    lowDemandEquation, ...
                    [0,upperBound]);

                if dLPost0 <= tol

                    lPost0Reg = "0l";
                    dLPost0 = 0;

                elseif abs(dLPost0 - nL) <= tol

                    lPost0Reg = "Fl";
                    dLPost0 = nL;

                else

                    lPost0Reg = "Pl";

                end

            end

        end

    end

    %% =====================================================
    % 4.5 Post-merger high-tier market: eta = etaBar
    %% =====================================================

    eta = etaBar;
    sH = (1 - eta) * SG;

    utilityAtZero = ...
        Vh - ghBench / sH;

    if utilityAtZero <= tol

        hPostEtaReg = "0h";
        dHPostEta = 0;

    else

        fullCoverageFeasible = ...
            sH - nH > tol;

        if fullCoverageFeasible

            utilityAtFull = ...
                Vh ...
                - tH ...
                - ghBench / (sH - nH);

        else

            utilityAtFull = -Inf;

        end

        if fullCoverageFeasible ...
                && utilityAtFull >= -tol

            hPostEtaReg = "Fh";
            dHPostEta = nH;

        else

            capacityBuffer = ...
                max(1.0e-8,1.0e-8 * sH);

            upperBound = ...
                min(nH,sH - capacityBuffer);

            if upperBound <= 0

                error(['Invalid post-merger high-tier demand ' ...
                       'interval at eta = etaBar and v = %.4f.'],v);

            end

            highDemandEquation = @(d) ...
                Vh ...
                - tH * d / nH ...
                - ghBench / (sH - d);

            utilityAtUpper = ...
                highDemandEquation(upperBound);

            if utilityAtUpper >= -tol ...
                    && abs(upperBound - nH) <= tol

                hPostEtaReg = "Fh";
                dHPostEta = nH;

            elseif utilityAtUpper >= 0

                error(['Unable to bracket the post-merger ' ...
                       'high-tier equilibrium at eta = etaBar ' ...
                       'and v = %.4f.'],v);

            else

                dHPostEta = fzero( ...
                    highDemandEquation, ...
                    [0,upperBound]);

                if dHPostEta <= tol

                    hPostEtaReg = "0h";
                    dHPostEta = 0;

                elseif abs(dHPostEta - nH) <= tol

                    hPostEtaReg = "Fh";
                    dHPostEta = nH;

                else

                    hPostEtaReg = "Ph";

                end

            end

        end

    end

    %% =====================================================
    % 4.6 Post-merger low-tier market: eta = etaBar
    %% =====================================================

    sL = SC + eta * SG;

    utilityAtZero = ...
        Vl - glBench / sL;

    if utilityAtZero <= tol

        lPostEtaReg = "0l";
        dLPostEta = 0;

    else

        fullCoverageFeasible = ...
            sL - nL > tol;

        if fullCoverageFeasible

            utilityAtFull = ...
                Vl ...
                - tL / 2 ...
                - glBench / (sL - nL);

        else

            utilityAtFull = -Inf;

        end

        if fullCoverageFeasible ...
                && utilityAtFull >= -tol

            lPostEtaReg = "Fl";
            dLPostEta = nL;

        else

            capacityBuffer = ...
                max(1.0e-8,1.0e-8 * sL);

            upperBound = ...
                min(nL,sL - capacityBuffer);

            if upperBound <= 0

                error(['Invalid post-merger low-tier demand ' ...
                       'interval at eta = etaBar and v = %.4f.'],v);

            end

            lowDemandEquation = @(d) ...
                Vl ...
                - tL * d / (2 * nL) ...
                - glBench / (sL - d);

            utilityAtUpper = ...
                lowDemandEquation(upperBound);

            if utilityAtUpper >= -tol ...
                    && abs(upperBound - nL) <= tol

                lPostEtaReg = "Fl";
                dLPostEta = nL;

            elseif utilityAtUpper >= 0

                error(['Unable to bracket the post-merger ' ...
                       'low-tier equilibrium at eta = etaBar ' ...
                       'and v = %.4f.'],v);

            else

                dLPostEta = fzero( ...
                    lowDemandEquation, ...
                    [0,upperBound]);

                if dLPostEta <= tol

                    lPostEtaReg = "0l";
                    dLPostEta = 0;

                elseif abs(dLPostEta - nL) <= tol

                    lPostEtaReg = "Fl";
                    dLPostEta = nL;

                else

                    lPostEtaReg = "Pl";

                end

            end

        end

    end

    %% =====================================================
    % 4.7 Platform profit
    %% =====================================================

    profitPre = ...
        pH * dHPre ...
        + pL * dLPre;

    profitPost0 = ...
        pH * dHPost0 ...
        + pL * dLPost0;

    profitPostEta = ...
        pH * dHPostEta ...
        + pL * dLPostEta;

    %% =====================================================
    % 4.8 High-tier consumer surplus
    %% =====================================================

    if hPreReg == "0h"

        cshPre = 0;

    elseif hPreReg == "Fh"

        cshPre = ...
            nH ...
            * (Vh ...
            - ghBench / (sHPre - nH) ...
            - tH / 2);

    else

        cshPre = ...
            tH * dHPre^2 ...
            / (2 * nH);

    end

    if hPost0Reg == "0h"

        cshPost0 = 0;

    elseif hPost0Reg == "Fh"

        sH0 = SG;

        cshPost0 = ...
            nH ...
            * (Vh ...
            - ghBench / (sH0 - nH) ...
            - tH / 2);

    else

        cshPost0 = ...
            tH * dHPost0^2 ...
            / (2 * nH);

    end

    if hPostEtaReg == "0h"

        cshPostEta = 0;

    elseif hPostEtaReg == "Fh"

        sHEta = ...
            (1 - etaBar) * SG;

        cshPostEta = ...
            nH ...
            * (Vh ...
            - ghBench / (sHEta - nH) ...
            - tH / 2);

    else

        cshPostEta = ...
            tH * dHPostEta^2 ...
            / (2 * nH);

    end

    %% =====================================================
    % 4.9 Low-tier consumer surplus
    %% =====================================================

    if lPreReg == "0l"

        cslPre = 0;

    elseif lPreReg == "Fl"

        cslPre = ...
            nL ...
            * (Vl ...
            - glBench / (sLPre - nL / 2) ...
            - tL / 4);

    else

        cslPre = ...
            tL * dLPre^2 ...
            / (4 * nL);

    end

    if lPost0Reg == "0l"

        cslPost0 = 0;

    elseif lPost0Reg == "Fl"

        sL0 = SC;

        cslPost0 = ...
            nL ...
            * (Vl ...
            - glBench / (sL0 - nL) ...
            - tL / 4);

    else

        cslPost0 = ...
            tL * dLPost0^2 ...
            / (4 * nL);

    end

    if lPostEtaReg == "0l"

        cslPostEta = 0;

    elseif lPostEtaReg == "Fl"

        sLEta = ...
            SC + etaBar * SG;

        cslPostEta = ...
            nL ...
            * (Vl ...
            - glBench / (sLEta - nL) ...
            - tL / 4);

    else

        cslPostEta = ...
            tL * dLPostEta^2 ...
            / (4 * nL);

    end

    %% Numerical safeguards

    cshPre = max(cshPre,0);
    cshPost0 = max(cshPost0,0);
    cshPostEta = max(cshPostEta,0);

    cslPre = max(cslPre,0);
    cslPost0 = max(cslPost0,0);
    cslPostEta = max(cslPostEta,0);

    %% =====================================================
    % 4.10 Total social welfare
    %% =====================================================

    welfarePre = ...
        profitPre ...
        + cshPre ...
        + cslPre;

    welfarePost0 = ...
        profitPost0 ...
        + cshPost0 ...
        + cslPost0;

    welfarePostEta = ...
        profitPostEta ...
        + cshPostEta ...
        + cslPostEta;

    %% =====================================================
    % 4.11 Store results
    %% =====================================================

    DHPre(idx) = dHPre;
    DHPost0(idx) = dHPost0;
    DHPostEta(idx) = dHPostEta;

    DLPre(idx) = dLPre;
    DLPost0(idx) = dLPost0;
    DLPostEta(idx) = dLPostEta;

    ProfitPre(idx) = profitPre;
    ProfitPost0(idx) = profitPost0;
    ProfitPostEta(idx) = profitPostEta;

    CSHPre(idx) = cshPre;
    CSHPost0(idx) = cshPost0;
    CSHPostEta(idx) = cshPostEta;

    CSLPre(idx) = cslPre;
    CSLPost0(idx) = cslPost0;
    CSLPostEta(idx) = cslPostEta;

    WelfarePre(idx) = welfarePre;
    WelfarePost0(idx) = welfarePost0;
    WelfarePostEta(idx) = welfarePostEta;

    HPreRegime(idx) = hPreReg;
    HPost0Regime(idx) = hPost0Reg;
    HPostEtaRegime(idx) = hPostEtaReg;

    LPreRegime(idx) = lPreReg;
    LPost0Regime(idx) = lPost0Reg;
    LPostEtaRegime(idx) = lPostEtaReg;

end

%% =========================================================
% 5. Separate merger and backfilling effects
%
% Merger effect:
%   Post-merger eta=0 minus pre-merger
%
% Backfilling effect:
%   Post-merger eta=etaBar minus post-merger eta=0
%% =========================================================

DeltaMergerProfit = ...
    ProfitPost0 - ProfitPre;

DeltaBackfillProfit = ...
    ProfitPostEta - ProfitPost0;

DeltaMergerCSH = ...
    CSHPost0 - CSHPre;

DeltaBackfillCSH = ...
    CSHPostEta - CSHPost0;

DeltaMergerCSL = ...
    CSLPost0 - CSLPre;

DeltaBackfillCSL = ...
    CSLPostEta - CSLPost0;

DeltaMergerWelfare = ...
    WelfarePost0 - WelfarePre;

DeltaBackfillWelfare = ...
    WelfarePostEta - WelfarePost0;

%% =========================================================
% 6. Display regime transitions
%% =========================================================

disp(' ');
disp('Market-regime transitions over v:');

regimeNames = { ...
    'Pre-merger high tier', ...
    'Pre-merger low tier', ...
    'Post-merger high tier, eta=0', ...
    'Post-merger low tier, eta=0', ...
    'Post-merger high tier, eta=etaBar', ...
    'Post-merger low tier, eta=etaBar' ...
    };

regimeArrays = { ...
    HPreRegime, ...
    LPreRegime, ...
    HPost0Regime, ...
    LPost0Regime, ...
    HPostEtaRegime, ...
    LPostEtaRegime ...
    };

for r = 1:length(regimeArrays)

    currentRegime = regimeArrays{r};

    changeIndex = ...
        [true, ...
        currentRegime(2:end) ~= ...
        currentRegime(1:end-1)];

    transitionLocations = ...
        find(changeIndex);

    fprintf('\n%s: ',regimeNames{r});

    for k = 1:length(transitionLocations)

        currentIndex = ...
            transitionLocations(k);

        fprintf('%s at v = %.2f', ...
            char(currentRegime(currentIndex)), ...
            vGrid(currentIndex));

        if k < length(transitionLocations)
            fprintf(' -> ');
        end

    end

    fprintf('\n');

end

%% =========================================================
% 7. Display signs of merger and backfilling effects
%% =========================================================

disp(' ');
disp('Effect definitions:');
disp('Delta merger = Post(eta=0) - Pre');
disp('Delta backfilling = Post(eta=etaBar) - Post(eta=0)');

fprintf('\nProfit merger effect range:     [%.4f, %.4f]\n', ...
    min(DeltaMergerProfit), ...
    max(DeltaMergerProfit));

fprintf('Profit backfilling effect range: [%.4f, %.4f]\n', ...
    min(DeltaBackfillProfit), ...
    max(DeltaBackfillProfit));

fprintf('Welfare merger effect range:     [%.4f, %.4f]\n', ...
    min(DeltaMergerWelfare), ...
    max(DeltaMergerWelfare));

fprintf('Welfare backfilling effect range: [%.4f, %.4f]\n', ...
    min(DeltaBackfillWelfare), ...
    max(DeltaBackfillWelfare));

%% =========================================================
% 8. Plot settings
%% =========================================================

fontName = 'Times New Roman';

lineWidth = 2.2;

colorPre = [0.0000,0.0000,0.0000];
colorPost0 = [0.0000,0.4470,0.7410];
colorPostEta = [0.8500,0.3250,0.0980];

legendText = { ...
    '$\mathrm{Pre}$', ...
    '$\mathrm{Post}\;(\eta=0)$', ...
    '$\mathrm{Post}\;(\eta=\bar{\eta})$' ...
    };

%% =========================================================
% 9. Figure 5(a): Total platform profit
%% =========================================================

fig5a = figure( ...
    'Color','w', ...
    'Position',[100,100,850,620], ...
    'Renderer','painters');

hold on;

plot( ...
    vGrid, ...
    ProfitPre, ...
    '-', ...
    'Color',colorPre, ...
    'LineWidth',lineWidth);

plot( ...
    vGrid, ...
    ProfitPost0, ...
    '-', ...
    'Color',colorPost0, ...
    'LineWidth',lineWidth);

plot( ...
    vGrid, ...
    ProfitPostEta, ...
    '--', ...
    'Color',colorPostEta, ...
    'LineWidth',lineWidth);

xlabel( ...
    'Common original service value  $v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',15);

ylabel( ...
    'Total platform profit', ...
    'FontName',fontName, ...
    'FontSize',15);

set( ...
    gca, ...
    'FontName',fontName, ...
    'FontSize',13, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','on', ...
    'Layer','top', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'GridAlpha',0.15);

legend( ...
    legendText, ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',12, ...
    'Location','northwest', ...
    'Box','on');

text( ...
    0.025, ...
    0.965, ...
    '(a)', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');

xlim([vMin,vMax]);

drawnow;

exportgraphics( ...
    fig5a, ...
    'Figure5a.pdf', ...
    'ContentType','vector');

exportgraphics( ...
    fig5a, ...
    'Figure5a.png', ...
    'Resolution',600);

%% =========================================================
% 10. Figure 5(b): High-tier consumer surplus
%% =========================================================

fig5b = figure( ...
    'Color','w', ...
    'Position',[120,120,850,620], ...
    'Renderer','painters');

hold on;

plot( ...
    vGrid, ...
    CSHPre, ...
    '-', ...
    'Color',colorPre, ...
    'LineWidth',lineWidth);

plot( ...
    vGrid, ...
    CSHPost0, ...
    '-', ...
    'Color',colorPost0, ...
    'LineWidth',lineWidth);

plot( ...
    vGrid, ...
    CSHPostEta, ...
    '--', ...
    'Color',colorPostEta, ...
    'LineWidth',lineWidth);

xlabel( ...
    'Common original service value  $v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',15);

ylabel( ...
    'High-tier consumer surplus', ...
    'FontName',fontName, ...
    'FontSize',15);

set( ...
    gca, ...
    'FontName',fontName, ...
    'FontSize',13, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','on', ...
    'Layer','top', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'GridAlpha',0.15);

legend( ...
    legendText, ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',12, ...
    'Location','northwest', ...
    'Box','on');

text( ...
    0.025, ...
    0.965, ...
    '(b)', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');

xlim([vMin,vMax]);

drawnow;

exportgraphics( ...
    fig5b, ...
    'Figure5b.pdf', ...
    'ContentType','vector');

exportgraphics( ...
    fig5b, ...
    'Figure5b.png', ...
    'Resolution',600);

%% =========================================================
% 11. Figure 5(c): Low-tier consumer surplus
%% =========================================================

fig5c = figure( ...
    'Color','w', ...
    'Position',[140,140,850,620], ...
    'Renderer','painters');

hold on;

plot( ...
    vGrid, ...
    CSLPre, ...
    '-', ...
    'Color',colorPre, ...
    'LineWidth',lineWidth);

plot( ...
    vGrid, ...
    CSLPost0, ...
    '-', ...
    'Color',colorPost0, ...
    'LineWidth',lineWidth);

plot( ...
    vGrid, ...
    CSLPostEta, ...
    '--', ...
    'Color',colorPostEta, ...
    'LineWidth',lineWidth);

xlabel( ...
    'Common original service value  $v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',15);

ylabel( ...
    'Low-tier consumer surplus', ...
    'FontName',fontName, ...
    'FontSize',15);

set( ...
    gca, ...
    'FontName',fontName, ...
    'FontSize',13, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','on', ...
    'Layer','top', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'GridAlpha',0.15);

legend( ...
    legendText, ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',12, ...
    'Location','northwest', ...
    'Box','on');

text( ...
    0.025, ...
    0.965, ...
    '(c)', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');

xlim([vMin,vMax]);

drawnow;

exportgraphics( ...
    fig5c, ...
    'Figure5c.pdf', ...
    'ContentType','vector');

exportgraphics( ...
    fig5c, ...
    'Figure5c.png', ...
    'Resolution',600);

%% =========================================================
% 12. Figure 5(d): Total social welfare
%% =========================================================

fig5d = figure( ...
    'Color','w', ...
    'Position',[160,160,850,620], ...
    'Renderer','painters');

hold on;

plot( ...
    vGrid, ...
    WelfarePre, ...
    '-', ...
    'Color',colorPre, ...
    'LineWidth',lineWidth);

plot( ...
    vGrid, ...
    WelfarePost0, ...
    '-', ...
    'Color',colorPost0, ...
    'LineWidth',lineWidth);

plot( ...
    vGrid, ...
    WelfarePostEta, ...
    '--', ...
    'Color',colorPostEta, ...
    'LineWidth',lineWidth);

xlabel( ...
    'Common original service value  $v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',15);

ylabel( ...
    'Total social welfare', ...
    'FontName',fontName, ...
    'FontSize',15);

set( ...
    gca, ...
    'FontName',fontName, ...
    'FontSize',13, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','on', ...
    'Layer','top', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'GridAlpha',0.15);

legend( ...
    legendText, ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',12, ...
    'Location','northwest', ...
    'Box','on');

text( ...
    0.025, ...
    0.965, ...
    '(d)', ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');

xlim([vMin,vMax]);

drawnow;

exportgraphics( ...
    fig5d, ...
    'Figure5d.pdf', ...
    'ContentType','vector');

exportgraphics( ...
    fig5d, ...
    'Figure5d.png', ...
    'Resolution',600);