clear; clc; close all;

%% =========================================================
% Figure 7: Endogenous backfilling in Extended Model 2
%
% This is a completely self-contained MATLAB script.
% It does not call any external or local function.
%
% Decision variable:
%   q = actual capacity transferred from the high-tier pool
%       to the low-tier pool.
%
% Feasible set:
%   0 <= q <= qBar(v),
%   qBar(v) = max{q: SG-q-DH(q;v) >= bH}.
% The common SLA reserve is bH = 1/WbarH. If the allocation at q=0
% already violates the SLA, that parameter point is reported separately
% as "SLA infeasible" and is excluded from endogenous optimization.
%
% Reported endogenous intensity:
%   etaStar = qStar / SG.
%
% The script produces one publication-oriented three-panel figure:
%   (a) welfare-optimal backfilling regime;
%   (b) welfare-optimal actual backfilling amount qSWStar;
%   (c) profit-welfare gap qPiStar-qSWStar.
%
% IMPORTANT:
%   scenarioMode = 'benchmark' reproduces the manuscript benchmark.
%   scenarioMode = 'low_cost_sensitivity' is an explicitly labelled
%   sensitivity case in which a full-feasible-backfilling region can
%   arise economically. The classification colors are never altered
%   merely to manufacture a full-backfilling region.
%% =========================================================

%% 1. Unified calibration used in the revised manuscript

N = 100;
alpha = 0.30;

NH = alpha * N;          % 30
NL = (1 - alpha) * N;    % 70; do not multiply by 2 again

tH = 1;
tL = 1;

pH = 2;
pL = 1.5;

gammaH = 20;
gammaL = 12;

SC = 100;

% Common high-tier SLA calibration used in Extended Models 1 and 2
WbarH = 0.10;
bH = 1 / WbarH;         % required high-tier spare-capacity reserve

% Benchmark endogenous-backfilling cost and externality parameters
kappaBenchmark = 0.10;
tau0Benchmark = 0.35;
tau1Benchmark = 0.25;
omegaBenchmark = 0.10;

% ---------------------------------------------------------
% Scenario selection
% ---------------------------------------------------------
% 'benchmark'
%   Uses the manuscript benchmark cost coefficients. Under this calibration,
%   the absence of a full-feasible region is a numerical result.
%
% 'low_cost_sensitivity'
%   Retains every market parameter above, but scales all endogenous-
%   backfilling cost/externality coefficients by costScale = 0.05.
%   This is the recommended transparent sensitivity scenario for
%   displaying the full-feasible corner predicted by Proposition 8.

scenarioMode = 'benchmark';

switch lower(scenarioMode)

    case 'benchmark'

        costScale = 1.00;
        scenarioLabel = 'Benchmark calibration';
        scenarioFileTag = 'Benchmark';

    case 'low_cost_sensitivity'

        costScale = 0.05;
        scenarioLabel = 'Low-cost sensitivity, $s_c=0.05$';
        scenarioFileTag = 'LowCostSensitivity';

    otherwise

        error(['Unknown scenarioMode. Use ''benchmark'' or ' ...
            '''low_cost_sensitivity''.']);

end

kappa = costScale * kappaBenchmark;
tau0 = costScale * tau0Benchmark;
tau1 = costScale * tau1Benchmark;
omega = costScale * omegaBenchmark;

%% 2. Parameter grids

% Common gross service value in the revised manuscript
vMin = 0;
vMax = 6;
% A finer service-value grid reduces staircase artifacts at regime
% boundaries while retaining the manuscript ratio calibration.
vStep = 0.02;
vGrid = vMin:vStep:vMax;

% Relative abundance of high-tier capacity
ratioMin = 0.20;
ratioMax = 2.20;
ratioStep = 0.02;
ratioGrid = ratioMin:ratioStep:ratioMax;

nV = numel(vGrid);
nRatio = numel(ratioGrid);

%% 3. Numerical-search settings

% Adaptive grid search. The first pass scans the full feasible interval.
% Each subsequent pass refines neighborhoods around both optima.
nCoarse = 401;
nFine = 121;
nRefinePasses = 3;

tolValue = 1.0e-10;
tolCapacity = 1.0e-10;
tolMarginal = 1.0e-8;
tolSLA = 1.0e-9;
nBisection = 80;

%% 4. Storage

qBarAll = nan(nRatio,nV);
BaselineSLAFeasible = false(nRatio,nV);
SLAResidualAtQBar = nan(nRatio,nV);
SLASlackAtPi = nan(nRatio,nV);
SLASlackAtSW = nan(nRatio,nV);

qPiStar = nan(nRatio,nV);
qSWStar = nan(nRatio,nV);

etaPiStar = nan(nRatio,nV);
etaSWStar = nan(nRatio,nV);
etaBarFeasible = nan(nRatio,nV);

PiStar = nan(nRatio,nV);
WStar = nan(nRatio,nV);

% Welfare-optimal regimes:
% 0 = baseline allocation is SLA infeasible
% 1 = no backfilling
% 2 = partial backfilling
% 3 = full feasible backfilling
RegimeSW = zeros(nRatio,nV);

% Modified Proposition 9 diagnostics
Prop9ConditionValue = nan(nRatio,nV);
Prop9RegularPoint = false(nRatio,nV);
Prop9Condition = false(nRatio,nV);
Prop9Ordering = false(nRatio,nV);
Prop9Violation = nan(nRatio,nV);

%% =========================================================
% 5. Main simulation
%% =========================================================

for iRatio = 1:nRatio

    ratio = ratioGrid(iRatio);
    SG = ratio * SC;

    for iV = 1:nV

        v = vGrid(iV);

        Vh = v - pH;
        Vl = v - pL;

        tauV = tau0 + tau1 * v;

        %% -------------------------------------------------
        % 5.1 High-tier demand without backfilling, DH0(v)
        %% -------------------------------------------------

        SH0 = SG;

        if SH0 <= tolCapacity || ...
                Vh <= tolValue || ...
                gammaH >= SH0 * Vh - tolValue

            DH0 = 0;

        elseif SH0 > NH + tolCapacity && ...
                Vh > tH + tolValue && ...
                gammaH <= (SH0 - NH) * (Vh - tH) + tolValue

            DH0 = NH;

        else

            Ah = NH / tH;

            % Corrected partial-coverage demand root:
            % D = [S + A*V - sqrt((S-A*V)^2+4*A*gamma)]/2
            discH0 = (SH0 - Ah * Vh)^2 + 4 * Ah * gammaH;

            DH0 = (SH0 + Ah * Vh - sqrt(max(discH0,0))) / 2;

            DH0 = max(0,min(NH,DH0));
            DH0 = min(DH0,max(SH0 - tolCapacity,0));

        end

        %% -------------------------------------------------
        % 5.2 SLA-feasible upper bound
        %
        % The boundary must be solved implicitly because high-tier
        % demand changes with q. The shortcut SG-DH0-bH is generally
        % incorrect. Define g(q)=SG-q-DH(q;v)-bH and locate the largest
        % q for which g(q)>=0.
        %% -------------------------------------------------

        reserveAtZero = SG - DH0;
        gAtZero = reserveAtZero - bH;

        if gAtZero < -tolSLA

            % Even the no-backfilling allocation violates the common
            % high-tier SLA. Keep NaNs in all endogenous outcomes and
            % retain regime code 0 for this parameter point.
            continue;

        end

        BaselineSLAFeasible(iRatio,iV) = true;

        qUpper = max(SG - bH,0);

        if qUpper <= tolCapacity

            qBar = 0;
            gAtQBar = gAtZero;

        else

            qLow = 0;
            qHigh = qUpper;
            gLow = gAtZero;

            for iBisect = 1:nBisection

                qMid = 0.5 * (qLow + qHigh);
                SHMid = SG - qMid;

                if SHMid <= tolCapacity || ...
                        Vh <= tolValue || ...
                        gammaH >= SHMid * Vh - tolValue

                    DHMid = 0;

                elseif SHMid > NH + tolCapacity && ...
                        Vh > tH + tolValue && ...
                        gammaH <= (SHMid - NH) * (Vh - tH) + tolValue

                    DHMid = NH;

                else

                    Ah = NH / tH;
                    discHMid = (SHMid - Ah * Vh)^2 + 4 * Ah * gammaH;
                    DHMid = (SHMid + Ah * Vh ...
                        - sqrt(max(discHMid,0))) / 2;
                    DHMid = max(0,min(NH,DHMid));
                    DHMid = min(DHMid,max(SHMid - tolCapacity,0));

                end

                gMid = SG - qMid - DHMid - bH;

                if gMid >= 0
                    qLow = qMid;
                    gLow = gMid;
                else
                    qHigh = qMid;
                end

                if qHigh - qLow <= tolSLA
                    break;
                end

            end

            % Use the feasible-side endpoint to prevent a floating-point
            % overshoot of the SLA boundary.
            qBar = qLow;
            gAtQBar = gLow;

        end

        qBarAll(iRatio,iV) = qBar;
        etaBarFeasible(iRatio,iV) = qBar / SG;
        SLAResidualAtQBar(iRatio,iV) = gAtQBar;

        %% -------------------------------------------------
        % 5.3 Adaptive maximization of profit and welfare
        %% -------------------------------------------------

        if qBar <= tolCapacity

            qLeftPi = 0;
            qRightPi = 0;

            qLeftSW = 0;
            qRightSW = 0;

        else

            qLeftPi = 0;
            qRightPi = qBar;

            qLeftSW = 0;
            qRightSW = qBar;

        end

        % At pass 0, scan the entire interval. Later passes refine
        % neighborhoods around the profit and welfare maximizers.
        for refinePass = 0:nRefinePasses

            if qBar <= tolCapacity

                qEval = 0;

            elseif refinePass == 0

                qEval = linspace(0,qBar,nCoarse);

            else

                qEvalPi = linspace(qLeftPi,qRightPi,nFine);
                qEvalSW = linspace(qLeftSW,qRightSW,nFine);

                % Both endpoints are always retained, so corner solutions
                % are compared directly with all interior candidates.
                qEval = unique([0,qEvalPi,qEvalSW,qBar]);

            end

            SH = SG - qEval;
            SL = SC + qEval;

            nQEval = numel(qEval);

            DH = zeros(1,nQEval);
            DL = zeros(1,nQEval);

            CSH = zeros(1,nQEval);
            CSL = zeros(1,nQEval);

            % Marginal consumer surplus with respect to effective capacity.
            % These values are needed for the sufficient condition in the
            % modified Proposition 9.
            CSHs = zeros(1,nQEval);
            CSLs = zeros(1,nQEval);

            %% ---------------------------------------------
            % High-tier market for every candidate q
            %% ---------------------------------------------

            emptyH = (SH <= tolCapacity) | ...
                (Vh <= tolValue) | ...
                (gammaH >= SH .* Vh - tolValue);

            fullH = (~emptyH) & ...
                (SH > NH + tolCapacity) & ...
                (Vh > tH + tolValue) & ...
                (gammaH <= (SH - NH) .* (Vh - tH) + tolValue);

            partialH = (~emptyH) & (~fullH);

            if any(partialH)

                Ah = NH / tH;

                SHp = SH(partialH);

                discH = (SHp - Ah * Vh).^2 + 4 * Ah * gammaH;

                DHp = (SHp + Ah * Vh - sqrt(max(discH,0))) / 2;

                DHp = max(DHp,0);
                DHp = min(DHp,NH);
                DHp = min(DHp,max(SHp - tolCapacity,0));

                DH(partialH) = DHp;

                CSH(partialH) = tH .* DHp.^2 ./ (2 * NH);

                uH = SHp - DHp;

                dDHdSH = Ah * gammaH ./ ...
                    (uH.^2 + Ah * gammaH);

                CSHs(partialH) = ...
                    (tH .* DHp ./ NH) .* dDHdSH;

            end

            if any(fullH)

                SHf = SH(fullH);

                DH(fullH) = NH;

                CSH(fullH) = NH .* (...
                    Vh ...
                    - gammaH ./ (SHf - NH) ...
                    - tH / 2);

                CSH(fullH) = max(CSH(fullH),0);

                CSHs(fullH) = ...
                    NH * gammaH ./ (SHf - NH).^2;

            end

            %% ---------------------------------------------
            % Low-tier market for every candidate q
            %% ---------------------------------------------

            emptyL = (SL <= tolCapacity) | ...
                (Vl <= tolValue) | ...
                (gammaL >= SL .* Vl - tolValue);

            fullL = (~emptyL) & ...
                (SL > NL + tolCapacity) & ...
                (Vl > tL / 2 + tolValue) & ...
                (gammaL <= (SL - NL) .* (Vl - tL / 2) + tolValue);

            partialL = (~emptyL) & (~fullL);

            if any(partialL)

                Al = 2 * NL / tL;

                SLp = SL(partialL);

                discL = (SLp - Al * Vl).^2 + 4 * Al * gammaL;

                DLp = (SLp + Al * Vl - sqrt(max(discL,0))) / 2;

                DLp = max(DLp,0);
                DLp = min(DLp,NL);
                DLp = min(DLp,max(SLp - tolCapacity,0));

                DL(partialL) = DLp;

                CSL(partialL) = tL .* DLp.^2 ./ (4 * NL);

                uL = SLp - DLp;

                dDLdSL = Al * gammaL ./ ...
                    (uL.^2 + Al * gammaL);

                CSLs(partialL) = ...
                    (tL .* DLp ./ (2 * NL)) .* dDLdSL;

            end

            if any(fullL)

                SLf = SL(fullL);

                DL(fullL) = NL;

                CSL(fullL) = NL .* (...
                    Vl ...
                    - gammaL ./ (SLf - NL) ...
                    - tL / 4);

                CSL(fullL) = max(CSL(fullL),0);

                CSLs(fullL) = ...
                    NL * gammaL ./ (SLf - NL).^2;

            end

            %% ---------------------------------------------
            % Profit and welfare objectives
            %% ---------------------------------------------

            revenue = pH .* DH + pL .* DL;

            PiEval = revenue - 0.5 * kappa .* qEval.^2;

            WEval = PiEval + CSH + CSL ...
                - tauV .* qEval ...
                - 0.5 * omega * v .* qEval.^2;

            [PiBestPass,idxPi] = max(PiEval);
            [WBestPass,idxSW] = max(WEval);

            qPiBestPass = qEval(idxPi);
            qSWBestPass = qEval(idxSW);

            % Neighborhood for the next profit refinement
            idxPiLeft = max(idxPi - 1,1);
            idxPiRight = min(idxPi + 1,nQEval);

            qLeftPi = qEval(idxPiLeft);
            qRightPi = qEval(idxPiRight);

            % Neighborhood for the next welfare refinement
            idxSWLeft = max(idxSW - 1,1);
            idxSWRight = min(idxSW + 1,nQEval);

            qLeftSW = qEval(idxSWLeft);
            qRightSW = qEval(idxSWRight);

        end

        %% -------------------------------------------------
        % 5.4 Store final optima
        %% -------------------------------------------------

        qPi = qPiBestPass;
        qSW = qSWBestPass;

        qPiStar(iRatio,iV) = qPi;
        qSWStar(iRatio,iV) = qSW;

        etaPiStar(iRatio,iV) = qPi / SG;
        etaSWStar(iRatio,iV) = qSW / SG;

        PiStar(iRatio,iV) = PiBestPass;
        WStar(iRatio,iV) = WBestPass;

        % Directly verify that both reported optima satisfy the same SLA.
        DHAtPi = DH(idxPi);
        DHAtSW = DH(idxSW);

        SLASlackAtPi(iRatio,iV) = SG - qPi - DHAtPi - bH;
        SLASlackAtSW(iRatio,iV) = SG - qSW - DHAtSW - bH;

        if SLASlackAtPi(iRatio,iV) < -10 * tolSLA || ...
                SLASlackAtSW(iRatio,iV) < -10 * tolSLA
            error('An endogenous optimum violates the high-tier SLA.');
        end

        %% -------------------------------------------------
        % 5.5 Welfare-optimal regime: Proposition 8
        %% -------------------------------------------------

        tolQ = max(1.0e-7,1.0e-6 * max(1,qBar));

        if qSW <= tolQ

            RegimeSW(iRatio,iV) = 1;

        elseif qBar - qSW <= tolQ

            RegimeSW(iRatio,iV) = 3;

        else

            RegimeSW(iRatio,iV) = 2;

        end

        %% -------------------------------------------------
        % 5.6 Modified Proposition 9
        %
        % Sufficient condition evaluated at qPiStar:
        %
        % CSH_s + tau(v) + omega*v*qPiStar >= CSL_s.
        %
        % Under the assumptions of Proposition 9, this condition
        % implies qSWStar <= qPiStar.
        %% -------------------------------------------------

        CSHsAtPi = CSHs(idxPi);
        CSLsAtPi = CSLs(idxPi);

        % Proposition 9 uses a differentiable interior first-order
        % condition. Do not evaluate that sufficient condition exactly
        % at q = 0, q = qBar, or at a 0/P/F market-state kink.
        stateH = zeros(1,nQEval);
        stateH(partialH) = 2;
        stateH(fullH) = 3;

        stateL = zeros(1,nQEval);
        stateL(partialL) = 2;
        stateL(fullL) = 3;

        idxPiLeft = max(idxPi - 1,1);
        idxPiRight = min(idxPi + 1,nQEval);

        profitOptimumIsInterior = ...
            (qPi > tolQ) && ...
            (qBar - qPi > tolQ);

        sameLocalMarketState = ...
            (stateH(idxPiLeft) == stateH(idxPi)) && ...
            (stateH(idxPiRight) == stateH(idxPi)) && ...
            (stateL(idxPiLeft) == stateL(idxPi)) && ...
            (stateL(idxPiRight) == stateL(idxPi));

        regularPoint = ...
            profitOptimumIsInterior && sameLocalMarketState;

        conditionValue = CSHsAtPi ...
            + tauV ...
            + omega * v * qPi ...
            - CSLsAtPi;

        Prop9ConditionValue(iRatio,iV) = conditionValue;
        Prop9RegularPoint(iRatio,iV) = regularPoint;
        Prop9Condition(iRatio,iV) = ...
            regularPoint && ...
            (conditionValue >= -tolMarginal);

        Prop9Ordering(iRatio,iV) = qSW <= qPi + tolQ;

        Prop9Violation(iRatio,iV) = max(qSW - qPi,0);

    end
end

%% =========================================================
% 6. Derived quantities and diagnostics
%% =========================================================

GapQ = qPiStar - qSWStar;
GapEta = etaPiStar - etaSWStar;

conditionAndOrdering = Prop9Condition & Prop9Ordering;
conditionButViolation = Prop9Condition & (~Prop9Ordering);

% Condition-map codes:
% 0 = baseline SLA infeasible
% 1 = sufficient condition not satisfied
% 2 = condition satisfied and ordering verified
% 3 = condition satisfied but ordering violated
Prop9Map = zeros(nRatio,nV);
Prop9Map(BaselineSLAFeasible) = 1;
Prop9Map(conditionAndOrdering) = 2;
Prop9Map(conditionButViolation) = 3;

nSLAInfeasible = nnz(~BaselineSLAFeasible);
nNo = nnz(RegimeSW == 1);
nPartial = nnz(RegimeSW == 2);
nFull = nnz(RegimeSW == 3);

if nFull > 0

    [fullRatioIndex,fullVIndex] = find(RegimeSW == 3);

    fullVRange = [ ...
        min(vGrid(fullVIndex)), ...
        max(vGrid(fullVIndex))];

    fullRatioRange = [ ...
        min(ratioGrid(fullRatioIndex)), ...
        max(ratioGrid(fullRatioIndex))];

else

    fullVRange = [NaN,NaN];
    fullRatioRange = [NaN,NaN];

end

nCondition = nnz(Prop9Condition);
nRegular = nnz(Prop9RegularPoint);
nConditionalViolations = nnz(conditionButViolation);

if any(Prop9Condition(:))
    maxConditionalViolation = max(Prop9Violation(Prop9Condition));
else
    maxConditionalViolation = NaN;
end

gapQValid = GapQ(BaselineSLAFeasible);
gapEtaValid = GapEta(BaselineSLAFeasible);
boundaryResidualValid = SLAResidualAtQBar(BaselineSLAFeasible);
slackAtOptimaValid = [ ...
    SLASlackAtPi(BaselineSLAFeasible); ...
    SLASlackAtSW(BaselineSLAFeasible)];

if isempty(gapQValid)
    gapQRange = [0,0];
    gapEtaRange = [0,0];
    maxBoundaryResidual = NaN;
    minOptimumSLASlack = NaN;
else
    gapQRange = [min(gapQValid),max(gapQValid)];
    gapEtaRange = [min(gapEtaValid),max(gapEtaValid)];
    maxBoundaryResidual = max(abs(boundaryResidualValid));
    minOptimumSLASlack = min(slackAtOptimaValid);
end

fprintf('\n');
fprintf('============================================================\n');
fprintf('Extended Model 2: Figure 8 numerical diagnostics\n');
fprintf('============================================================\n');
fprintf('Scenario: %s\n',scenarioMode);
fprintf('Common cost scale: %.4f\n',costScale);
fprintf('High-tier SLA: WbarH=%.4f, bH=%.4f\n',WbarH,bH);
fprintf(['Active coefficients: kappa=%.6f, tau0=%.6f, ' ...
    'tau1=%.6f, omega=%.6f\n'], ...
    kappa,tau0,tau1,omega);
fprintf('Grid points: %d\n',nRatio * nV);
fprintf('Baseline SLA-infeasible points:           %d\n',nSLAInfeasible);
fprintf('Welfare optimum -- no backfilling:       %d\n',nNo);
fprintf('Welfare optimum -- partial backfilling:  %d\n',nPartial);
fprintf('Welfare optimum -- full feasible:        %d\n',nFull);

if nFull > 0
    fprintf(['Full-feasible region -- service value:   ' ...
        '[%.4f, %.4f]\n'],fullVRange(1),fullVRange(2));
    fprintf(['Full-feasible region -- SG/SC ratio:     ' ...
        '[%.4f, %.4f]\n'],fullRatioRange(1),fullRatioRange(2));
end
fprintf('\n');
fprintf('Regular interior points for Proposition 9: %d\n',nRegular);
fprintf('Modified Proposition 9 condition holds:  %d points\n',nCondition);
fprintf('Conditional ordering violations:         %d points\n',nConditionalViolations);
fprintf('Maximum conditional violation qSW-qPi:   %.12e\n',maxConditionalViolation);
fprintf('Global range of qPi-qSW:                  [%.8f, %.8f]\n', ...
    gapQRange(1),gapQRange(2));
fprintf('Global range of etaPi-etaSW:              [%.8f, %.8f]\n', ...
    gapEtaRange(1),gapEtaRange(2));
fprintf('Maximum absolute qBar boundary residual:  %.12e\n', ...
    maxBoundaryResidual);
fprintf('Minimum SLA slack at reported optima:     %.12e\n', ...
    minOptimumSLASlack);

if nFull == 0
    fprintf('\n');
    fprintf(['NOTE: No full-feasible welfare optimum appears under the\n' ...
        'active scenario calibration. This is a numerical result,\n' ...
        'not a program error.\n']);
end

if nConditionalViolations == 0
    fprintf('\nModified Proposition 9 is numerically consistent on all\n');
    fprintf('grid points where its sufficient condition is satisfied.\n');
else
    fprintf('\nWARNING: Proposition 9 ordering is violated at some points\n');
    fprintf('where the coded sufficient condition is satisfied. Check\n');
    fprintf('concavity, numerical tolerance, and manuscript formulas.\n');
end

fprintf('============================================================\n\n');

%% 7. Save numerical results

resultsFile = sprintf( ...
    'Figure8_endogenous_results_%s.mat', ...
    scenarioFileTag);

save(resultsFile, ...
    'scenarioMode','scenarioLabel','scenarioFileTag','costScale', ...
    'WbarH','bH','tolSLA','nBisection', ...
    'kappaBenchmark','tau0Benchmark','tau1Benchmark','omegaBenchmark', ...
    'kappa','tau0','tau1','omega', ...
    'vGrid','ratioGrid', ...
    'BaselineSLAFeasible','SLAResidualAtQBar', ...
    'SLASlackAtPi','SLASlackAtSW', ...
    'qBarAll','qPiStar','qSWStar', ...
    'etaBarFeasible','etaPiStar','etaSWStar', ...
    'PiStar','WStar','RegimeSW', ...
    'GapQ','GapEta', ...
    'Prop9ConditionValue','Prop9RegularPoint','Prop9Condition', ...
    'Prop9Ordering','Prop9Violation','Prop9Map', ...
    'fullVRange','fullRatioRange', ...
    'nSLAInfeasible','maxBoundaryResidual','minOptimumSLASlack');
%% =========================================================
% 8. Publication-oriented three-panel Figure 8
%
% Panel (a): welfare-optimal regime under Proposition 8
% Panel (b): welfare-optimal actual backfilling amount qSWStar
% Panel (c): actual profit-welfare gap qPiStar - qSWStar
%
% The modified Proposition 9 condition remains checked in the
% command-window diagnostics; the uninformative all-green
% condition map is intentionally not plotted.
%% =========================================================

fontName = 'Times New Roman';

axisFontSize = 10.5;
labelFontSize = 12;
titleFontSize = 11.5;
colorbarFontSize = 9.5;

% Panel (a) reports the global phase structure on v in [0,6].
% Panels (b) and (c) use the same enlarged detail range. It covers
% market entry, endogenous backfilling, the platform-welfare gap,
% and the return to a corner regime without wasting space on the
% long zero-valued tail of the benchmark calibration.

useDetailRange = true;
detailVLim = [1.20,3.20];

% Sequential color maps
nColor = 256;

blueMap = [ ...
    linspace(1.00,0.12,nColor)', ...
    linspace(1.00,0.38,nColor)', ...
    linspace(1.00,0.72,nColor)'];

% Zero-centered diverging map for qPiStar-qSWStar. Negative values
% are blue, zero is white, and positive values are red.
nHalf = nColor / 2;

blueToWhite = [ ...
    linspace(0.15,1.00,nHalf)', ...
    linspace(0.38,1.00,nHalf)', ...
    linspace(0.72,1.00,nHalf)'];

whiteToRed = [ ...
    linspace(1.00,0.72,nHalf)', ...
    linspace(1.00,0.12,nHalf)', ...
    linspace(1.00,0.12,nHalf)'];

gapMap = [blueToWhite;whiteToRed];

redSequentialMap = [ ...
    linspace(1.00,0.72,nColor)', ...
    linspace(1.00,0.12,nColor)', ...
    linspace(1.00,0.12,nColor)'];

negativeSequentialMap = [ ...
    linspace(0.15,1.00,nColor)', ...
    linspace(0.38,1.00,nColor)', ...
    linspace(0.72,1.00,nColor)'];

regimeColors = [ ...
    0.50 0.50 0.50; ...   % baseline SLA infeasible
    0.86 0.86 0.86; ...   % no backfilling
    0.94 0.64 0.20; ...   % partial backfilling
    0.30 0.62 0.82];      % full feasible backfilling

%% 9. Create the three-panel figure

fig8 = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[1.5,2.0,31.0,10.5], ...
    'Renderer','painters');

layout8 = tiledlayout( ...
    fig8, ...
    1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');



%% ---------------------------------------------------------
% Panel (a): Welfare-optimal regime
%% ---------------------------------------------------------

ax1 = nexttile(layout8,1);

imagesc(ax1,vGrid,ratioGrid,RegimeSW);
set(ax1,'YDir','normal');
hold(ax1,'on');

colormap(ax1,regimeColors);
caxis(ax1,[-0.5,3.5]);

cb1 = colorbar(ax1,'eastoutside');
cb1.Ticks = [0,1,2,3];
cb1.TickLabels = {'SLA infeasible','No','Partial','Full feasible'};
cb1.FontName = fontName;
cb1.FontSize = colorbarFontSize;
cb1.TickDirection = 'out';

% Explicit boundaries improve readability without changing the
% underlying classification. Dotted: SLA-feasibility boundary;
% solid: no/partial boundary; dashed: partial/full boundary.
contour( ...
    ax1, ...
    vGrid,ratioGrid,RegimeSW, ...
    [0.5,0.5], ...
    'k:', ...
    'LineWidth',1.0);

contour( ...
    ax1, ...
    vGrid,ratioGrid,RegimeSW, ...
    [1.5,1.5], ...
    'k-', ...
    'LineWidth',0.9);

if nFull > 0

    contour( ...
        ax1, ...
        vGrid,ratioGrid,RegimeSW, ...
        [2.5,2.5], ...
        'k--', ...
        'LineWidth',1.1);

end

% Dotted vertical lines identify the enlarged service-value interval
% displayed in panels (b) and (c).
if useDetailRange

    xline( ...
        ax1,detailVLim(1), ...
        ':', ...
        'Color',[0.25,0.25,0.25], ...
        'LineWidth',0.8);

    xline( ...
        ax1,detailVLim(2), ...
        ':', ...
        'Color',[0.25,0.25,0.25], ...
        'LineWidth',0.8);

end

title(ax1, ...
    '(a) Welfare-optimal regime', ...
    'FontName',fontName, ...
    'FontSize',titleFontSize, ...
    'FontWeight','bold');

xlabel(ax1, ...
    'Service value  $v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',labelFontSize);

ylabel(ax1, ...
    '$S_G/S_C$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',labelFontSize);

xlim(ax1,[vMin,vMax]);
ylim(ax1,[ratioMin,ratioMax]);

if nFull == 0

    text(ax1, ...
        0.97,0.045, ...
        'No full-feasible region', ...
        'Units','normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom', ...
        'FontName',fontName, ...
        'FontSize',8.5, ...
        'FontAngle','italic', ...
        'Color',[0.25,0.25,0.25]);

end

%% ---------------------------------------------------------
% Panel (b): Welfare-optimal actual backfilling amount
%% ---------------------------------------------------------

ax2 = nexttile(layout8,2);

img2 = imagesc(ax2,vGrid,ratioGrid,qSWStar);
set(img2,'AlphaData',double(BaselineSLAFeasible));
set(ax2,'YDir','normal');
set(ax2,'Color',regimeColors(1,:));
hold(ax2,'on');

colormap(ax2,blueMap);

qSWValid = qSWStar(BaselineSLAFeasible);

if isempty(qSWValid)
    qSWColorMax = 1;
else
    qSWColorMax = max(qSWValid);
end

if qSWColorMax <= tolValue
    qSWColorMax = 1;
end

caxis(ax2,[0,qSWColorMax]);

contour(ax2,vGrid,ratioGrid,double(BaselineSLAFeasible), ...
    [0.5,0.5],'k:','LineWidth',1.0);

cb2 = colorbar(ax2,'eastoutside');
cb2.FontName = fontName;
cb2.FontSize = colorbarFontSize;
cb2.TickDirection = 'out';
cb2.Label.String = '$q_{SW}^*$';
cb2.Label.Interpreter = 'latex';
cb2.Label.FontName = fontName;
cb2.Label.FontSize = labelFontSize;

title(ax2, ...
    '(b) Welfare-optimal amount', ...
    'FontName',fontName, ...
    'FontSize',titleFontSize, ...
    'FontWeight','bold');

xlabel(ax2, ...
    'Service value  $v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',labelFontSize);

ylim(ax2,[ratioMin,ratioMax]);

if useDetailRange
    xlim(ax2,detailVLim);
else
    xlim(ax2,[vMin,vMax]);
end

%% ---------------------------------------------------------
% Panel (c): Profit-welfare gap in actual backfilling amount
%% ---------------------------------------------------------

ax3 = nexttile(layout8,3);

img3 = imagesc(ax3,vGrid,ratioGrid,GapQ);
set(img3,'AlphaData',double(BaselineSLAFeasible));
set(ax3,'YDir','normal');
set(ax3,'Color',regimeColors(1,:));
hold(ax3,'on');

gapQMin = gapQRange(1);
gapQMax = gapQRange(2);

% Use a sequential scale when the ordering has one sign, and a
% zero-centered diverging scale only when both signs are present.
% This avoids compressing the benchmark's informative positive range.
if gapQMin < -tolValue && gapQMax > tolValue

    colormap(ax3,gapMap);

    gapQAbsMax = max(abs([gapQMin,gapQMax]));
    caxis(ax3,[-gapQAbsMax,gapQAbsMax]);

elseif gapQMin >= -tolValue

    colormap(ax3,redSequentialMap);

    if gapQMax <= tolValue
        gapQMax = 1;
    end

    caxis(ax3,[0,gapQMax]);

else

    colormap(ax3,negativeSequentialMap);

    if abs(gapQMin) <= tolValue
        gapQMin = -1;
    end

    caxis(ax3,[gapQMin,0]);

end

contour(ax3,vGrid,ratioGrid,double(BaselineSLAFeasible), ...
    [0.5,0.5],'k:','LineWidth',1.0);

cb3 = colorbar(ax3,'eastoutside');
cb3.FontName = fontName;
cb3.FontSize = colorbarFontSize;
cb3.TickDirection = 'out';
cb3.Label.String = '$q_\Pi^*-q_{SW}^*$';
cb3.Label.Interpreter = 'latex';
cb3.Label.FontName = fontName;
cb3.Label.FontSize = labelFontSize;

title(ax3, ...
    '(c) Platform--welfare gap', ...
    'FontName',fontName, ...
    'FontSize',titleFontSize, ...
    'FontWeight','bold');

xlabel(ax3, ...
    'Service value  $v$', ...
    'Interpreter','latex', ...
    'FontName',fontName, ...
    'FontSize',labelFontSize);

ylim(ax3,[ratioMin,ratioMax]);

if useDetailRange
    xlim(ax3,detailVLim);
else
    xlim(ax3,[vMin,vMax]);
end

%% 10. Shared axis formatting

allAxes = [ax1,ax2,ax3];

for iAxis = 1:numel(allAxes)

    set(allAxes(iAxis), ...
        'FontName',fontName, ...
        'FontSize',axisFontSize, ...
        'LineWidth',0.8, ...
        'TickDir','out', ...
        'Box','on', ...
        'Layer','top');

end

% Panels share the same vertical variable. Retaining y tick labels
% only in panel (a) reduces repetition in the final manuscript.
set(ax2,'YTickLabel',[]);
set(ax3,'YTickLabel',[]);

linkaxes(allAxes,'y');

%% 11. Export

drawnow;

figureFileStem = sprintf( ...
    'Figure8_EndogenousBackfilling_ThreePanel_%s', ...
    scenarioFileTag);

exportgraphics( ...
    fig8, ...
    [figureFileStem,'.pdf'], ...
    'ContentType','vector');

exportgraphics( ...
    fig8, ...
    [figureFileStem,'.png'], ...
    'Resolution',600);

exportgraphics( ...
    fig8, ...
    [figureFileStem,'.tif'], ...
    'Resolution',600);

fprintf('\nThree-panel Figure 8 exported successfully.\n');
fprintf('Panel (a): full service-value range [%.2f, %.2f].\n',vMin,vMax);

if useDetailRange
    fprintf('Panels (b)-(c): detailed range [%.2f, %.2f].\n', ...
        detailVLim(1),detailVLim(2));
else
    fprintf('Panels (b)-(c): full service-value range [%.2f, %.2f].\n', ...
        vMin,vMax);
end

fprintf('Files: PDF, PNG, and TIFF.\n\n');
fprintf('Output stem: %s\n',figureFileStem);
fprintf('Numerical results: %s\n\n',resultsFile);

%% =========================================================
% End of the self-contained script
%% =========================================================
