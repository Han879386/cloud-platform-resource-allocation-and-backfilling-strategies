clear; clc; close all;

%% =========================================================
% Revised Figure 9 and Appendix H (H-1 and H-2 only)
% Extended Model 2
%
% Deleted: Appendix Figure H-3 and both of its subplots.
% Retained: Figure 9, Figure H-1, Figure H-2.
%% =========================================================

fprintf('\nRunning revised Figure 9 and Appendix H simulations (H-3 deleted).\n');

%% 1. Benchmark calibration
N = 100;
alpha = 0.30;
NH = alpha * N;
NL = (1 - alpha) * N;

tH = 1;
tL = 1;

pH = 2;
pL = 1.5;

gammaH = 20;
gammaL = 12;

SC = 100;

WbarH = 0.10;
bH = 1 / WbarH;

kappa = 0.10;
tau0 = 0.35;
tau1 = 0.25;
omega = 0.10;

etaStatic = 0.30;

mainRatio = 1.00;
ratioSlices = [0.50,1.00,1.50,2.00];

vMin = 0;
vMax = 6;
vStep = 0.02;
vGrid = vMin:vStep:vMax;

nV = numel(vGrid);
nRatio = numel(ratioSlices);

detailVLim = [1.20,3.20];
zoomVLim = [1.55,2.45];

%% 2. Numerical settings and parameter checks
nCoarse = 401;
nFine = 161;
nRefinePasses = 4;

tolValue = 1.0e-10;
tolCapacity = 1.0e-10;
tolMarginal = 1.0e-8;
tolAcceptance = 1.0e-6;
tolSlope = 1.0e-4;
tolSLA = 1.0e-9;
nBisection = 80;

assert(N > 0 && NH > 0 && NL > 0, 'Market sizes must be positive.');
assert(tH > 0 && tL > 0, 'Matching-cost parameters must be positive.');
assert(pH >= 0 && pL >= 0, 'Prices must be nonnegative.');
assert(gammaH >= 0 && gammaL >= 0, 'Congestion parameters must be nonnegative.');
assert(SC > 0 && all(ratioSlices > 0), 'Capacities and capacity ratios must be positive.');
assert(etaStatic >= 0 && etaStatic < 1, 'etaStatic must lie in [0,1).');
assert(WbarH > 0 && bH > 0, 'The high-tier SLA parameters must be positive.');
assert(any(abs(ratioSlices - mainRatio) <= 1.0e-12), 'mainRatio must be included in ratioSlices.');

%% 3. Storage
qBarAll = nan(nRatio,nV);
qStaticAll = zeros(nRatio,nV);
qIdleAll = nan(nRatio,nV);

BaselineSLAFeasible = false(nRatio,nV);
StaticSLAFeasible = false(nRatio,nV);
StaticSLASlack = nan(nRatio,nV);
SLAResidualAtQBar = nan(nRatio,nV);
SLASlackAtPi = nan(nRatio,nV);
SLASlackAtSW = nan(nRatio,nV);

qPiStar = nan(nRatio,nV);
qSWStar = nan(nRatio,nV);

PiStar = nan(nRatio,nV);
WStar = nan(nRatio,nV);
WAtPiStar = nan(nRatio,nV);

PiStatic = nan(nRatio,nV);
PiIdle = nan(nRatio,nV);
PiNoBackfill = nan(nRatio,nV);

WStatic = nan(nRatio,nV);
WIdle = nan(nRatio,nV);
WNoBackfill = nan(nRatio,nV);

StateHPi = zeros(nRatio,nV);
StateLPi = zeros(nRatio,nV);
StateHSW = zeros(nRatio,nV);
StateLSW = zeros(nRatio,nV);

RegimePi = zeros(nRatio,nV);
RegimeSW = zeros(nRatio,nV);

RegularPi = false(nRatio,nV);
RegularSW = false(nRatio,nV);

Prop9ConditionValue = nan(nRatio,nV);
Prop9Condition = false(nRatio,nV);
Prop9Ordering = false(nRatio,nV);
Prop9Violation = nan(nRatio,nV);

PiCrossV = nan(nRatio,nV);
WCrossV = nan(nRatio,nV);

%% 4. Simulation over capacity ratios and service values
for iRatio = 1:nRatio
    ratio = ratioSlices(iRatio);
    SG = ratio * SC;
    qStatic = etaStatic * SG;

    for iV = 1:nV
        v = vGrid(iV);
        Vh = v - pH;
        Vl = v - pL;
        tauV = tau0 + tau1 * v;

        %% 4.1 Static-rule SLA diagnostic and endogenous qBar(v)
        SH0 = SG;

        if SH0 <= tolCapacity || Vh <= tolValue || gammaH >= SH0 * Vh - tolValue
            DH0 = 0;
        elseif SH0 > NH + tolCapacity && Vh > tH + tolValue && ...
                gammaH <= (SH0 - NH) * (Vh - tH) + tolValue
            DH0 = NH;
        else
            Ah = NH / tH;
            discH0 = (SH0 - Ah * Vh)^2 + 4 * Ah * gammaH;
            DH0 = (SH0 + Ah * Vh - sqrt(max(discH0,0))) / 2;
            DH0 = max(0,min(NH,DH0));
            DH0 = min(DH0,max(SH0 - tolCapacity,0));
        end

        qStaticAll(iRatio,iV) = qStatic;

        SHStatic = SG - qStatic;

        if SHStatic <= tolCapacity || Vh <= tolValue || gammaH >= SHStatic * Vh - tolValue
            DHStatic = 0;
        elseif SHStatic > NH + tolCapacity && Vh > tH + tolValue && ...
                gammaH <= (SHStatic - NH) * (Vh - tH) + tolValue
            DHStatic = NH;
        else
            Ah = NH / tH;
            discHStatic = (SHStatic - Ah * Vh)^2 + 4 * Ah * gammaH;
            DHStatic = (SHStatic + Ah * Vh - sqrt(max(discHStatic,0))) / 2;
            DHStatic = max(0,min(NH,DHStatic));
            DHStatic = min(DHStatic,max(SHStatic - tolCapacity,0));
        end

        StaticSLASlack(iRatio,iV) = SHStatic - DHStatic - bH;
        StaticSLAFeasible(iRatio,iV) = StaticSLASlack(iRatio,iV) >= -tolSLA;

        gAtZero = SG - DH0 - bH;

        if gAtZero < -tolSLA
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

                if SHMid <= tolCapacity || Vh <= tolValue || ...
                        gammaH >= SHMid * Vh - tolValue
                    DHMid = 0;
                elseif SHMid > NH + tolCapacity && Vh > tH + tolValue && ...
                        gammaH <= (SHMid - NH) * (Vh - tH) + tolValue
                    DHMid = NH;
                else
                    Ah = NH / tH;
                    discHMid = (SHMid - Ah * Vh)^2 + 4 * Ah * gammaH;
                    DHMid = (SHMid + Ah * Vh - sqrt(max(discHMid,0))) / 2;
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

            qBar = qLow;
            gAtQBar = gLow;
        end

        qIdle = min(qStatic,qBar);

        qBarAll(iRatio,iV) = qBar;
        qIdleAll(iRatio,iV) = qIdle;
        SLAResidualAtQBar(iRatio,iV) = gAtQBar;

        %% 4.2 Adaptive maximization on [0,qBar]
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

        for refinePass = 0:nRefinePasses
            if qBar <= tolCapacity
                qEval = 0;
            elseif refinePass == 0
                qEval = unique([0,qIdle,linspace(0,qBar,nCoarse),qBar]);
            else
                qEvalPi = linspace(qLeftPi,qRightPi,nFine);
                qEvalSW = linspace(qLeftSW,qRightSW,nFine);
                qEval = unique([0,qIdle,qEvalPi,qEvalSW,qBar]);
            end

            SH = SG - qEval;
            SL = SC + qEval;
            nQEval = numel(qEval);

            DH = zeros(1,nQEval);
            DL = zeros(1,nQEval);
            CSH = zeros(1,nQEval);
            CSL = zeros(1,nQEval);
            CSHs = zeros(1,nQEval);
            CSLs = zeros(1,nQEval);

            RHSv = zeros(1,nQEval);
            RLSv = zeros(1,nQEval);
            FHSv = zeros(1,nQEval);
            FLSv = zeros(1,nQEval);

            %% High-tier market
            emptyH = (SH <= tolCapacity) | (Vh <= tolValue) | ...
                (gammaH >= SH .* Vh - tolValue);

            fullH = (~emptyH) & (SH > NH + tolCapacity) & ...
                (Vh > tH + tolValue) & ...
                (gammaH <= (SH - NH) .* (Vh - tH) + tolValue);

            partialH = (~emptyH) & (~fullH);

            if any(partialH)
                Ah = NH / tH;
                SHp = SH(partialH);

                discH = (SHp - Ah * Vh).^2 + 4 * Ah * gammaH;
                rootH = sqrt(max(discH,0));

                DHp = (SHp + Ah * Vh - rootH) / 2;
                DHp = max(DHp,0);
                DHp = min(DHp,NH);
                DHp = min(DHp,max(SHp - tolCapacity,0));

                DH(partialH) = DHp;
                CSH(partialH) = tH .* DHp.^2 ./ (2 * NH);

                uH = SHp - DHp;
                denH = uH.^2 + Ah * gammaH;

                dDHdSH = Ah * gammaH ./ denH;
                dDHdv = Ah / 2 .* (1 + (SHp - Ah * Vh) ./ rootH);
                dDHdSHdv = 2 * Ah * gammaH .* uH .* dDHdv ./ denH.^2;

                CSHs(partialH) = (tH .* DHp ./ NH) .* dDHdSH;
                CSHsv = (tH / NH) .* (dDHdv .* dDHdSH + DHp .* dDHdSHdv);

                RHSv(partialH) = pH .* dDHdSHdv;
                FHSv(partialH) = RHSv(partialH) + CSHsv;
            end

            if any(fullH)
                SHf = SH(fullH);
                DH(fullH) = NH;
                CSH(fullH) = NH .* (Vh - gammaH ./ (SHf - NH) - tH / 2);
                CSH(fullH) = max(CSH(fullH),0);
                CSHs(fullH) = NH * gammaH ./ (SHf - NH).^2;
            end

            %% Low-tier market
            emptyL = (SL <= tolCapacity) | (Vl <= tolValue) | ...
                (gammaL >= SL .* Vl - tolValue);

            fullL = (~emptyL) & (SL > NL + tolCapacity) & ...
                (Vl > tL / 2 + tolValue) & ...
                (gammaL <= (SL - NL) .* (Vl - tL / 2) + tolValue);

            partialL = (~emptyL) & (~fullL);

            if any(partialL)
                Al = 2 * NL / tL;
                SLp = SL(partialL);

                discL = (SLp - Al * Vl).^2 + 4 * Al * gammaL;
                rootL = sqrt(max(discL,0));

                DLp = (SLp + Al * Vl - rootL) / 2;
                DLp = max(DLp,0);
                DLp = min(DLp,NL);
                DLp = min(DLp,max(SLp - tolCapacity,0));

                DL(partialL) = DLp;
                CSL(partialL) = tL .* DLp.^2 ./ (4 * NL);

                uL = SLp - DLp;
                denL = uL.^2 + Al * gammaL;

                dDLdSL = Al * gammaL ./ denL;
                dDLdv = Al / 2 .* (1 + (SLp - Al * Vl) ./ rootL);
                dDLdSLdv = 2 * Al * gammaL .* uL .* dDLdv ./ denL.^2;

                CSLs(partialL) = (tL .* DLp ./ (2 * NL)) .* dDLdSL;
                CSLsv = (tL / (2 * NL)) .* (dDLdv .* dDLdSL + DLp .* dDLdSLdv);

                RLSv(partialL) = pL .* dDLdSLdv;
                FLSv(partialL) = RLSv(partialL) + CSLsv;
            end

            if any(fullL)
                SLf = SL(fullL);
                DL(fullL) = NL;
                CSL(fullL) = NL .* (Vl - gammaL ./ (SLf - NL) - tL / 4);
                CSL(fullL) = max(CSL(fullL),0);
                CSLs(fullL) = NL * gammaL ./ (SLf - NL).^2;
            end

            %% Profit and welfare objectives
            revenue = pH .* DH + pL .* DL;
            PiEval = revenue - 0.5 * kappa .* qEval.^2;
            WEval = PiEval + CSH + CSL - tauV .* qEval ...
                - 0.5 * omega * v .* qEval.^2;

            [PiBestPass,idxPi] = max(PiEval);
            [WBestPass,idxSW] = max(WEval);

            qPiBestPass = qEval(idxPi);
            qSWBestPass = qEval(idxSW);

            idxPiLeft = max(idxPi - 1,1);
            idxPiRight = min(idxPi + 1,nQEval);
            qLeftPi = qEval(idxPiLeft);
            qRightPi = qEval(idxPiRight);

            idxSWLeft = max(idxSW - 1,1);
            idxSWRight = min(idxSW + 1,nQEval);
            qLeftSW = qEval(idxSWLeft);
            qRightSW = qEval(idxSWRight);
        end

        %% 4.3 Store optima and regimes
        qPi = qPiBestPass;
        qSW = qSWBestPass;

        qPiStar(iRatio,iV) = qPi;
        qSWStar(iRatio,iV) = qSW;
        PiStar(iRatio,iV) = PiBestPass;
        WStar(iRatio,iV) = WBestPass;
        WAtPiStar(iRatio,iV) = WEval(idxPi);

        SLASlackAtPi(iRatio,iV) = SG - qPi - DH(idxPi) - bH;
        SLASlackAtSW(iRatio,iV) = SG - qSW - DH(idxSW) - bH;

        if SLASlackAtPi(iRatio,iV) < -10 * tolSLA || ...
                SLASlackAtSW(iRatio,iV) < -10 * tolSLA
            error('An endogenous optimum violates the high-tier SLA.');
        end

        tolQ = max(1.0e-7,1.0e-6 * max(1,qBar));

        if qPi <= tolQ
            RegimePi(iRatio,iV) = 1;
        elseif qBar - qPi <= tolQ
            RegimePi(iRatio,iV) = 3;
        else
            RegimePi(iRatio,iV) = 2;
        end

        if qSW <= tolQ
            RegimeSW(iRatio,iV) = 1;
        elseif qBar - qSW <= tolQ
            RegimeSW(iRatio,iV) = 3;
        else
            RegimeSW(iRatio,iV) = 2;
        end

        %% 4.4 Regularity
        stateH = zeros(1,nQEval);
        stateH(partialH) = 2;
        stateH(fullH) = 3;

        stateL = zeros(1,nQEval);
        stateL(partialL) = 2;
        stateL(fullL) = 3;

        StateHPi(iRatio,iV) = stateH(idxPi);
        StateLPi(iRatio,iV) = stateL(idxPi);
        StateHSW(iRatio,iV) = stateH(idxSW);
        StateLSW(iRatio,iV) = stateL(idxSW);

        profitOptimumIsInterior = (qPi > tolQ) && (qBar - qPi > tolQ);
        welfareOptimumIsInterior = (qSW > tolQ) && (qBar - qSW > tolQ);

        sameLocalStatePi = ...
            (stateH(idxPiLeft) == stateH(idxPi)) && ...
            (stateH(idxPiRight) == stateH(idxPi)) && ...
            (stateL(idxPiLeft) == stateL(idxPi)) && ...
            (stateL(idxPiRight) == stateL(idxPi));

        sameLocalStateSW = ...
            (stateH(idxSWLeft) == stateH(idxSW)) && ...
            (stateH(idxSWRight) == stateH(idxSW)) && ...
            (stateL(idxSWLeft) == stateL(idxSW)) && ...
            (stateL(idxSWRight) == stateL(idxSW));

        RegularPi(iRatio,iV) = profitOptimumIsInterior && sameLocalStatePi;
        RegularSW(iRatio,iV) = welfareOptimumIsInterior && sameLocalStateSW;

        %% 4.5 Proposition 9
        conditionValue = CSHs(idxPi) + tauV + omega * v * qPi - CSLs(idxPi);
        Prop9ConditionValue(iRatio,iV) = conditionValue;
        Prop9Condition(iRatio,iV) = RegularPi(iRatio,iV) && ...
            (conditionValue >= -tolMarginal);
        Prop9Ordering(iRatio,iV) = qSW <= qPi + tolQ;
        Prop9Violation(iRatio,iV) = max(qSW - qPi,0);

        %% 4.6 Proposition 10 cross partials
        PiCrossV(iRatio,iV) = RLSv(idxPi) - RHSv(idxPi);
        WCrossV(iRatio,iV) = FLSv(idxSW) - FHSv(idxSW) ...
            - tau1 - omega * qSW;

        %% 4.7 Fixed, idle, and no-backfilling rules
        qPolicy = [qStatic,qIdle,0];
        SHPolicy = SG - qPolicy;
        SLPolicy = SC + qPolicy;
        nPolicy = numel(qPolicy);

        DHPolicy = zeros(1,nPolicy);
        DLPolicy = zeros(1,nPolicy);
        CSHPolicy = zeros(1,nPolicy);
        CSLPolicy = zeros(1,nPolicy);

        emptyHPolicy = (SHPolicy <= tolCapacity) | (Vh <= tolValue) | ...
            (gammaH >= SHPolicy .* Vh - tolValue);

        fullHPolicy = (~emptyHPolicy) & (SHPolicy > NH + tolCapacity) & ...
            (Vh > tH + tolValue) & ...
            (gammaH <= (SHPolicy - NH) .* (Vh - tH) + tolValue);

        partialHPolicy = (~emptyHPolicy) & (~fullHPolicy);

        if any(partialHPolicy)
            Ah = NH / tH;
            SHp = SHPolicy(partialHPolicy);
            discH = (SHp - Ah * Vh).^2 + 4 * Ah * gammaH;

            DHp = (SHp + Ah * Vh - sqrt(max(discH,0))) / 2;
            DHp = max(DHp,0);
            DHp = min(DHp,NH);
            DHp = min(DHp,max(SHp - tolCapacity,0));

            DHPolicy(partialHPolicy) = DHp;
            CSHPolicy(partialHPolicy) = tH .* DHp.^2 ./ (2 * NH);
        end

        if any(fullHPolicy)
            SHf = SHPolicy(fullHPolicy);
            DHPolicy(fullHPolicy) = NH;
            CSHPolicy(fullHPolicy) = NH .* (Vh - gammaH ./ (SHf - NH) - tH / 2);
            CSHPolicy(fullHPolicy) = max(CSHPolicy(fullHPolicy),0);
        end

        emptyLPolicy = (SLPolicy <= tolCapacity) | (Vl <= tolValue) | ...
            (gammaL >= SLPolicy .* Vl - tolValue);

        fullLPolicy = (~emptyLPolicy) & (SLPolicy > NL + tolCapacity) & ...
            (Vl > tL / 2 + tolValue) & ...
            (gammaL <= (SLPolicy - NL) .* (Vl - tL / 2) + tolValue);

        partialLPolicy = (~emptyLPolicy) & (~fullLPolicy);

        if any(partialLPolicy)
            Al = 2 * NL / tL;
            SLp = SLPolicy(partialLPolicy);
            discL = (SLp - Al * Vl).^2 + 4 * Al * gammaL;

            DLp = (SLp + Al * Vl - sqrt(max(discL,0))) / 2;
            DLp = max(DLp,0);
            DLp = min(DLp,NL);
            DLp = min(DLp,max(SLp - tolCapacity,0));

            DLPolicy(partialLPolicy) = DLp;
            CSLPolicy(partialLPolicy) = tL .* DLp.^2 ./ (4 * NL);
        end

        if any(fullLPolicy)
            SLf = SLPolicy(fullLPolicy);
            DLPolicy(fullLPolicy) = NL;
            CSLPolicy(fullLPolicy) = NL .* (Vl - gammaL ./ (SLf - NL) - tL / 4);
            CSLPolicy(fullLPolicy) = max(CSLPolicy(fullLPolicy),0);
        end

        revenuePolicy = pH .* DHPolicy + pL .* DLPolicy;
        PiPolicy = revenuePolicy - 0.5 * kappa .* qPolicy.^2;
        WPolicy = PiPolicy + CSHPolicy + CSLPolicy ...
            - tauV .* qPolicy - 0.5 * omega * v .* qPolicy.^2;

        PiStatic(iRatio,iV) = PiPolicy(1);
        PiIdle(iRatio,iV) = PiPolicy(2);
        PiNoBackfill(iRatio,iV) = PiPolicy(3);

        WStatic(iRatio,iV) = WPolicy(1);
        WIdle(iRatio,iV) = WPolicy(2);
        WNoBackfill(iRatio,iV) = WPolicy(3);
    end
end

%% 5. Derived series, Proposition 10 derivatives, and checks
GapQ = qPiStar - qSWStar;

DeltaW_ProfitStatic = WAtPiStar - WStatic;
DeltaW_WelfareStatic = WStar - WStatic;
DeltaW_EndoStatic = DeltaW_WelfareStatic;
DeltaW_EndoIdle = WStar - WIdle;
DeltaW_EndoNoBackfill = WStar - WNoBackfill;

dqPiDv = nan(nRatio,nV);
dqSWDv = nan(nRatio,nV);

for iRatio = 1:nRatio
    for iV = 2:(nV - 1)
        localPiRegular = all(RegularPi(iRatio,(iV - 1):(iV + 1)));
        samePiStatesAcrossV = ...
            all(StateHPi(iRatio,(iV - 1):(iV + 1)) == StateHPi(iRatio,iV)) && ...
            all(StateLPi(iRatio,(iV - 1):(iV + 1)) == StateLPi(iRatio,iV));

        if localPiRegular && samePiStatesAcrossV
            dqPiDv(iRatio,iV) = ...
                (qPiStar(iRatio,iV + 1) - qPiStar(iRatio,iV - 1)) ...
                / (vGrid(iV + 1) - vGrid(iV - 1));
        end

        localSWRegular = all(RegularSW(iRatio,(iV - 1):(iV + 1)));
        sameSWStatesAcrossV = ...
            all(StateHSW(iRatio,(iV - 1):(iV + 1)) == StateHSW(iRatio,iV)) && ...
            all(StateLSW(iRatio,(iV - 1):(iV + 1)) == StateLSW(iRatio,iV));

        if localSWRegular && sameSWStatesAcrossV
            dqSWDv(iRatio,iV) = ...
                (qSWStar(iRatio,iV + 1) - qSWStar(iRatio,iV - 1)) ...
                / (vGrid(iV + 1) - vGrid(iV - 1));
        end
    end
end

Prop10PiCondition = isfinite(dqPiDv) & (PiCrossV < -tolMarginal);
Prop10SWCondition = isfinite(dqSWDv) & (WCrossV < -tolMarginal);
Prop10PiViolation = Prop10PiCondition & (dqPiDv > tolSlope);
Prop10SWViolation = Prop10SWCondition & (dqSWDv > tolSlope);
Prop9ConditionalViolation = Prop9Condition & (~Prop9Ordering);

validOptima = BaselineSLAFeasible & isfinite(qBarAll) & ...
    isfinite(qPiStar) & isfinite(qSWStar);

piLowerViolation = max(-qPiStar(validOptima),0);
swLowerViolation = max(-qSWStar(validOptima),0);
piUpperViolation = max(qPiStar(validOptima) - qBarAll(validOptima),0);
swUpperViolation = max(qSWStar(validOptima) - qBarAll(validOptima),0);
idleUpperViolation = max(qIdleAll(validOptima) - qBarAll(validOptima),0);

welfareAdvantageOverIdle = WStar - WIdle;
welfarePlannerAdvantage = WStar - WAtPiStar;

maxPiLowerViolation = max(piLowerViolation);
maxSWLowerViolation = max(swLowerViolation);
maxPiUpperViolation = max(piUpperViolation);
maxSWUpperViolation = max(swUpperViolation);
maxIdleUpperViolation = max(idleUpperViolation);

minWelfareAdvantageOverIdle = min(welfareAdvantageOverIdle(validOptima));
minWelfarePlannerAdvantage = min(welfarePlannerAdvantage(validOptima));

boundaryResidualValid = SLAResidualAtQBar(validOptima);
slackAtOptimaValid = [SLASlackAtPi(validOptima); SLASlackAtSW(validOptima)];

maxBoundaryResidual = max(abs(boundaryResidualValid));
minOptimumSLASlack = min(slackAtOptimaValid);
nBaselineSLAInfeasible = nnz(~BaselineSLAFeasible);
nStaticSLAInfeasible = nnz(~StaticSLAFeasible);

assert(maxPiLowerViolation <= tolAcceptance, 'A profit optimum violates q>=0.');
assert(maxSWLowerViolation <= tolAcceptance, 'A welfare optimum violates q>=0.');
assert(maxPiUpperViolation <= tolAcceptance, 'A profit optimum violates q<=qBar.');
assert(maxSWUpperViolation <= tolAcceptance, 'A welfare optimum violates q<=qBar.');
assert(maxIdleUpperViolation <= tolAcceptance, 'The idle-constrained benchmark violates q<=qBar.');
assert(minWelfareAdvantageOverIdle >= -tolAcceptance, 'W(qSWStar) is below W(qIdle) beyond tolerance.');
assert(minWelfarePlannerAdvantage >= -tolAcceptance, 'W(qSWStar) is below W(qPiStar) beyond tolerance.');
assert(maxBoundaryResidual <= 10 * tolSLA, 'The computed qBar boundary residual is too large.');
assert(minOptimumSLASlack >= -10 * tolSLA, 'An endogenous optimum violates the high-tier SLA.');

[~,mainIndex] = min(abs(ratioSlices - mainRatio));
mainRatioUsed = ratioSlices(mainIndex);

mainStaticIdleDifferenceSeries = abs(WStatic(mainIndex,:) - WIdle(mainIndex,:));
mainStaticIdleDifferenceSeries = mainStaticIdleDifferenceSeries(isfinite(mainStaticIdleDifferenceSeries));
mainStaticIdleDifference = max(mainStaticIdleDifferenceSeries);
staticIdleCoincide = mainStaticIdleDifference <= tolAcceptance;

%% 6. Command-window diagnostics
fprintf('\n============================================================\n');
fprintf('Extended Model 2: Figure 9 and Appendix H-1/H-2 diagnostics\n');
fprintf('============================================================\n');
fprintf('Benchmark coefficients: kappa=%.2f, tau(v)=%.2f+%.2fv, omega=%.2f\n', ...
    kappa,tau0,tau1,omega);
fprintf('High-tier SLA: WbarH=%.4f, bH=%.4f\n',WbarH,bH);
fprintf('Static benchmark ratio: etaStatic=%.2f\n',etaStatic);
fprintf('Main-text capacity ratio: SG/SC=%.2f\n',mainRatioUsed);
fprintf('Service-value grid: [%.2f,%.2f], step %.3f\n',vMin,vMax,vStep);
fprintf('Baseline SLA-infeasible points: %d\n',nBaselineSLAInfeasible);
fprintf('Static-rule SLA-infeasible points: %d\n',nStaticSLAInfeasible);

fprintf('\nAcceptance checks\n');
fprintf('Maximum lower-bound violation, qPiStar: %.3e\n',maxPiLowerViolation);
fprintf('Maximum lower-bound violation, qSWStar:  %.3e\n',maxSWLowerViolation);
fprintf('Maximum upper-bound violation, qPiStar: %.3e\n',maxPiUpperViolation);
fprintf('Maximum upper-bound violation, qSWStar:  %.3e\n',maxSWUpperViolation);
fprintf('Maximum qIdle-qBar violation:            %.3e\n',maxIdleUpperViolation);
fprintf('Minimum W(qSWStar)-W(qIdle):              %.8f\n',minWelfareAdvantageOverIdle);
fprintf('Maximum absolute qBar boundary residual: %.3e\n',maxBoundaryResidual);
fprintf('Minimum SLA slack at endogenous optima:  %.3e\n',minOptimumSLASlack);

fprintf('\nMain-text Figure 9\n');
mainProfitGain = DeltaW_ProfitStatic(mainIndex,:);
mainWelfareGain = DeltaW_WelfareStatic(mainIndex,:);
mainProfitGain = mainProfitGain(isfinite(mainProfitGain));
mainWelfareGain = mainWelfareGain(isfinite(mainWelfareGain));
fprintf('Platform-optimal gain vs static: [%.8f, %.8f]\n',min(mainProfitGain),max(mainProfitGain));
fprintf('Welfare-optimal gain vs static:  [%.8f, %.8f]\n',min(mainWelfareGain),max(mainWelfareGain));
fprintf('Minimum W(qSWStar)-W(qPiStar):     %.8f\n',minWelfarePlannerAdvantage);

fprintf('\nProposition 8 regime counts by SG/SC\n');
for iRatio = 1:nRatio
    fprintf(['SG/SC=%.2f: baseline SLA-infeasible=%d, ' ...
        'static SLA-infeasible=%d; profit [No=%d, Partial=%d, Full=%d]; ' ...
        'welfare [No=%d, Partial=%d, Full=%d]\n'], ...
        ratioSlices(iRatio), ...
        nnz(RegimePi(iRatio,:) == 0), ...
        nnz(~StaticSLAFeasible(iRatio,:)), ...
        nnz(RegimePi(iRatio,:) == 1), ...
        nnz(RegimePi(iRatio,:) == 2), ...
        nnz(RegimePi(iRatio,:) == 3), ...
        nnz(RegimeSW(iRatio,:) == 1), ...
        nnz(RegimeSW(iRatio,:) == 2), ...
        nnz(RegimeSW(iRatio,:) == 3));
end

fprintf('\nProposition 9 checks\n');
for iRatio = 1:nRatio
    gapValid = GapQ(iRatio,isfinite(GapQ(iRatio,:)));
    fprintf(['SG/SC=%.2f: regular=%d, condition holds=%d, ' ...
        'conditional violations=%d, gap range=[%.8f,%.8f]\n'], ...
        ratioSlices(iRatio), ...
        nnz(RegularPi(iRatio,:)), ...
        nnz(Prop9Condition(iRatio,:)), ...
        nnz(Prop9ConditionalViolation(iRatio,:)), ...
        min(gapValid),max(gapValid));
end

fprintf('\nProposition 10 checks\n');
for iRatio = 1:nRatio
    fprintf(['SG/SC=%.2f: finite d(qPi*)/dv=%d, condition=%d, ' ...
        'violations=%d; finite d(qSW*)/dv=%d, condition=%d, violations=%d\n'], ...
        ratioSlices(iRatio), ...
        nnz(isfinite(dqPiDv(iRatio,:))), ...
        nnz(Prop10PiCondition(iRatio,:)), ...
        nnz(Prop10PiViolation(iRatio,:)), ...
        nnz(isfinite(dqSWDv(iRatio,:))), ...
        nnz(Prop10SWCondition(iRatio,:)), ...
        nnz(Prop10SWViolation(iRatio,:)));
end

if nnz(Prop9ConditionalViolation) == 0
    fprintf('\nProposition 9 has no conditional ordering violation.\n');
else
    warning('Proposition 9 has conditional ordering violations.');
end

if nnz(Prop10PiViolation) == 0 && nnz(Prop10SWViolation) == 0
    fprintf('Proposition 10 has no conditional slope-sign violation.\n');
else
    warning('Proposition 10 has conditional slope-sign violations.');
end
fprintf('============================================================\n\n');

%% 7. Save numerical results
resultsFile = 'Figure9_AppendixH12_EndogenousBackfilling_Benchmark.mat';

save(resultsFile, ...
    'N','alpha','NH','NL','tH','tL','pH','pL', ...
    'gammaH','gammaL','SC','WbarH','bH','tolSLA','nBisection', ...
    'kappa','tau0','tau1','omega','etaStatic', ...
    'mainRatio','mainRatioUsed','ratioSlices', ...
    'vMin','vMax','vStep','vGrid','detailVLim','zoomVLim', ...
    'BaselineSLAFeasible','StaticSLAFeasible','StaticSLASlack', ...
    'SLAResidualAtQBar','SLASlackAtPi','SLASlackAtSW', ...
    'qBarAll','qStaticAll','qIdleAll','qPiStar','qSWStar', ...
    'PiStar','WStar','WAtPiStar', ...
    'PiStatic','PiIdle','PiNoBackfill', ...
    'WStatic','WIdle','WNoBackfill', ...
    'DeltaW_ProfitStatic','DeltaW_WelfareStatic', ...
    'DeltaW_EndoStatic','DeltaW_EndoIdle', ...
    'DeltaW_EndoNoBackfill','GapQ', ...
    'StateHPi','StateLPi','StateHSW','StateLSW', ...
    'RegimePi','RegimeSW','RegularPi','RegularSW', ...
    'Prop9ConditionValue','Prop9Condition','Prop9Ordering', ...
    'Prop9Violation','Prop9ConditionalViolation', ...
    'PiCrossV','WCrossV','dqPiDv','dqSWDv', ...
    'Prop10PiCondition','Prop10SWCondition', ...
    'Prop10PiViolation','Prop10SWViolation', ...
    'staticIdleCoincide','mainStaticIdleDifference', ...
    'welfarePlannerAdvantage','minWelfarePlannerAdvantage', ...
    'nBaselineSLAInfeasible','nStaticSLAInfeasible', ...
    'maxBoundaryResidual','minOptimumSLASlack');

%% 8. Publication formatting
fontName = 'Times New Roman';
axisFontSize = 11;
labelFontSize = 13;
titleFontSize = 12;
legendFontSize = 10;

blue = [0.00,0.45,0.74];
orange = [0.85,0.33,0.10];
green = [0.47,0.67,0.19];
purple = [0.49,0.18,0.56];
gray = [0.30,0.30,0.30];

ratioColors = [blue; orange; green; purple];
ratioLineStyles = {'-','--','-.',':'};

%% 9. Main-text Figure 9: welfare gains
fig9 = figure('Color','w','Units','centimeters', ...
    'Position',[2.0,2.0,17.0,12.0],'Renderer','painters');

ax9 = axes(fig9,'Position',[0.12,0.13,0.82,0.79]);
hold(ax9,'on');

hPlatform9 = plot(ax9,vGrid,DeltaW_ProfitStatic(mainIndex,:), ...
    '-','Color',blue,'LineWidth',2.2);

hWelfare9 = plot(ax9,vGrid,DeltaW_WelfareStatic(mainIndex,:), ...
    '--','Color',orange,'LineWidth',2.0, ...
    'Marker','o','MarkerIndices',1:20:nV, ...
    'MarkerSize',4.2,'MarkerFaceColor','w');

yline(ax9,0,':','Color',[0.45,0.45,0.45],'LineWidth',1.0, ...
    'HandleVisibility','off');

xlabel(ax9,'Service value  $v$','Interpreter','latex', ...
    'FontName',fontName,'FontSize',labelFontSize);

ylabel(ax9,'Welfare improvement','FontName',fontName, ...
    'FontSize',labelFontSize);

etaLegend9 = sprintf('%.2f',etaStatic);
legendLabels9 = { ...
    ['$W(q_\Pi^*)-W(q^{\mathrm{Static}}),\ ' ...
        '\eta^{\mathrm{Static}}=',etaLegend9,'$'], ...
    ['$W(q_{SW}^*)-W(q^{\mathrm{Static}}),\ ' ...
        '\eta^{\mathrm{Static}}=',etaLegend9,'$']};

legend(ax9,[hPlatform9,hWelfare9],legendLabels9, ...
    'Interpreter','latex','FontName',fontName, ...
    'FontSize',legendFontSize,'Location','northwest','Box','on');

xlim(ax9,[vMin,vMax]);

yData9 = [DeltaW_ProfitStatic(mainIndex,:), ...
    DeltaW_WelfareStatic(mainIndex,:)];
yData9 = yData9(isfinite(yData9));

if isempty(yData9)
    yMin9 = 0;
    yMax9 = 1;
else
    yMin9 = min(0,min(yData9));
    yMax9 = max(yData9);
end

if yMax9 <= yMin9 + tolValue
    yMax9 = yMin9 + 1;
end

ylim(ax9,[yMin9,yMin9 + 1.10 * (yMax9 - yMin9)]);

set(ax9,'FontName',fontName,'FontSize',axisFontSize, ...
    'LineWidth',0.8,'TickDir','out','Box','on','Layer','top');
grid(ax9,'on');
drawnow;

exportgraphics(fig9,'Figure9_WelfareGains_Endogenous_Benchmark.pdf', ...
    'ContentType','vector');
exportgraphics(fig9,'Figure9_WelfareGains_Endogenous_Benchmark.png', ...
    'Resolution',600);
exportgraphics(fig9,'Figure9_WelfareGains_Endogenous_Benchmark.tif', ...
    'Resolution',600);

%% 10. Appendix Figure H-1: optimal backfilling amounts
figH1 = figure('Color','w','Units','centimeters', ...
    'Position',[2.0,2.0,17.0,12.5],'Renderer','painters');

mainAxH1 = axes(figH1,'Position',[0.12,0.13,0.82,0.79]);
hold(mainAxH1,'on');

plot(mainAxH1,vGrid,qBarAll(mainIndex,:),'k--','LineWidth',2.1);
plot(mainAxH1,vGrid,qPiStar(mainIndex,:),'-','Color',blue,'LineWidth',2.2);
plot(mainAxH1,vGrid,qSWStar(mainIndex,:),'-','Color',orange,'LineWidth',2.2);

xlabel(mainAxH1,'Service value  $v$','Interpreter','latex', ...
    'FontName',fontName,'FontSize',labelFontSize);

ylabel(mainAxH1,'Actual backfilling amount  $q$','Interpreter','latex', ...
    'FontName',fontName,'FontSize',labelFontSize);

title(mainAxH1,sprintf(['Optimal endogenous backfilling, ' ...
    '$S_G/S_C=%.2f$, $b_h=%.2f$'],mainRatioUsed,bH), ...
    'Interpreter','latex','FontName',fontName, ...
    'FontSize',titleFontSize,'FontWeight','normal');

legend(mainAxH1,{'$\bar q(v)$','$q_\Pi^*(v)$','$q_{SW}^*(v)$'}, ...
    'Interpreter','latex','FontName',fontName, ...
    'FontSize',legendFontSize,'Location','northeast','Box','on');

xlim(mainAxH1,[vMin,vMax]);

qBarMainValues = qBarAll(mainIndex,:);
qBarMainValues = qBarMainValues(isfinite(qBarMainValues));

if isempty(qBarMainValues)
    qBarMainMax = 1;
else
    qBarMainMax = max(qBarMainValues);
end

if qBarMainMax <= tolValue
    qBarMainMax = 1;
end

ylim(mainAxH1,[0,1.08 * qBarMainMax]);

zoomIndex = vGrid >= zoomVLim(1) & vGrid <= zoomVLim(2);

zoomYMax = 1.15 * max([qPiStar(mainIndex,zoomIndex), ...
    qSWStar(mainIndex,zoomIndex)]);

if ~isfinite(zoomYMax)
    zoomYMax = 1;
end

if zoomYMax <= tolValue
    zoomYMax = 1;
end

rectangle(mainAxH1,'Position',[zoomVLim(1),0,diff(zoomVLim),zoomYMax], ...
    'EdgeColor',gray,'LineStyle','--','LineWidth',1.1);

insetAxH1 = axes(figH1,'Position',[0.54,0.43,0.32,0.30]);
hold(insetAxH1,'on');

plot(insetAxH1,vGrid,qPiStar(mainIndex,:),'-','Color',blue,'LineWidth',1.9);
plot(insetAxH1,vGrid,qSWStar(mainIndex,:),'-','Color',orange,'LineWidth',1.9);

xlim(insetAxH1,zoomVLim);
ylim(insetAxH1,[0,zoomYMax]);

title(insetAxH1,'Interior detail','FontName',fontName, ...
    'FontSize',10,'FontWeight','normal');

legend(insetAxH1,{'$q_\Pi^*$','$q_{SW}^*$'}, ...
    'Interpreter','latex','FontName',fontName, ...
    'FontSize',8.5,'Location','northeast','Box','off');

set([mainAxH1,insetAxH1],'FontName',fontName, ...
    'LineWidth',0.8,'TickDir','out','Box','on','Layer','top');
set(mainAxH1,'FontSize',axisFontSize);
set(insetAxH1,'FontSize',9.5);
grid(mainAxH1,'on');
grid(insetAxH1,'on');
drawnow;

exportgraphics(figH1,'FigureH1_OptimalBackfilling_Benchmark.pdf', ...
    'ContentType','vector');
exportgraphics(figH1,'FigureH1_OptimalBackfilling_Benchmark.png', ...
    'Resolution',600);
exportgraphics(figH1,'FigureH1_OptimalBackfilling_Benchmark.tif', ...
    'Resolution',600);

%% 11. Appendix Figure H-2: Proposition 9 gap
figH2 = figure('Color','w','Units','centimeters', ...
    'Position',[2.0,2.0,17.0,12.0],'Renderer','painters');

axH2 = axes(figH2,'Position',[0.12,0.13,0.82,0.79]);
hold(axH2,'on');

hGap = gobjects(1,nRatio);
ratioLabels = cell(1,nRatio);

for iRatio = 1:nRatio
    hGap(iRatio) = plot(axH2,vGrid,GapQ(iRatio,:), ...
        'LineStyle',ratioLineStyles{iRatio}, ...
        'Color',ratioColors(iRatio,:),'LineWidth',2.0);
    ratioLabels{iRatio} = sprintf('$S_G/S_C=%.2f$',ratioSlices(iRatio));
end

yline(axH2,0,':','Color',[0.45,0.45,0.45],'LineWidth',1.0, ...
    'HandleVisibility','off');

xlabel(axH2,'Service value  $v$','Interpreter','latex', ...
    'FontName',fontName,'FontSize',labelFontSize);

ylabel(axH2,'$q_\Pi^*-q_{SW}^*$','Interpreter','latex', ...
    'FontName',fontName,'FontSize',labelFontSize);

title(axH2,'Platform--welfare gap in endogenous backfilling', ...
    'FontName',fontName,'FontSize',titleFontSize,'FontWeight','normal');

legend(axH2,hGap,ratioLabels,'Interpreter','latex', ...
    'FontName',fontName,'FontSize',legendFontSize, ...
    'Location','northeast','Box','on');

xlim(axH2,detailVLim);

gapDetailIndex = vGrid >= detailVLim(1) & vGrid <= detailVLim(2);
gapDetailValues = GapQ(:,gapDetailIndex);
gapDetailValues = gapDetailValues(isfinite(gapDetailValues));

if isempty(gapDetailValues)
    gapDetailMin = 0;
    gapDetailMax = 1;
else
    gapDetailMax = max(gapDetailValues);
    gapDetailMin = min(gapDetailValues);
end

if gapDetailMax <= gapDetailMin + tolValue
    gapDetailMax = gapDetailMin + 1;
end

gapPadding = 0.08 * (gapDetailMax - min(0,gapDetailMin));
ylim(axH2,[min(0,gapDetailMin) - gapPadding, ...
    gapDetailMax + gapPadding]);

text(axH2,0.02,0.96, ...
    ['Proposition 9 is tested only at SLA-feasible, ' ...
    'differentiable interior profit optima'], ...
    'Units','normalized','HorizontalAlignment','left', ...
    'VerticalAlignment','top','FontName',fontName, ...
    'FontSize',9,'FontAngle','italic','Color',gray);

set(axH2,'FontName',fontName,'FontSize',axisFontSize, ...
    'LineWidth',0.8,'TickDir','out','Box','on','Layer','top');
grid(axH2,'on');
drawnow;

exportgraphics(figH2,'FigureH2_ProfitWelfareGap_Benchmark.pdf', ...
    'ContentType','vector');
exportgraphics(figH2,'FigureH2_ProfitWelfareGap_Benchmark.png', ...
    'Resolution',600);
exportgraphics(figH2,'FigureH2_ProfitWelfareGap_Benchmark.tif', ...
    'Resolution',600);

fprintf('Figure 9, Figure H-1, and Figure H-2 were exported successfully.\n');
fprintf('Figure H-3 has been deleted.\n');
fprintf('Numerical output: %s\n\n',resultsFile);

%% =========================================================
% End of self-contained script
%% =========================================================