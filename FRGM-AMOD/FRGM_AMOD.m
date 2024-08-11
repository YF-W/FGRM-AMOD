function [FS] = FRGM_AMOD(MV)
    nV = length(MV);
    for i = 1:nV
        MV{i} = MV{i}';
    end
    nS = size(MV{1}, 1);

    G = [0.3, 0.5, 0.7];
    SM = cell(nV, length(G));
    oP = zeros(nS, nV, length(G));
    pP = zeros(nS, nV, length(G));
   
    CES = zeros(nS, nV, nV, length(G));
   
    for gI = 1:length(G)
        g = G(gI);
        for v = 1:nV
            D = MV{v};
            SMtx = FuzzySimilarityMatrix(D, g);
            nSMtx = SMtx;
            SM{v, gI} = nSMtx;

            oP(:, v, gI) = max(unionSets(RB(SMtx, g), SMtx), [], 2);

            pP(:, v, gI) = min(intersectSets(RB(SMtx, g), SMtx), [], 2);
        end
    end

    for gI = 1:length(G)
        for v1 = 1:nV
            SMtx1 = SM{v1, gI};
            for v2 = 1:nV
                if v1 ~= v2
                    SMtx2 = SM{v2, gI};
                    for i = 1:nS
                        Wo1 = oP(i, v1, gI);
                        Wp1 = pP(i, v1, gI);
                        Wo2 = oP(i, v2, gI);
                        Wp2 = pP(i, v2, gI);

                        P = (SMtx1(i, :) - Wp1) / (Wo1 - Wp1 + eps); 
                        Q = (SMtx2(i, :) - Wp2) / (Wo2 - Wp2 + eps); 
                        P = max(P, 0); P = P / sum(P + eps);
                        Q = max(Q, 0); Q = Q / sum(Q + eps);
                        CES(i, v1, v2, gI) = FuzzyCrossEntropy(P, Q);
                    end
                end
            end
        end
    end
    
    FA = (mean(oP, [2, 3]) + mean(pP, [2, 3])) / 2;
    
    nFA = sqrt(sum(FA .^ 2));
    
    th = 1e-10 * nFA;
    
    sA = FA < th;

    W = 1 ./ (1 + exp(-CES)); 
    WCES = CES .* W;
    ACES = sum(WCES, [2, 3, 4]) ./ sum(W, [2, 3, 4]);

    NACES = mat2gray(ACES);

    AOS = NACES;
    AOS(sA) = max(NACES);
    
    FS = (AOS - min(AOS)) / (max(AOS) - min(AOS));
end
    
function SMtx = FuzzySimilarityMatrix(D, G)

    sD = zscore(D);
    D = EuDist2(sD, sD, 1);
    SMtx = exp(-G * D);
    SMtx(1:size(D, 1)+1:end) = 0;
end

function CE = FuzzyCrossEntropy(P, Q)
    P(P == 0) = eps;
    Q(Q == 0) = eps;
    CE = -sum(P .* log(Q) + (1 - P) .* log(1 - Q));
end


function result = RB(SMtx, g)
    result = SMtx * g;
end

function result = unionSets(A, B)
    result = 1 - (1 - A) .* (1 - B);
end

function result = intersectSets(A, B)
    result = A .* B;
end
