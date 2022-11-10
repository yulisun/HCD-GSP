function [Z, delt,RelDiff] = SDA_Regression(Y,Hx,opt)
% X----> Y
    Niter = opt.Niter;
    lambda = opt.lambda;
    mu = 0.4;
    [M,N] = size(Y);
    delt = zeros(M,N);% initialization
    W = zeros(M,N);% initialization
    inv_TmuI = inv(2*Hx + mu * eye(N));% (4*Lx+u*I)^-1
    for i=1:Niter
        delt_old = delt;
        Z = (mu * (Y + delt) - W) * inv_TmuI; % Z update
        Q = Z - Y + W/mu;
        delt = delt_Update(Q,lambda/mu);% delt update£»
        W = W + mu * (Z - Y - delt); % W update
        RelDiff(i) = norm(delt - delt_old,'fro')/norm(delt,'fro');
        if i > 3 && RelDiff(i) < 1e-2
            break
        end
    end
end

function [delt] = delt_Update(Q,alfa)
    Q1 = sqrt(sum(Q.^2,1));
    Q1(Q1==0) = alfa;
    Q2 = (Q1 - alfa) ./ Q1;
    delt = Q * diag((Q1 > alfa) .* Q2);
end

