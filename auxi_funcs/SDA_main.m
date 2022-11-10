function [CM,DI,RegImg] = SDA_main(image_t1,image_t2,opt)
    % pre-processing
    Compactness = 1;
    [sup_img,Ns] =  SuperpixelSegmentation(image_t1,image_t2,opt.Ns,Compactness,'SDA');
    [t1_feature,t2_feature,norm_par] = MMfeatureExtraction(sup_img,image_t1,image_t2) ; %feature extraction

    % Graph Matrix
    Kmax =round(size(t1_feature,2).^0.5);
    [Sx] = SDA_Graph_Matrix (t1_feature,Kmax);

    % transfer function
    Npolynomial = length(opt.h)-1;
    Lx =  LaplacianMatrix(Sx);
    Lx = sparse(Lx);
    Hx = zeros(size(Lx));
    Hx = sparse(Hx);
    for t = 1:Npolynomial
        Hx = opt.h(t)*Lx^(Npolynomial+1-t)+Hx;
    end
    Hx = full(Hx);

    % SDA Regression
    [regression_t1, delt] = SDA_Regression(t2_feature,Hx,opt);% t1--->t2
    DI_tmep = sum(delt.^2,1);
    DI  = suplabel2DI(sup_img,DI_tmep);
    [RegImg,~,~] = suplabel2ImFeature(sup_img,regression_t1,size(image_t2,3));% t1--->t2
    RegImg = DenormImage(RegImg,norm_par(size(image_t1,3)+1:end));

    % MRF segmentation
    [CM,~] = MRFsegmentation(sup_img,opt.alpha,delt);
end


function [S] = SDA_Graph_Matrix (X,kmax)
    X = X';
    kmax = kmax+1;
    kmin = round(kmax/10)+1;
    [idx, distX] = knnsearch(X,X,'k',kmax);
    [N,~] = size(X);
    degree_x = tabulate(idx(:));
    kmat = degree_x(:,2);
    kmat(kmat >= kmax)=kmax;
    kmat(kmat <= kmin)=kmin;
    if length(kmat) < N
        kmat(length(kmat)+1:N) = kmin;
    end
    S = zeros(N,N);
    for i = 1:N
        K = kmat(i);
        k = K-1;
        id_x = idx(i,1:K);
        di = distX(i,1:K);
        W = (di(K)-di)/(k*di(K)-sum(di(1:k))+eps);
        S(i,id_x) = W;
    end
end

function [Lx] = LaplacianMatrix(Sx)
    N = size(Sx,2);
    Lx_temp1 = -(Sx + Sx')/2;
    Lx_temp2 = sum(Lx_temp1,2);
    Lx = Lx_temp1;
    for i = 1:N
        Lx(i,i) = - Lx_temp2(i) + Lx_temp1(i,i);
    end
end

function [DI] = suplabel2DI(sup_cog,suplabel)
    [h,w]   = size(sup_cog);
    B = size(suplabel,1);
    nbr_sp  = max(sup_cog(:));
    idx_t1 = label2idx(sup_cog);
    for b = 1:B
        for i = 1:nbr_sp
            index_vector = idx_t1{i};
            DI_temp(index_vector) = suplabel(b,i);
        end
    DI(:,:,b) =reshape(DI_temp,[h w]);
    end
end

function [Im,feature,sigband] = suplabel2ImFeature(sup_img,X,b)
    feature = suplabel2DI(sup_img,X);
    Im = feature(:,:,1:b);
    sigband = sum(feature.^2,3);
end

function [Img] = DenormImage(X,norm_par)
    b=size(X,3);
    Img = zeros(size(X));
    for i=1:b
        Img(:,:,i) = X(:,:,i)*norm_par(i);
    end
end



