function [CM_map_result,DI_x_result,DI_y_result] = VDF_main(image_t1,image_t2,opt)
    % pre-processing
    Compactness = 1;
    [Cosup,Ns] = SuperpixelSegmentation(image_t1,image_t2,opt.Ns,Compactness,'VDF');
    CoNs = max(Cosup(:));
    [t1_feature,t2_feature] = MMfeatureExtraction(Cosup,image_t1,image_t2); 
    
    % transfer function
    x = linspace(-1,1,Ns);
    g=@(x) (sign(x-opt.Th)+1)/2;
    [opt.h] = VDF_polynomial_coefficient(g,opt.Norder,x);
    
    iter = 1;
    Kmax =round(size(t1_feature,2).^0.5);
    Kmin = round(Kmax/10);
    labels = zeros(size(t1_feature,2),1);
    [Px_o] = VDF_Graph_Matrix(t1_feature,t1_feature,Kmax,opt);
    [Py_o] = VDF_Graph_Matrix(t2_feature,t2_feature,Kmax,opt);
    [Gx_o] = VDF_Filter_Matrix(Px_o,opt);
    [Gy_o] = VDF_Filter_Matrix(Py_o,opt);
    dist_t1 = pdist2(t1_feature',t1_feature','squaredeuclidean');
    dist_t2 = pdist2(t2_feature',t2_feature','squaredeuclidean');

    idx_Co = label2idx(Cosup);
    while iter <= opt.Niter
        idex_unchange = labels==0;
        idex_changed = labels==1;
        t1_feature_lib = t1_feature;
        t2_feature_lib = t2_feature;
        t1_feature_lib (:,idex_changed) =  inf;
        t2_feature_lib (:,idex_changed) =  inf;
        %--------------------- CalculateChangeLevel----------------------%
        if iter == 1
            Gx = Gx_o;
            Gy = Gy_o;
        else
            [Px] = VDF_Graph_Matrix(t1_feature_lib,t1_feature,Kmax,opt);
            [Py] = VDF_Graph_Matrix(t2_feature_lib,t2_feature,Kmax,opt);
            [Gx] = VDF_Filter_Matrix(Px,opt);
            [Gy] = VDF_Filter_Matrix(Py,opt);
        end
        fx_gsp = sum((Gy-Gx_o).*dist_t1,2);
        fy_gsp = sum((Gx-Gy_o).*dist_t2,2);
        fx_result(:,iter) = remove_outlier(fx_gsp);
        fy_result(:,iter) = remove_outlier(fy_gsp);
        %--------------------- MRF CoSegmentation----------------------%
        if iter == 1
            [CM_map,labels] = InitialOtsu(Cosup,fx_result(:,iter),fy_result(:,iter));
        else
            [CM_map,labels] = MRF_CoSegmentation(Cosup,opt.alpha,t1_feature,t2_feature,fx_result(:,iter),fy_result(:,iter));
        end
        CM_map_result(:,:,iter) = CM_map;
        for i = 1:size(t1_feature,2)
            index_vector = idx_Co{i};
            DI_x(index_vector) = fx_gsp(i);
            DI_y(index_vector) = fy_gsp(i);
        end
        DI_x_result(:,:,iter) = reshape(DI_x,[size(Cosup,1) size(Cosup,2)]);
        DI_y_result(:,:,iter) = reshape(DI_y,[size(Cosup,1) size(Cosup,2)]);
        iter = iter+1;
    end
end

function [P] = VDF_Graph_Matrix(X_library,X,kmax,opt)
    X = X';
    X_library = X_library';
    k = kmax+1;
    [idx, distX] = knnsearch(X_library,X,'k',k);
    [N,~] = size(X);
    A = zeros(N);
    P = zeros(N);
    for i = 1:N
        idx_temp = idx(i,1:k);
        A(i,idx_temp) = 1;
        A(i,i) = 0;
        P(i,idx_temp) = A(i,idx_temp)/sum(A(i,idx_temp));
    end
end

function [Gx] = VDF_Filter_Matrix(Fx,opt)
    Npolynomial = length(opt.h)-1;
    Gx = zeros(size(Fx));
    Gx = sparse(Gx);
    Fx = sparse(Fx);
    for i = 1:Npolynomial
        Gx = opt.h(i)*Fx^(Npolynomial+1-i)+Gx;
    end
    Gx = full(Gx);
end

function [CM_map,labels] = InitialOtsu(Cosup,fx,fy)
    Ic1 = fx/max(fx);
    Ic2 = fy/max(fy);
    T_theory1 = graythresh(Ic1);
    T_theory2 = graythresh(Ic2);
    labels_1 = Ic1>T_theory1;
    labels_2 = Ic2>T_theory2;
    labels = labels_1 | labels_2;
    idx_Co = label2idx(Cosup);
    for i = 1:length(fx)
        index_vector = idx_Co{i};
        CM_map(index_vector) = labels(i);
    end
    CM_map =reshape(CM_map,[size(Cosup,1) size(Cosup,2)]);
end