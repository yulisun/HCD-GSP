
clear;
close all
addpath(genpath(pwd))
%% Parameter setting
% #1-Italy, #2-Texas-L8, #3-Img7, #4-Img17, #5-Shuguang
% #6-California, #7-YellowRiver, #8-Img5, #9-Texas-ALI
% #1-#7 are used in the paper.
dataset = '#1-Italy'; % For other datasets, we recommend a similar pre-processing. 
opt.Ns = 5000;% the number of superpxiels
opt.Th = 0.9; % the cut-off frequency of transfer function
opt.Norder = 4; % the order of transfer function
opt.Niter = 6; % the maximum number of iterations
opt.alpha = 0.05; % the balance parameter of  MRF segmentation
if strcmp(dataset,'#1-Italy') == 1 || strcmp(dataset,'#3-Img7') == 1 || strcmp(dataset,'#6-California') == 1
    opt.alpha = 0.1;
end
VDF_Load_Dataset % Load Dataset

%% VDF-HCD
fprintf(['\n VDF-HCD is running...... ' '\n'])
time = clock;
[CMmap_result,DIx_result,DIy_result] = VDF_main(image_t1,image_t2,opt);
fprintf('\n');fprintf('The total computational time of VDF-HCD is %i \n', etime(clock,time));
fprintf(['\n' '====================================================================== ' '\n'])

%% Displaying results
fprintf(['\n Displaying the results...... ' '\n'])
figure
for i = 1:opt.Niter
    subplot(3,opt.Niter,i);imshow(remove_outlier(DIx_result(:,:,i)),[]);title(sprintf('%dth DIx',i))
    subplot(3,opt.Niter,i+opt.Niter);imshow(remove_outlier(DIy_result(:,:,i)),[]);title(sprintf('%dth DIy',i))
    subplot(3,opt.Niter,i+2*opt.Niter);imshow(CMmap_result(:,:,i));title(sprintf('%dth CM',i))
end
for iter = 1:opt.Niter
    [tp,fp,tn,fn,fplv,fnlv,~,~,pcc(iter),kappa(iter),imw]=performance(CMmap_result(:,:,iter),Ref_gt);
    F1(iter) = 2*tp/(2*tp+fp+fn);
end
[~,iter_idx] = max(pcc+F1+kappa);
CM_map = CMmap_result(:,:,iter_idx);
DIx = DIx_result(:,:,iter_idx);
DIy = DIy_result(:,:,iter_idx);
figure;
subplot(131);imshow(remove_outlier(DIx),[]);title('DIx')
subplot(132);imshow(remove_outlier(DIy),[]);title('DIy')
subplot(133);imshow(CM_map,[]);title('CM')

pcc = pcc(iter_idx);
kappa = kappa(iter_idx);
F1 = F1(iter_idx);
result = 'PCC is %4.3f; KC is %4.3f; F1 is %4.3f\n';
fprintf('\n');fprintf(result,pcc,kappa,F1)

if F1 < 0.3
   fprintf('\n');disp('Please check the opt.type of input images and select the appropriate opt.alpha for VDF-HCD!')
end 


