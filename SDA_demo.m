clear;
close all
addpath(genpath(pwd))
%% Parameter setting
% #1-Italy, #2-Texas-L8, #3-Img7, #4-Img17, #5-Shuguang
% #6-California, #7-YellowRiver, #8-Img5, #9-Texas-ALI
% #1-#7 are used in the paper.
dataset = '#1-Italy'; % For other datasets, we recommend a similar pre-processing. 
opt.Ns = 5000; % the number of superpxiels
opt.h = [1,1,1,0]; % transfer function
opt.Niter = 10; % the maximum iteration number of SDA Regression
opt.lambda = 0.05; % the weight of sparse penalty term 
opt.alpha = 0.01; % the balance parameter of  MRF segmentation
if strcmp(dataset,'#1-Italy') == 1 || strcmp(dataset,'#5-Shuguang') == 1 
    opt.alpha = 0.1;
elseif strcmp(dataset,'#6-California') == 1 
    opt.alpha = 0.5;
end
SDA_Load_Dataset

%% SDA-HCD
fprintf(['\n SDA-HCD is running...... ' '\n'])
time = clock;
[CM,DI,RegImg] = SDA_main(image_t1,image_t2,opt);
fprintf('\n');fprintf('The total computational time of SDA-HCD is %i \n', etime(clock,time));
fprintf(['\n' '====================================================================== ' '\n'])

%% Displaying results
fprintf(['\n Displaying the results...... ' '\n'])

[tp,fp,tn,fn,fplv,fnlv,~,~,pcc,kappa,imw]=performance(CM,Ref_gt);
F1 = 2*tp/(2*tp + fp + fn);
result = 'PCC is %4.3f; KC is %4.3f; F1 is %4.3f \n';
fprintf('\n');fprintf(result,pcc,kappa,F1)

figure;
if strcmp(dataset,'#2-Texas-L8') == 1
    subplot(131);imshow(2*uint16(RegImg(:,:,[7 5 4])));title('Regression image')
elseif strcmp(dataset,'#6-California') == 1
    subplot(131);imshow(RegImg);title('Regression image')
elseif strcmp(dataset,'#8-Img5') == 1
    subplot(131);imshow(uint8(exp(RegImg*3.75+1.8)));title('Regression image')
elseif strcmp(dataset,'#9-Texas-ALI') == 1
    subplot(131);imshow(6*uint16(RegImg));title('Regression image')        
else
    subplot(131);imshow(uint8(RegImg));title('Regression image')
end
subplot(132);imshow(remove_outlier(DI),[]);title('Difference image')
subplot(133);imshow(CM,[]);title('Change mape')

if F1 < 0.3
   fprintf('\n');disp('Please exchange the order of the input images OR select the appropriate opt.alfa for SDA-HCD!')
end   

