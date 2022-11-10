load(dataset)
%% show datasets
figure
if strcmp(dataset,'#2-Texas-L8') == 1
    subplot(131);imshow(2*uint8(image_t1(:,:,[4 3 2])));title('imaget1')
    subplot(132);imshow(2*uint16(image_t2(:,:,[7 5 4])));title('imaget2')
elseif strcmp(dataset,'#6-California') == 1
    subplot(131);imshow(image_t1(:,:,[4 3 2])+1,[]);title('imaget1')
    subplot(132);imshow(image_t2,[-1 1]);title('imaget2')
elseif strcmp(dataset,'#9-Texas-ALI') == 1
    subplot(131);imshow(2*uint8(image_t1));title('imaget1')
    subplot(132);imshow(6*uint16(image_t2));title('imaget2')
else
    subplot(131);imshow(uint8(image_t1));title('imaget1')
    subplot(132);imshow(uint8(image_t2));title('imaget2')
end
subplot(133);imshow(Ref_gt,[]);title('Refgt')
%% Normalization
opt.type_t1 = 'optical';
opt.type_t2 = 'optical';
if strcmp(dataset,'#5-Shuguang') == 1
    image_t1 = remove_outlier(double(image_t1)); % for SAR image of Shuguang dataset
elseif strcmp(dataset,'#7-YellowRiver') == 1 % Heterogeneous CD of SAR VS. Optical
    opt.type_t1 = 'sar';
    opt.type_t2 = 'optical';
elseif strcmp(dataset,'#8-Img5') == 1 % Heterogeneous CD of Optical VS. SAR
    opt.type_t1 = 'optical';
    opt.type_t2 = 'sar';
    image_t1 = image_t1(1:4:end,1:4:end,:);
    image_t2 = image_t2(1:4:end,1:4:end,:);
    Ref_gt = Ref_gt(1:4:end,1:4:end);
end
image_t1 = image_normlized(image_t1,opt.type_t1);
image_t2 = image_normlized(image_t2,opt.type_t2);


