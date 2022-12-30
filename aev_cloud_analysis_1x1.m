% Matlab script to:
%   [1] upload raw Poisson tiff images;
%   [2] correct bg by imopen;
%   [3] threshold by Otsu;
%   [4] centroid estimation
% Andreas E. Vasdekis, University of Idaho, 6/18/2021

clear; clc;

tStart = cputime;

% path for uploading images
srcFiles = dir('E:\ls\5_P_R\4_CL_P_comp\5_mag\220218_2x2\2_output\PS_200nm_63x\*.tif'); 
imin_dir = 'E:\ls\5_P_R\4_CL_P_comp\5_mag\220218_2x2\2_output\PS_200nm_63x\';

% path for saving output images
figdir_out = 'E:\ls\5_P_R\4_CL_P_comp\5_mag\220218_2x2\2_output\PS_200nm_63x\1x1\';    

se = strel('disk', 1);       % morphological operation with 1 pixel radius

for i = 1 : length(srcFiles)

filename = strcat(imin_dir, srcFiles(i).name);
FileTif1 = filename;
InfoImage = imfinfo(FileTif1);
mImage = InfoImage(1).Width;
nImage = InfoImage(1).Height;
z_planes = length(InfoImage);

m1 = zeros(nImage, mImage, z_planes, 'uint16');
m2 = zeros(nImage, mImage, z_planes, 'uint16');
m3 = zeros(nImage, mImage, z_planes, 'uint16');
m3_2x2 = zeros(2*nImage, 2*mImage, z_planes, 'uint16');

TT = Tiff(FileTif1, 'r+');
TT.setDirectory(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTSU adaptive thresholding 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for ii = 1:z_planes
        disp(ii);
        TT.setDirectory(ii);
        m1(:,:,ii) = TT.read();
        m1(:,:,ii) = double(imbinarize(m1(:,:,ii),'adaptive','Sensitivity',0.2));
        m1(:,:,ii) = imopen(m1(:,:,ii), se);
        D = -bwdist(~m1(:,:,ii)); 
        L = watershed(D);
        L(~m1(:,:,ii)) = -Inf; 
        m2(:,:,ii) = L;
        clear D L
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% centroid estimation and output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for ii = 1:z_planes
        disp(ii);
        c(ii) = bwconncomp(m2(:,:,ii));
        sc = regionprops(m2(:,:,ii),'centroid');
        centroids = round(cat(1,sc.Centroid),0);
        [c1, ~] = size(centroids);
        if c1 > 0
            for x = 1:c1
                A = centroids(x,:);
                m3(A(1), A(2), ii) = 1;
            end
        end

        clear c1 c2 x y A m4
    end

IM = sum(m3,3);

% loop to save the stack of increasing <N>

    for ii = 1:z_planes
        disp('final step');
        disp(ii);
        m4(:,:,1:ii) = m3(:,:,1:ii);
        m5 = sum(m4,3);   
        IMi(:,:,ii) = m5;

        clear m4 m5
    end

IMi = imrotate(IMi,90,'bilinear','crop');
  
% saving the output image
baseFileName =  sprintf(srcFiles(i).name, 1);  
fullFileName = fullfile(figdir_out, baseFileName);
clear options;
saveastiff(IMi, fullFileName);

clearvars -except srcFiles imin_dir figdir_out i se tStart centroids

end

disp(cputime - tStart);