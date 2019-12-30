function [recov,input] = realDataDemo(dataset)
%REALDATADEMO example code for Fourier Ptychography
%   [recov,input] = realDataDemo(dataset)
%   recovers a high resolution image from a series of low resolution inputs
%   saved in dataset.mat. The recovery parameters have been hardcoded to
%   provide the results shown in the paper.
%   NOTE: For real data, the squared magnitude of the recovered image is
%   shown as the recorded input image is also the squared magnitude.
% 
%   - Outputs
%   recov - The recovered m x m complex field, take the inverse fourier
%   transform and view the squared magnitude to compare to the input images
%   
%   input - This is the central view of the N x N sampling grid
% 
%   - Inputs (all necessary parameters have been hard coded)
%   dataset - a string specifying which dataset to use ('USAF',
%   'fingerprint', or 'dasani'). Default ['fingerprint'].

clc, close all force

% make sure the functions are located on MATLAB's path
setupPtych;


if ~exist('dataset','var') || isempty(dataset)
    dataset = 'realData/USAF';
end

fprintf('Loading data\n');
data = load([dataset '.mat']);

% get the HDR images
if ndims(data.ims)==5
    fprintf('Forming the HDR images\n');
    data.ims = createHDR(data.ims);
end

y= data.ims;
[h,w,nX,nY] = size(y);
subsampling = 'randpix'; % 'randcam' for random cameras, 'randpix' for random pixels
if strcmp(subsampling,'randpix')
    samplingPattern = ones(nX,nY); %pattern 1
elseif strcmp (subsampling,'randcam')
%% random camera subsampling
    samplingPattern = randi([0 1],nY,nX) ; %pattern 2, not to be used with random pixel sampling
	samplingPattern(ceil(nX*nY/2)) = samplingPattern(ceil(nX*nY/2)) || 1; %pattern 2
end

%% random pixel subsampling
m = h*w*nX*nY;
f = 0.1;%0.02 for block.mat; %fraction of samples to be used
if f==1 || strcmp(subsampling,'randcam')
    Num = nnz(samplingPattern);
    Cen = nnz(samplingPattern(1:ceil(nX*nY/2)));    
    P_op = ones(h,w,Num);
    y_sub = y;
    f=1;
elseif f<1 || strcmp(subsampling,'randpix')
    fn = ceil(f*m);
    yvec = y(:);
    ind = randperm(m,fn);
    y_subvec = zeros(m,1); y_subvec(ind) = yvec(ind); y_sub = reshape(y_subvec,[h w nX*nY]);
    P_op = zeros(h,w,nX*nY); P_op(ind) = 1;
    Cen = ceil(nY*nX/2);
end

fprintf('Recovering the high resolution image\n');
recov = ptychMain(y_sub,'fourier',data.apDia,data.spacing,data.nIts,samplingPattern,[],[],P_op,data.tau);

% compare the input center image and the recovered image
dispRecov = ifft2(ifftshift(recov));
dispRecov = abs(dispRecov).^2;

input = data.ims(:,:,ceil(end/2)); % extract the center image

h = figure(9);
set(h,'name',sprintf('Recovery of %s',dataset),'numbertitle','off');
set(h,'units','normalized','OuterPosition',[0 0 1 1]);
drawnow;
subplot(121)
imagesc(input), colormap(gray), axis image off
title('Center input image')
subplot(122)
imagesc(dispRecov), colormap(gray), axis image off
title('Recovered image')

end