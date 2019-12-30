%%this code has been adapted from: https://github.com/jason-holloway/towardCCA
clear all;
clc 
close all force

%% make sure the functions are located on MATLAB's path
setupPtych;
addpath('data')
addpath('Utils')

%% check to see which input parameters have been provided
if ~exist('dataset','var') || isempty(dataset)
    dataset = 'resChart2'; %'block'
end
if ~exist('apDia','var') || isempty(apDia)
    apDia = 57.5;
end
if ~exist('overlap','var') || isempty(overlap)
    %overlap = 0.12;
    overlap = 0.72;
end
if ~exist('N','var') || isempty(N)
    N = 9;
end
if ~exist('SNR','var') || isempty(SNR)
    SNR = 30;
end
if ~exist('nIts','var') || isempty(nIts)
    nIts = 10;
end

%% load the ground truth image
try
    load([dataset '.mat'],'im');
catch
    error('Dataset must contain a variable called ''im''');
end
if strcmp(dataset,'MITlogo')
    im = imcomplement(im);
end
% convert to floating point (use singles to save memory)
im = im2single(im); 
if ~ismatrix(im) % provided images are grayscale, just to double check
    im = rgb2gray(im);
end

[h,w] = size(im);
gt = im;

%% determine the spacing between adjacent apertures (in pixels)
spacing = apDia * (1-overlap);

%% sampling setup params
opts = struct();
opts.imHeight = h;
opts.imWidth = w; 
opts.nX = N; 
opts.nY = N;
opts.apertureShift = spacing; 
opts.apDia = apDia;
opts.pupilType = 'circle';
subsampling = 'randpix'; % 'randcam' for random cameras, 'randpix' for random pixels

if strcmp(subsampling,'randpix')
    opts.samplingPattern = ones(opts.nX,opts.nY); %pattern 1
elseif strcmp (subsampling,'randcam')
%% random camera subsampling
opts.samplingPattern = randi([0 1],opts.nY,opts.nX) ; %pattern 2, not to be used with random pixel sampling
opts.samplingPattern(ceil(opts.nX*opts.nY/2)) = opts.samplingPattern(ceil(opts.nX*opts.nY/2)) || 1; %pattern 2
end

%% create the observed images
fprintf('Creating the input data cube\n');
[y pupil samplingIndices] = forwardModel(im,opts); % y is the squared magnitude

%% random pixel subsampling
m = h*w*N*N;
f = 0.6;%0.02 for block.mat; %fraction of samples to be used
if f==1 || strcmp(subsampling,'randcam')
    Num = nnz(opts.samplingPattern);
    Cen = nnz(opts.samplingPattern(1:ceil(opts.nX*opts.nY/2)));    
    P_op = ones(h,w,Num);
    y_sub = y;
    f=1;
elseif f<1 || strcmp(subsampling,'randpix')
    fn = ceil(f*m);
    yvec = y(:);
    ind = randperm(m,fn);
    y_subvec = zeros(m,1); y_subvec(ind) = yvec(ind); y_sub = reshape(y_subvec,[h w N*N]);
    P_op = zeros(h,w,N*N); P_op(ind) = 1;
    Cen = ceil(N*N/2);
end

%% noise
fprintf('Adding noise\n');
% add noise
if ~isinf(SNR)
    y = addNoise(y,SNR);
end
y(y<0)=0; % input cannot be negative (avoid noise causing a negative signal)

%% recover the high resolution image
nn = (h+floor(spacing*(N-1)))*(w+floor(spacing*(N-1))); 
n = h*w; 
sf = 0.25; %0.1 for block.mat %estimated sparsity
s = floor(sf*nn);

%% 
fprintf('Recovering high resolution image - spatial sparsity\n');

[recovS, initS] = ptychMain_sparse(y_sub,'spatial',apDia,spacing,nIts,opts.samplingPattern,pupil,samplingIndices,P_op,s);

%% 
fprintf('Recovering high resolution image - block sparsity\n');

[recovB, initB] = ptychMain_sparse(y_sub,'block',apDia,spacing,nIts,opts.samplingPattern,pupil,samplingIndices,P_op,s);

%% 
fprintf('Recovering high resolution image - fourier sparsity\n');

[recovF, initF] = ptychMain_sparse(y_sub,'fourier',apDia,spacing,nIts,opts.samplingPattern,pupil,samplingIndices,P_op,s);


%% 
fprintf('Recovering high resolution image - no model\n');
model = 'fourier'; %fourier (default) or spatial
[recovAM initAM] = ptychMain(y_sub,model,apDia,spacing,nIts,opts.samplingPattern,pupil,samplingIndices,P_op);

%% display the ground truth, input, and recovered magnitudes
nn = sqrt(nn);
lb = (nn-w)/2+1;
ub = lb-1+w;

inputim = y_sub(:,:,Cen); % extract the center image
inputim = double(sqrt(inputim)); % display the magnitude of the observation (not the squared magnitude)

[dispRecovS, dispInitS, dispRecovFS, ssimmS] = display_params(gt,inputim,recovS,initS,n,lb,ub,'spatial','sparsity');
[dispRecov_F, dispInitF, dispRecovFF, ssimmF] = display_params(gt,inputim,recovF,initF,n,lb,ub,'fourier','sparsity');
[dispRecovB, dispInitB, dispRecovFB, ssimmB] = display_params(gt,inputim,recovB,initB,n,lb,ub,'block','block sparsity');
[dispRecovAM, dispInitAM, dispRecovFAM, ssimmAM] = display_params(gt,inputim,recovAM,initAM,n,lb,ub,model,'none');

display_result(gt,inputim,dispRecovS,dispRecovFS,'CoPRAM',ssimmS,sf,'spatial')
display_result(gt,inputim,dispRecov_F,dispRecovFF,'CoPRAM',ssimmF,sf,'fourier')
display_result(gt,inputim,dispRecovB,dispRecovFB,'Block CoPRAM',ssimmB,sf,'spatial')
display_result(gt,inputim,dispRecovAM,dispRecovFAM,'AltMin',ssimmAM)
