function [x xInit] = ptychMain(y,apDia,spacing,nIts,samplingPattern,pupil,P_op,tau)

% input images are the squared magnitudes, take the square root
y = sqrt(y);

% check to make sure the Fourier domain sampling is known
if ~exist('samplingPattern','var') || isempty(samplingPattern)
    if ndims(y)~=4
        error('If no sampling pattern is specified, the input images must be 4 dimensional');
    end
    
    % all images are on a square grid
    [~,~,nY,nX] = size(y);
    samplingPattern = ones(nY,nX); 
end

[nY,nX] = size(samplingPattern);
[h,w,~] = size(y);
y = reshape(y,h,w,[]);

% make sure the input is square
if h<w
    pad = (w-h)/2;
    padTop = floor(pad);
    padBottom = ceil(pad);
    y = padarray(y,[padTop,0,0],'pre');
    y = padarray(y,[padBottom,0,0],'post');
elseif w<h
    pad = (h-w)/2;
    padLeft = floor(pad);
    padRight = ceil(pad);
    y = padarray(y,[0,padLeft,0],'pre');
    y = padarray(y,[0,padRight,0],'post');
end
[h,w,~] = size(y);

% get the Fourier domain sampling pattern
opts.imHeight = h; opts.imWidth = w; opts.nX = nX; opts.nY = nY;
opts.samplingPattern = samplingPattern;
opts.apertureShift = spacing; opts.apDia = apDia;
opts.pupilType = 'circle';
[samplingIndices,pupil,hROW,hCOL] = getSampling(opts);

%% new
N = nX;
m = h*w*N*N;
f = 1;%0.02 for block.mat; %fraction of samples to be used
subsampling = 'randpix';
if f==1 || strcmp(subsampling,'randcam')
    Num = nnz(opts.samplingPattern);
    Cen = nnz(opts.samplingPattern(1:ceil(opts.nX*opts.nY/2)));
    P_op = ones(h,w,Num);
    y_sub = y;
    f = 1;
elseif f<1 || strcmp(subsampling,'randpix')
    fn = ceil(f*m);
    yvec = y(:);
    ind = randperm(m,fn);
    y_subvec = zeros(m,1); y_subvec(ind) = yvec(ind); y_sub = reshape(y_subvec,[h w N*N]);
    P_op = zeros(h,w,N*N); P_op(ind) = 1;
    Cen = ceil(N*N/2);
end
y = y_sub;
%% 

% get initial estimate
lowresMean = sqrt(mean(mean(y.^2,3),4));
center = fftshift(fft2(lowresMean))/sqrt(h*w);
% centerscale = sqrt(nX*nY)*center; %scale to approximately maintain overall norm (this will work best for 0.2 overlap)
centerscale = center; %scale to approximately maintain overall norm (this will work best for 0.2 overlap)
xInit = padarray(centerscale,floor([(hROW-h)/2 (hCOL-w)/2]));

% centerView = imresize(lowresMean,[hROW hCOL],'bilinear');
% centerView = centerView/norm(centerView(:));
% xInit = fftshift(fft2(centerView))/sqrt(h*w); %this is like signal proxy
% xInit = imresize(xInit,[hROW hCOL],'bilinear'); %dunno why this is required

% set up the annoyomous helper functions
f_h = @(z) F_LENS2SENSOR((1/sqrt(h*w))*fftshift(fft2(z)),samplingIndices,pupil,h,w);
f_hp = @(z) P_op.*f_h(z);
ft_h = @(z) F_SENSOR2LENS(z,samplingIndices,hROW,hCOL,h,w,pupil,'spatial');
ft_hp = @(z) sqrt(hROW*hCOL)*ifft2(ifftshift(ft_h(P_op.*z)));
% intialize x
x = sqrt(hROW*hCOL)*ifft2(ifftshift(xInit));

% regularization parameter to ensure no dividing by zero
if ~exist('tau','var') || isempty(tau)
    tau = 2;
end

AtA = real(ft_hp(f_hp(ones(size(x)))))+tau;
iAtA = tau./AtA;

% set up the display window
f1 = figure('numbertitle','off','units','normalized','outerposition',[0 0 1 1]);
drawnow;
nn=hROW*hCOL;
s = round(0.25*nn);
for ii = 1:nIts
    
    if mod(ii,2)==1 || ii==nIts
        fprintf('N: %02d Iteration: %04d\n',opts.nX,ii);
    end
    
    % compute y = A*x
    y0 = f_hp(x); 
    
    % enforce magnitude measurements
    y0 = y0./abs(y0+eps).*y;
    
    % update x estimate: x = A'*y
    %x = iAtA.*ft_hp(y0); 
    x = cosamp_fun(y0, f_hp, ft_hp, nn, s, 10);

    if mod(ii,2)==1 || ii==nIts
        set(f1,'name',sprintf('ApDia: %02d Spacing: %02.2f Iteration: %04d',apDia,spacing,ii));
        subplot(2,2,1)
        imagesc(log10(abs(x))); colormap(gray); axis image
        title('recovered Fourier magnitude');
        
        recov = sqrt(hCOL*hROW)*ifft2(ifftshift(x));
        subplot(2,2,[2 4]);
        imagesc((abs(recov))) 
        title('spatial magntude'); axis image;
        subplot(2,2,3);
        imagesc(angle(recov)); axis image
        title('spatial phase')
        drawnow;
    end
end

end