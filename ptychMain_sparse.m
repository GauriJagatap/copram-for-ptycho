function [x, xInit] = ptychMain_sparse(y,basis,apDia,spacing,nIts,samplingPattern,pupil,samplingIndices,P_op,s)

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

hROW = h+floor(spacing*(nX-1));
hCOL = w+floor(spacing*(nY-1));


sl = ceil(s*(h*w)/(hROW*hCOL));
if strcmp(basis,'block')
    Jb = 4;%48; %breadth of block for padded image, set to 48 for block.mat
    Jh = 4;%48; %height of block
    J = Jb*Jh;
end

%% initialize and set up the annoyomous helper functions
switch basis
case 'fourier'
    % get initial estimate
    lowresMean = sqrt(mean(mean(y.^2,3),4));
    center = fftshift(fft2(lowresMean))/sqrt(h*w);
    %%
    xInit = padarray(center,floor([(hROW-h)/2 (hCOL-w)/2]));
    f_h = @(z) F_LENS2SENSOR(z,samplingIndices,pupil,h,w);
    f_hp = @(z)P_op.*f_h(z);
    ft_h = @(z) F_SENSOR2LENS(z,samplingIndices,hROW,hCOL,h,w,pupil,basis);
    ft_hp = @(z) ft_h(P_op.*z);
    % intialize x
    x = xInit;
case {'spatial','block'}  
    % get initial estimate
    lowresMean = sqrt(mean(mean(y.^2,3),4));
    %%
    center = fftshift(fft2(lowresMean))/sqrt(h*w);
    xInit = padarray(center,floor([(hROW-h)/2 (hCOL-w)/2]));
    xInit = sqrt(hROW*hCOL)*ifft2(ifftshift(xInit));
    f_h = @(z) F_LENS2SENSOR(z,samplingIndices,pupil,h,w);
    f_hp = @(z)P_op.*f_h(fftshift(fft2(z))/sqrt(h*w));
    ft_h = @(z) F_SENSOR2LENS(z,samplingIndices,hROW,hCOL,h,w,pupil,basis);
    ft_hp = @(z) ft_h(P_op.*z);
    % intialize x
    x = xInit;
end

% set up the display window
f1 = figure('numbertitle','off','units','normalized','outerposition',[0 0 1 1]);
drawnow;

for ii = 1:nIts
    
    if mod(ii,2)==1 || ii==nIts
        fprintf('N: %02d Iteration: %04d\n',nX,ii);
    end
    
    % compute y = A*x
    y0 = f_hp(x); 
    
    % enforce magnitude measurements
    y0 = y0./abs(y0+eps).*y;
    
    % update x estimate: x = A'*y
    nn = hROW*hCOL;
    x_hist = x;
    
    switch basis
        case {'spatial','fourier'}
            x = cosamp_fun(y0, f_hp, ft_hp, nn, s, 10);
        case 'block'
            x = jsmp_fun(y0, f_hp, ft_hp, nn, floor(s/(Jb*Jh)),Jh,Jb, 10); %height x breadth of block
    end
    
    rel_err = 1e-4;
    if norm(x-x_hist)/norm(x) < rel_err
        fprintf('met relative error tolerance condition\n');
        break;
    end
    if mod(ii,2)==1 || ii==nIts
        set(f1,'name',sprintf('ApDia: %02d Spacing: %02.2f Iteration: %04d',apDia,spacing,ii));
        subplot(2,2,1)
        switch basis
        case 'fourier'
            imagesc(log10(abs(x))); colormap(gray); axis image
            title('recovered Fourier magnitude');
            recov = sqrt(hCOL*hROW)*ifft2(ifftshift(x));
            subplot(2,2,[2 4]);
            imagesc((abs(recov))) 
            title('spatial magnitude'); axis image;
            subplot(2,2,3);
            imagesc(angle(recov)); axis image      
        case {'spatial','block'}
            imagesc(log10(abs(fftshift(fft2(x))/sqrt(hCOL*hROW)))); colormap(gray); axis image
            title('recovered Fourier magnitude');
            recov = x;
            subplot(2,2,[2 4]);
            imagesc((abs(recov))) 
            title('spatial magnitude'); axis image;
            subplot(2,2,3);
            imagesc(angle(fftshift(fft2(recov)))); axis image    
        end       
        title('spatial phase')
        drawnow;
    end
end

end