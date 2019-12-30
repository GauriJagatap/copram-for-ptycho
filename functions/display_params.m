function [dispRecov, dispInit, dispRecovF, ssimm] = display_params(gt,inputim,recov,init,n,lb,ub,basis,expt)

    switch basis
        case {'fourier','none'}
            dispRecov = sqrt(n)*abs(ifft2(ifftshift(recov(lb:ub,lb:ub))));
            dispInit = sqrt(n)*abs(ifft2(ifftshift(init(lb:ub,lb:ub))));
            dispRecovF = log10(abs(recov(lb:ub,lb:ub)));
        case {'spatial','block'}
            fftrecov = fftshift(fft2(recov))/sqrt(n);
            fftrim = fftrecov(lb:ub,lb:ub);
            dispRecovF = log10(abs(fftrim));
            recovtrim = sqrt(n)*ifft2(ifftshift(fftrim));
            dispRecov = abs(recovtrim);
            dispInit = abs(init);
    end
    
    %%
    dgt = double(gt);
    ssimmC = ssim(inputim,dgt);
    ssimm = ssim(double(dispRecov),dgt);

    fp = ['SSIM - basis:',basis,', model:',expt,':%2.8f\n'];
    fprintf(fp,ssimm);
end