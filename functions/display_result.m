%% displays results
function [] = display_results(gt,inputim,dispRecov,dispRecovF,method,ssimm,sf,basis)

%%
h2 = figure;
set(h2,'name',sprintf('Simulation of %s',dataset),'numbertitle','off');
set(h2,'units','normalized','OuterPosition',[0 0 1 1]);
subplot(221)
imagesc(gt), colormap(gray), axis image off
title('Ground truth image')
subplot(222)
imagesc(abs(inputim)), colormap(gray), axis image off
title('Center image')
subplot(223)
imagesc(dispRecov), colormap(gray), axis image off
if nargin > 6
    tt = [method,': Recovered image: SSIM - ',num2str(ssimm),', ',basis,' sparsity s=',num2str(sf),'n'];
else
tt = [method,': Recovered image: SSIM - ',num2str(ssimm)];
end
title(tt)
subplot(224)
imagesc(dispRecovF), colormap(gray), axis image off
title('Fourier magnitude of recovered image')

end