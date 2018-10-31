function y = F_LENS2SENSOR(z,indices,pupil,ROW,COL)
z = repmat(z,1,1,size(indices,3));
x = z(indices);
 %induces phase shift on y, but we don't care about the phase of y
if isequal(size(x),size(pupil))
    y = ifft2(x.*pupil);
else
    y = ifft2(bsxfun(@times,x,pupil));
end
y = y*sqrt(ROW*COL);
end