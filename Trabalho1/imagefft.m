function imagefft(I)
    F=fft2(double(I));
    S=fftshift(F);
    L=log2(S);
    A=abs(L);
    imagesc(A)
end