function sigma=fsigma(x,expon,Xsigma)

sigma=0;
for i=1:length(expon)
    sigma=sigma+Xsigma(i)*x.^expon(i);
end
