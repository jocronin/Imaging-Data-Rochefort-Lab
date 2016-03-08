% denoise

SS_dn=[];

for i=1:74
    
S_dn = func_denoise_dw1d(AA(i,:));
X_dn = linspace(1,length(S_dn),length(S_dn));

SS_dn = [SS_dn; S_dn];

end

for k = 1:10
    
   subplot(2,5,k)
   plot(SS_dn(k,:))
    
end