% Spetrogram stuff
figure
for i = 1:10
Nx = length(AA(i,:));
nsc = floor(Nx/4.5);
wind=hamming(nsc);
nov = floor(nsc/2);
nff = max(256,2^nextpow2(nsc));
subplot(2,5,i)
spectrogram(AA(i,:),wind,nov,nff,'yaxis') % 1528 time bins
end

spectrogram(AA(3,:),300,100,328,'yaxis') % 1528 time bins
ff=psdfreqvec(AA(3,:));

figure
for i = 1:10
subplot(2,5,i)
spectrogram(VVft(i,:),'yaxis')
end

s=3;
SS = sgolayfilt(AA(s,:),1,41);
subplot(3,1,1)
plot(SS)
subplot(3,1,2)
plot(VVft(s,:))
subplot(3,1,3)
plot(AA(s,:))


[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('db1');
[WW, L] = wavedec(AA(s,:),3,'db1');
plot(WW)

% denoising with wavlets

% save data file 

dat=AA(3,:);
save('test.mat','dat')