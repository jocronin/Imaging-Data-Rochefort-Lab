% Get Data

%% PV cre 106,117,121,122 (recs,days)

% 106
clear
days = ['08';'09';'10';'11';'12'];
Data_pvcre106 = cell(5,34);


for d=1:length(days)

files = dir(['/Users/josephcronin/Documents/Imaging Data Rochefort Lab/201412',days(d,:),'_PVcre106/CaActivity']);

nn=[];
for i=3:length(files)
    
n=files(i,1).name;

t=findstr('deltaf',n);

tt = isempty(t);

if (tt==1)
    
else
    
nn = [nn; n];    
    
end

end

length(nn)

recs = nn;


for j=1:length(recs)
Data_pvcre106{d,j}=load(['/../201412',days(d,:),'_PVcre106/CaActivity/',nn(j,:)]);
end

end


% 117
days = ['09';'10';'11';'12';'14';'21'];
Data_pvcre117 = cell(6,46);


for d=1:length(days)

files = dir(['/Users/josephcronin/Documents/Imaging Data Rochefort Lab/201502',days(d,:),'_PVcre117/CaActivity']);

nn=[];
for i=3:length(files)
    
n=files(i,1).name;

t=findstr('deltaf',n);

tt = isempty(t);

if (tt==1)
    
else
    
nn = [nn; n];    
    
end

end

length(nn)

recs = nn;


for j=1:length(recs)
Data_pvcre117{d,j}=load(['/../201502',days(d,:),'_PVcre117/CaActivity/',nn(j,:)]);
end

end


% 121
days = ['18';'19';'20';'21';'22';'23';'09'];
Data_pvcre121 = cell(7,43);


for d=1:length(days)

if (d~=7)
files = dir(['/Users/josephcronin/Documents/Imaging Data Rochefort Lab/201502',days(d,:),'_PVcre121/CaActivity']);
else
files = dir(['/Users/josephcronin/Documents/Imaging Data Rochefort Lab/201503',days(d,:),'_PVcre121/CaActivity']);   
end

nn=[];
for i=3:length(files)
    
n=files(i,1).name;

t=findstr('deltaf',n);

tt = isempty(t);

if (tt==1)
    
else
    
nn = [nn; n];    
    
end

end

length(nn)

recs = nn;


for j=1:length(recs)
if (d~=7)    
Data_pvcre121{d,j}=load(['/../201502',days(d,:),'_PVcre121/CaActivity/',nn(j,:)]);
else
Data_pvcre121{d,j}=load(['/../201503',days(d,:),'_PVcre121/CaActivity/',nn(j,:)]);   
end

end

end


% 122
days = ['07';'08';'09';'10';'11';'12';'19'];
Data_pvcre122 = cell(7,87);


for d=1:length(days)

files = dir(['/Users/josephcronin/Documents/Imaging Data Rochefort Lab/201503',days(d,:),'_PVcre122/CaActivity']);

nn=[];
for i=3:length(files)
    
n=files(i,1).name;

t=findstr('deltaf',n);

tt = isempty(t);

if (tt==1)
    
else
    
nn = [nn; n];    
    
end

end

length(nn)

recs = nn;


for j=1:length(recs)
Data_pvcre122{d,j}=load(['/../201503',days(d,:),'_PVcre122/CaActivity/',nn(j,:)]);
end

end

% plot some

for j = 1: 10
figure(j)
for i =1:16

   subplot(4,4,i)
   plot(Data_pvcre117{1,i}.deltaf(j,:))
    
end
end

% peak det / filter

V = Data_pvcre117{1,i}.deltaf(1,:);
X = linspace(1,length(V),length(V));

Vf=filter([1/4 1/4 1/4 1/4 1/4],1,V);
Xf = linspace(1,length(Vf),length(Vf));


[ma, mi] = peakdet(V,0.15,X);
[maf, mif] = peakdet(Vf,0.15,Xf);

subplot(2,1,1)
plot(ma(:,1),ma(:,2),'ro',X,V,'b')
subplot(2,1,2)
plot(maf(:,1),maf(:,2),'ro',Xf,Vf,'b')

% fft filter

[Vft, f, y, y2] = fftf(X,V);
Xft = linspace(1,length(Vft),length(Vft));

[ma, mi] = peakdet(V,0.15,X);
[maft, mift] = peakdet(Vft,0.15,Xft);

subplot(2,1,1)
plot(ma(:,1),ma(:,2),'ro',X,V,'b')
subplot(2,1,2)
plot(maft(:,1),maft(:,2),'ro',Xft,Vft,'b')

% for example take PVcre117 0210 (mostly dark with stim in between)

VV=[];
for i = 1 : 30
Vv = Data_pvcre117{2,i}.deltaf;
VV=[VV; Vv];
end
X = linspace(1,length(VV),length(VV));

for k = 1:100
    
   subplot(10,10,k)
   plot(VV(k,:))
    
end

% are the channels similar across recordings???

for k = 1:16
    
   subplot(4,4,k)
   plot(VV(6+74*(k-1),:))
    
end

% average them
AA = zeros(74,2880);

for m = 1:74
for n = 1:30
    
  AA(m,:) = AA(m,:) + VV(m+74*(n-1),:);
    
end
end

AA = AA / 30;

for k = 1:64
    
   subplot(8,8,k)
   plot(AA(k,:))
    
end

% maybe not...try and catononante all same chanel recordings

AA = [];

for m = 1:74
  Aa = [];
for n = 1:30
    
  Aa = [Aa VV(m+74*(n-1),:)];
    
end
AA = [AA; Aa];

end
X=linspace(1,length(AA),length(AA));

for k = 3:10
    
   subplot(2,4,k-2)
   plot(AA(k,:))
    
end

% for a PVcre117 0210 filter A03 - A10 (nicest looking ones

VVft=[];
for i=1:10
    
[Vft, f, y, y2] = fftf(X,AA(i,:));
Xft = linspace(1,length(Vft),length(Vft));

VVft = [VVft; Vft];

end

for k = 1:8
    
   subplot(2,4,k)
   plot(VVft(k,:))
    
end

spi=[];
thresh=0.1;

[maft, mift] = peakdet(VVft(1,:),thresh,X);
[maft, ift1] = peakdet(VVft(1,:),thresh);

subplot(2,4,1)
plot(maft(:,1),maft(:,2),'ro',X,VVft(1,:),'b')

[maft, mift] = peakdet(VVft(2,:),thresh,X);
[maft, ift2] = peakdet(VVft(2,:),thresh);

subplot(2,4,2)
plot(maft(:,1),maft(:,2),'ro',X,VVft(2,:),'b')

[maft, mift] = peakdet(VVft(3,:),thresh,X);
[maft, ift3] = peakdet(VVft(3,:),thresh);
subplot(2,4,3)
plot(maft(:,1),maft(:,2),'ro',X,VVft(3,:),'b')

[maft, mift] = peakdet(VVft(4,:),thresh,X);
[maft, ift4] = peakdet(VVft(4,:),thresh);
subplot(2,4,4)
plot(maft(:,1),maft(:,2),'ro',X,VVft(4,:),'b')

[maft, mift] = peakdet(VVft(5,:),thresh,X);
[maft, ift5] = peakdet(VVft(5,:),thresh);

subplot(2,4,5)
plot(maft(:,1),maft(:,2),'ro',X,VVft(5,:),'b')

[maft, mift] = peakdet(VVft(6,:),thresh,X);
[maft, ift6] = peakdet(VVft(6,:),thresh);

subplot(2,4,6)
plot(maft(:,1),maft(:,2),'ro',X,VVft(6,:),'b')

[maft, mift] = peakdet(VVft(7,:),thresh,X);
[maft, ift7] = peakdet(VVft(7,:),thresh);
subplot(2,4,7)
plot(maft(:,1),maft(:,2),'ro',X,VVft(7,:),'b')

[maft, mift] = peakdet(VVft(8,:),thresh,X);
[maft, ift8] = peakdet(VVft(8,:),thresh);
subplot(2,4,8)
plot(maft(:,1),maft(:,2),'ro',X,VVft(8,:),'b')

% convert to bin matrixes

    
spii = zeros(1,86400);
spii(ift1(:,1))=1;
spi=[spi; spii];

spii = zeros(1,86400);
spii(ift2(:,1))=1;
spi=[spi; spii];

spii = zeros(1,86400);
spii(ift3(:,1))=1;
spi=[spi; spii];

spii = zeros(1,86400);
spii(ift4(:,1))=1;
spi=[spi; spii];

spii = zeros(1,86400);
spii(ift5(:,1))=1;
spi=[spi; spii];

spii = zeros(1,86400);
spii(ift6(:,1))=1;
spi=[spi; spii];

spii = zeros(1,86400);
spii(ift7(:,1))=1;
spi=[spi; spii];

spii = zeros(1,86400);
spii(ift8(:,1))=1;
spi=[spi; spii];

sum(sum(spi))

hold on
for i =1:8 
    
plot(spi(i,:)*i,'+')    
    
end
hold off