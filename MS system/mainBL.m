clear;
clc;

zm=2;
zxmin=-zm;
zxmax=zm;
zymin=-zm;
zymax=zm;
Nz=1e7;
zxf0=rand(1,Nz)*2*zm-zm;
zyf0=rand(1,Nz)*2*zm-zm;

% N=1;                   %%% 一个初始点给出的样本数
h=0.001;
alphax=1.05;
alphay=1.5;
betax=0.5;
betay=-0.5;

%%% Generate data
Mx=stblrnd(alphax,betax,h^(1/alphax),0,1,Nz);
My=stblrnd(alphay,betay,h^(1/alphay),0,1,Nz);
Bh=sqrt(h)*randn(2,Nz);
zxf=(zxf0-zxf0.^3-5*zxf0.*zyf0.^2)*h+(1+zxf0).*Bh(1,:)+Bh(2,:)+(1-zxf0+zxf0.^2).*Mx;
zyf=-(1+zxf0.^2).*zyf0*h+zyf0.*Bh(2,:)+1./(1+0.5*zyf0.^2).*My;

%%% Identify alphax betax sigmax
zx0=[];
NI=10:25;
nklp=[];
nkln=[];
Ml=[];
q=0.5;
m=5;
N=1;
for ini=1:length(NI)
    Nint=NI(ini);
    zlin=linspace(zxmin,zxmax,Nint+1);
    z0=0.5*(zlin(1:end-1)+zlin(2:end));
    dz=0.5*(z0(2)-z0(1));
    zx0=[zx0,z0];
    for i=1:Nint
        I=(zxf0>=z0(i)-dz)&(zxf0<z0(i)+dz);
        zxfi0=zxf0(I);
        zxfi=zxf(I);
        Ml=[Ml;length(zxfi)];
        nkp=zeros(1,N+1);
        nkn=zeros(1,N+1);
        for k=0:N
            I=((zxfi)>=m^k*q)&((zxfi)<m^(k+1)*q);
            nkp(k+1)=length(zxfi(I));
            I=((zxfi)<=-m^k*q)&((zxfi)>-m^(k+1)*q);
            nkn(k+1)=length(zxfi(I));
        end
        nklp=[nklp;nkp];
        nkln=[nkln;nkn];
    end
end
nkl=nklp+nkln;
nratio=sum(nkl(:,1)./Ml)/sum(nkl(:,2)./Ml);
alphax0=log(nratio)/(log(m));
npnratio=(sum(nkln(:,1)./Ml)+sum(nkln(:,2)./Ml))/(sum(nklp(:,1)./Ml)+sum(nklp(:,2)./Ml));
betax0=(1-npnratio)/(1+npnratio);
kalphax=alphax0*(1-alphax0)/(cos(pi*alphax0/2)*gamma(2-alphax0));
sigmax0=q*((nkl(:,1)+nkl(:,2))./Ml*alphax0/(h*kalphax*(1-m^(-(N+1)*alphax0)))).^(1/alphax0);
sigmax0=sigmax0';
errorx=sqrt(mean((sigmax0-(1-zx0+zx0.^2)).^2))/max(abs(1-zx0+zx0.^2));
% figure;
% plot(zx0,sigmax0,'*');
% hold on
% plot(zx0,1-zx0+zx0.^2,'o');
% hold off

%%% Identify alphay betay sigmay
zy0=[];
nklp=[];
nkln=[];
Ml=[];
for ini=1:length(NI)
    Nint=NI(ini);
    zlin=linspace(zymin,zymax,Nint+1);
    z0=0.5*(zlin(1:end-1)+zlin(2:end));
    dz=0.5*(z0(2)-z0(1));
    zy0=[zy0,z0];
    for i=1:Nint
        I=(zyf0>=z0(i)-dz)&(zyf0<z0(i)+dz);
        zyfi0=zyf0(I);
        zyfi=zyf(I);
        Ml=[Ml;length(zyfi)];
        nkp=zeros(1,N+1);
        nkn=zeros(1,N+1);
        for k=0:N
            I=((zyfi)>=m^k*q)&((zyfi)<m^(k+1)*q);
            nkp(k+1)=length(zyfi(I));
            I=((zyfi)<=-m^k*q)&((zyfi)>-m^(k+1)*q);
            nkn(k+1)=length(zyfi(I));
        end
        nklp=[nklp;nkp];
        nkln=[nkln;nkn];
    end
end
nkl=nklp+nkln;
nratio=sum(nkl(:,1)./Ml)/sum(nkl(:,2)./Ml);
alphay0=log(nratio)/(log(m));
npnratio=(sum(nkln(:,1)./Ml)+sum(nkln(:,2)./Ml))/(sum(nklp(:,1)./Ml)+sum(nklp(:,2)./Ml));
betay0=(1-npnratio)/(1+npnratio);
kalphay=alphay0*(1-alphay0)/(cos(pi*alphay0/2)*gamma(2-alphay0));
sigmay0=q*((nkl(:,1)+nkl(:,2))./Ml*alphay0/(h*kalphay*(1-m^(-(N+1)*alphay0)))).^(1/alphay0);
sigmay0=sigmay0';
errory=sqrt(mean((sigmay0-1./(1+0.5*zy0.^2)).^2))/max(abs(1./(1+0.5*zy0.^2)));

% %%% Identify sigmax and sigmay by hard-thresholding sparse regression
% Ncoef=4;
% nz=length(zx0);
% A=zeros(nz,Ncoef+1);
% A(:,1)=1;
% for i=1:Ncoef
%     A(:,i+1)=zx0'.^i;
% end
% % sigmax
% Delta=0.05;
% A2=A;
% Bx=sigmax0';
% NA=length(A(1,:));
% posx=1:NA;
% for k=1:NA
%     X=(A'*A)\(A'*Bx);
%     I=abs(X)<Delta;
%     A(:,I)=[];
%     posx(I)=[];
%     if isempty(X(I))
%         break;
%     end
% end
% Xsigma=(A'*A)\(A'*sigmax0');
% exponx=posx-1;
% fxs=fsigma(zx0,exponx,Xsigma);
% figure;
% plot(zx0,sigmax0,'*');
% hold on
% plot(zx0,1-zx0+zx0.^2,'o');
% hold on
% plot(zx0,fxs,'r.');
% hold off
% % sigmay
% A=A2;
% By=sigmay0';
% NA=length(A(1,:));
% posy=1:NA;
% for k=1:NA
%     Y=(A'*A)\(A'*By);
%     I=abs(Y)<Delta;
%     A(:,I)=[];
%     posy(I)=[];
%     if isempty(Y(I))
%         break;
%     end
% end
% Ysigma=(A'*A)\(A'*sigmay0');
% expony=posy-1;
% fys=fsigma(zy0,expony,Ysigma);
% figure;
% plot(zy0,sigmay0,'*');
% hold on
% plot(zy0,1./(1+0.5*zy0.^2),'o');
% hold on
% plot(zy0,fys,'r.');
% hold off


%%% Identify sigmax and sigmay using SSR sparsing
Ncoef=4;
nz=length(zx0);
A=zeros(nz,Ncoef+1);
A(:,1)=1;
for i=1:Ncoef
    A(:,i+1)=zx0'.^i;
end
Xsigma0=(A'*A)\(A'*sigmax0');
% SSR
kfold=2;
nblock=floor(nz/kfold);
n=kfold*nblock;
IA=zeros(kfold,nblock);
position=1:n;
for i=1:kfold-1
    u=rand(1,nblock);
    J=kfold+1-i;
    p=0:J:J*(nblock-1);
    IA(i,:)=position(p+floor(J*u)+1);
    position(p+floor(J*u)+1)=[];
end
IA(kfold,:)=position;
% sigmax
A3=A;
A2=A;
B=sigmax0';
NA=length(A(1,:));
deleteordersigmax=zeros(1,NA);
deltaSSRsigmax=zeros(1,NA);
delta0=1.8;
possx=1:1:NA;
for k=1:NA
    A=A2;
    Bx=B;
    X=(A'*A)\(A'*Bx);
    [m,I]=min(abs(X));
    deleteordersigmax(k)=possx(I);
    A(:,I)=[];
    possx(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRsigmax(k)=norm(Bx)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Bx=B;
        Ai=A(IA(i,:)',:);
        Bi=Bx(IA(i,:)');
        A(IA(i,:)',:)=[];
        Bx(IA(i,:)')=[];
        SX=(A'*A)\(A'*Bx);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRsigmax(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRsigmax(k)/deltaSSRsigmax(k-1)>delta0)&&(deltaSSRsigmax(k-1)/deltaSSRsigmax(k-2)<delta0)
        break;
    end
end
Ix=deleteordersigmax(1:k-1);
A=A3;
A(:,Ix)=[];
Xsigma=(A'*A)\(A'*sigmax0');
exponx=0:Ncoef;
exponx(Ix)=[];
fxs=fsigma(zx0,exponx,Xsigma);
figure;
plot(zx0,sigmax0,'*');
hold on
plot(zx0,1-zx0+zx0.^2,'o');
hold on
plot(zx0,fxs,'r.');
hold off
% sigmay
Ncoef=4;
nz=length(zy0);
A=zeros(nz,Ncoef+1);
A(:,1)=1;
for i=1:Ncoef
    A(:,i+1)=zy0'.^i;
end
Ysigma0=(A'*A)\(A'*sigmay0');
A2=A;
B=sigmay0';
NA=length(A(1,:));
deleteordersigmay=zeros(1,NA);
deltaSSRsigmay=zeros(1,NA);
possy=1:1:NA;
for k=1:NA
    A=A2;
    By=B;
    X=(A'*A)\(A'*By);
    [m,I]=min(abs(X));
    deleteordersigmay(k)=possy(I);
    A(:,I)=[];
    possy(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRsigmay(k)=norm(By)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        By=B;
        Ai=A(IA(i,:)',:);
        Bi=By(IA(i,:)');
        A(IA(i,:)',:)=[];
        By(IA(i,:)')=[];
        SX=(A'*A)\(A'*By);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRsigmay(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRsigmay(k)/deltaSSRsigmay(k-1)>delta0)&&(deltaSSRsigmay(k-1)/deltaSSRsigmay(k-2)<delta0)
        break;
    end
end
Iy=deleteordersigmay(1:k-1);
A=A3;
A(:,Iy)=[];
Ysigma=(A'*A)\(A'*sigmay0');
expony=0:Ncoef;
expony(Iy)=[];
fys=fsigma(zy0,expony,Ysigma);
figure;
plot(zy0,sigmay0,'*');
hold on
plot(zy0,1./(1+0.5*zy0.^2),'o');
hold on
plot(zy0,fys,'r.');
hold off


% %%% Identify drift and diffusion by hard-thresholding sparse regression
% epsilong=1;
% Ncoef=3;
% I=(abs(zxf)<epsilong)&(abs(zyf)<epsilong);
% zxinitial=zxf0(I);
% zyinitial=zyf0(I);
% x=zxf(I);
% y=zyf(I);
% n=length(zxinitial);
% sx=fsigma(zxinitial,exponx,Xsigma);
% sy=fsigma(zyinitial,expony,Ysigma);
% % % neural network
% % path = sprintf('zx0.mat');
% % save(path,'zx0');
% % path = sprintf('sigmax0.mat');
% % save(path,'sigmax0');
% % path = sprintf('zxinitial.mat');
% % save(path,'zxinitial');
% % path = sprintf('zy0.mat');
% % save(path,'zy0');
% % path = sprintf('sigmay0.mat');
% % save(path,'sigmay0');
% % path = sprintf('zyinitial.mat');
% % save(path,'zyinitial');
% % testdata=linspace(-2,2,1000);
% % path = sprintf('testdata.mat');
% % save(path,'testdata');
% % load('sxtest.mat');
% % load('sytest.mat');
% % figure;
% % plot(testdata,sxtest);
% % figure;
% % plot(testdata,sytest);
% % sx=sx';
% % sy=sy';
% A=zeros(n,(Ncoef+1)*(Ncoef+2)/2);
% A(:,1)=1;
% for i=1:Ncoef
%     for j=1:i+1
%         A(:,i*(i+1)/2+j)=zxinitial'.^(i+1-j).*zyinitial'.^(j-1);
%     end
% end
% if alphax0==1
%     Bx=n/Nz*x'/h;
% else
%     Bx=n/Nz*x'/h-sx'.^alphax0*betax0*epsilong^(1-alphax0)*kalphax/(1-alphax0);
% end
% if alphay0==1
%     By=n/Nz*y'/h;
% else
%     By=n/Nz*y'/h-sy'.^alphay0*betay0*epsilong^(1-alphay0)*kalphay/(1-alphay0);
% end
% Bxx=n/Nz*x.^2'/h-sx'.^alphax0*epsilong^(2-alphax0)*kalphax/(2-alphax0);
% Byy=n/Nz*y.^2'/h-sy'.^alphay0*epsilong^(2-alphay0)*kalphay/(2-alphay0);
% Bxy=n/Nz*(x.*y)'/h;
% X0=(A'*A)\(A'*Bx);
% Y0=(A'*A)\(A'*By);
% XX0=(A'*A)\(A'*Bxx);
% YY0=(A'*A)\(A'*Byy);
% XY0=(A'*A)\(A'*Bxy);
% %%% Identify drift term
% Delta=0.05;
% A2=A;
% NA=length(A(1,:));
% posx=1:NA;
% for k=1:NA
%     X=(A'*A)\(A'*Bx);
%     I=abs(X)<Delta;
%     A(:,I)=[];
%     posx(I)=[];
%     if isempty(X(I))
%         break;
%     end
% end
% A=A2;
% posy=1:NA;
% for k=1:NA
%     Y=(A'*A)\(A'*By);
%     I=abs(Y)<Delta;
%     A(:,I)=[];
%     posy(I)=[];
%     if isempty(Y(I))
%         break;
%     end
% end
% %%% Identify diffusion term
% A=A2;
% posxx=1:NA;
% for k=1:NA
%     XX=(A'*A)\(A'*Bxx);
%     I=abs(XX)<Delta;
%     A(:,I)=[];
%     posxx(I)=[];
%     if isempty(XX(I))
%         break;
%     end
% end
% A=A2;
% posyy=1:NA;
% for k=1:NA
%     YY=(A'*A)\(A'*Byy);
%     I=abs(YY)<Delta;
%     A(:,I)=[];
%     posyy(I)=[];
%     if isempty(YY(I))
%         break;
%     end
% end
% A=A2;
% posxy=1:NA;
% for k=1:NA
%     XY=(A'*A)\(A'*Bxy);
%     I=abs(XY)<Delta;
%     A(:,I)=[];
%     posxy(I)=[];
%     if isempty(XY(I))
%         break;
%     end
% end


%%% Identify drift and diffusion terms using SSR sparsing
epsilong=1;
Ncoef=3;
I=(abs(zxf)<epsilong)&(abs(zyf)<epsilong);
zxinitial=zxf0(I);
zyinitial=zyf0(I);
x=zxf(I);
y=zyf(I);
n=length(zxinitial);
sx=fsigma(zxinitial,exponx,Xsigma);
sy=fsigma(zyinitial,expony,Ysigma);
if alphax0==1
    Bx=n/Nz*x'/h;
else
    Bx=n/Nz*x'/h-sx'.^alphax0*betax0*epsilong^(1-alphax0)*kalphax/(1-alphax0);
end
if alphay0==1
    By=n/Nz*y'/h;
else
    By=n/Nz*y'/h-sy'.^alphay0*betay0*epsilong^(1-alphay0)*kalphay/(1-alphay0);
end
Bxx=n/Nz*x.^2'/h-sx'.^alphax0*epsilong^(2-alphax0)*kalphax/(2-alphax0);
Byy=n/Nz*y.^2'/h-sy'.^alphay0*epsilong^(2-alphay0)*kalphay/(2-alphay0);
Bxy=n/Nz*(x.*y)'/h;
nbinsq=30;
nbin=nbinsq^2;
xbin=zeros(1,nbin);
ybin=zeros(1,nbin);
Bxbin=zeros(1,nbin);
Bybin=zeros(1,nbin);
Bbinxx=zeros(1,nbin);
Bbinyy=zeros(1,nbin);
Bbinxy=zeros(1,nbin);
wbin=zeros(1,nbin);
xl=linspace(zxmin,zxmax,nbinsq+1);
yl=linspace(zymin,zymax,nbinsq+1);
for j=1:nbinsq
    for i=1:nbinsq
        Ix=(zxinitial>xl(i))&(zxinitial<xl(i+1));
        Iy=(zyinitial>yl(j))&(zyinitial<yl(j+1));
        I=Ix&Iy;
        n0=length(zxinitial(I));
        wbin((j-1)*nbinsq+i)=n0;
        xbin((j-1)*nbinsq+i)=sum(zxinitial(I))/n0;
        ybin((j-1)*nbinsq+i)=sum(zyinitial(I))/n0;
        Bxbin((j-1)*nbinsq+i)=sum(Bx(I))/n0;
        Bybin((j-1)*nbinsq+i)=sum(By(I))/n0;
        Bbinxx((j-1)*nbinsq+i)=sum(Bxx(I))/n0;
        Bbinyy((j-1)*nbinsq+i)=sum(Byy(I))/n0;
        Bbinxy((j-1)*nbinsq+i)=sum(Bxy(I))/n0;
    end
end
wbin=diag(wbin)/(n/nbin);
A=zeros(nbin,(Ncoef+1)*(Ncoef+2)/2);
A(:,1)=1;
for i=1:Ncoef
    for j=1:i+1
        A(:,i*(i+1)/2+j)=xbin'.^(i+1-j).*ybin'.^(j-1);
    end
end
A=wbin*A;
Bx=wbin*Bxbin';
By=wbin*Bybin';
Bxx=wbin*Bbinxx';
Byy=wbin*Bbinyy';
Bxy=wbin*Bbinxy';

X0=(A'*A)\(A'*Bx);
Y0=(A'*A)\(A'*By);
XX0=(A'*A)\(A'*Bxx);
YY0=(A'*A)\(A'*Byy);
XY0=(A'*A)\(A'*Bxy);

kfold=2;
nblock=floor(nbin/kfold);
n=kfold*nblock;
IA=zeros(kfold,nblock);
position=1:n;
for i=1:kfold-1
    u=rand(1,nblock);
    J=kfold+1-i;
    p=0:J:J*(nblock-1);
    IA(i,:)=position(p+floor(J*u)+1);
    position(p+floor(J*u)+1)=[];
end
IA(kfold,:)=position;

% drift b1(x)
A3=A;
A2=A;
Bx2=Bx;
NA=length(A(1,:));
deleteorderx=zeros(1,NA);
deltaSSRx=zeros(1,NA);
posx=1:1:NA;
for k=1:NA
    A=A2;
    Bx=Bx2;
    X=(A'*A)\(A'*Bx);
    [m,I]=min(abs(X));
    deleteorderx(k)=posx(I);
    A(:,I)=[];
    posx(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRx(k)=norm(Bx)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Bx=Bx2;
        Ai=A(IA(i,:)',:);
        Bi=Bx(IA(i,:)');
        A(IA(i,:)',:)=[];
        Bx(IA(i,:)')=[];
        SX=(A'*A)\(A'*Bx);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRx(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRx(k)/deltaSSRx(k-1)>delta0)&&(deltaSSRx(k-1)/deltaSSRx(k-2)<delta0)
        break;
    end
end

% drift b2(x)
A2=A3;
By2=By;
deleteordery=zeros(1,NA);
deltaSSRy=zeros(1,NA);
posy=1:1:NA;
for k=1:NA
    A=A2;
    By=By2;
    Y=(A'*A)\(A'*By);
    [m,I]=min(abs(Y));
    deleteordery(k)=posy(I);
    A(:,I)=[];
    posy(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRy(k)=norm(By)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        By=By2;
        Ai=A(IA(i,:)',:);
        Bi=By(IA(i,:)');
        A(IA(i,:)',:)=[];
        By(IA(i,:)')=[];
        SX=(A'*A)\(A'*By);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRy(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRy(k)/deltaSSRy(k-1)>delta0)&&(deltaSSRy(k-1)/deltaSSRy(k-2)<delta0)
        break;
    end
end

% diffusion sigma11(x)
A2=A3;
Bxx2=Bxx;
deleteorderxx=zeros(1,NA);
deltaSSRxx=zeros(1,NA);
posxx=1:1:NA;
for k=1:NA
    A=A2;
    Bxx=Bxx2;
    XX=(A'*A)\(A'*Bxx);
    [m,I]=min(abs(XX));
    deleteorderxx(k)=posxx(I);
    A(:,I)=[];
    posxx(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRxx(k)=norm(Bxx)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Bxx=Bxx2;
        Ai=A(IA(i,:)',:);
        Bi=Bxx(IA(i,:)');
        A(IA(i,:)',:)=[];
        Bxx(IA(i,:)')=[];
        SX=(A'*A)\(A'*Bxx);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRxx(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRxx(k)/deltaSSRxx(k-1)>delta0)&&(deltaSSRxx(k-1)/deltaSSRxx(k-2)<delta0)
        break;
    end
end

% diffusion sigma22(x)
A2=A3;
Byy2=Byy;
deleteorderyy=zeros(1,NA);
deltaSSRyy=zeros(1,NA);
posyy=1:1:NA;
for k=1:NA
    A=A2;
    Byy=Byy2;
    YY=(A'*A)\(A'*Byy);
    [m,I]=min(abs(YY));
    deleteorderyy(k)=posyy(I);
    A(:,I)=[];
    posyy(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRyy(k)=norm(Byy)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Byy=Byy2;
        Ai=A(IA(i,:)',:);
        Bi=Byy(IA(i,:)');
        A(IA(i,:)',:)=[];
        Byy(IA(i,:)')=[];
        SX=(A'*A)\(A'*Byy);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRyy(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRyy(k)/deltaSSRyy(k-1)>delta0)&&(deltaSSRyy(k-1)/deltaSSRyy(k-2)<delta0)
        break;
    end
end

% diffusion sigma12(x)
A2=A3;
Bxy2=Bxy;
deleteorderxy=zeros(1,NA);
deltaSSRxy=zeros(1,NA);
posxy=1:1:NA;
for k=1:NA
    A=A2;
    Bxy=Bxy2;
    XY=(A'*A)\(A'*Bxy);
    [m,I]=min(abs(XY));
    deleteorderxy(k)=posxy(I);
    A(:,I)=[];
    posxy(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRxy(k)=norm(Bxy)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Bxy=Bxy2;
        Ai=A(IA(i,:)',:);
        Bi=Bxy(IA(i,:)');
        A(IA(i,:)',:)=[];
        Bxy(IA(i,:)')=[];
        SX=(A'*A)\(A'*Bxy);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRxy(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRxy(k)/deltaSSRxy(k-1)>delta0)&&(deltaSSRxy(k-1)/deltaSSRxy(k-2)<delta0)
        break;
    end
end

Nlin=1000;
xlin=linspace(zxmin,zxmax,Nlin);
ylin=linspace(zymin,zymax,Nlin);
[XL,YL]=meshgrid(xlin,ylin);
b1xtrue=XL-XL.^3-5*XL.*YL.^2;
b1xlearn=1.0478*XL-1.0086*XL.^3-4.9929*XL.*YL.^2;
figure;
mesh(XL,YL,b1xtrue);
figure;
mesh(XL,YL,b1xlearn);
b2xtrue=-(1+XL.^2).*YL;
b2xlearn=-1.0031*YL-1.0043*XL.^2.*YL;
figure;
mesh(XL,YL,b2xtrue);
figure;
mesh(XL,YL,b2xlearn);
a11xtrue=2+2*XL+XL.^2;
a11xlearn=1.9682+1.9892*XL+1.1169*XL.^2;
figure;
mesh(XL,YL,a11xtrue);
figure;
mesh(XL,YL,a11xlearn);
a12xtrue=YL;
a12xlearn=0.9954*YL;
figure;
mesh(XL,YL,a12xtrue);
figure;
mesh(XL,YL,a12xlearn);
a22xtrue=YL.^2;
a22xlearn=0.9894*YL.^2;
figure;
mesh(XL,YL,a22xtrue);
figure;
mesh(XL,YL,a22xlearn);
