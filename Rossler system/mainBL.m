clear;
clc;

zm=2;
zxmin=-zm;
zxmax=zm;
zymin=-zm;
zymax=zm;
zzmin=-zm;
zzmax=zm;
Nz0=2e4;
node=5;

h=0.001;
NT=1e3;
T=NT*h;
Nz=NT*Nz0;
alphax=0.5;
alphay=1.01;
alphaz=1.5;
betax=0.5;
betay=0;
betaz=-0.5;

%%% initial points
zxf0s=rand(node,Nz0)*2*zm-zm;
zyf0s=rand(node,Nz0)*2*zm-zm;
zzf0s=rand(node,Nz0)*2*zm-zm;
zxf1s=zeros(node,Nz0);
zyf1s=zeros(node,Nz0);
zzf1s=zeros(node,Nz0);
zxf0=[];
zyf0=[];
zzf0=[];
zxf=[];
zyf=[];
zzf=[];
am=2;

%%% Generate data
Asys=[0,0.7161,0,0,0;
    0,0,-1.2678,0,0;
    0,0,0,-1.0141,0;
    0,0,0,0,2.3633;
    2.7307,0,0,0,0];
for k=1:NT
    zxf0=[zxf0,zxf0s];
    zyf0=[zyf0,zyf0s];
    zzf0=[zzf0,zzf0s];
    Mx=stblrnd(alphax,betax,h^(1/alphax),0,node,Nz0);
    My=stblrnd(alphay,betay,h^(1/alphay),0,node,Nz0);
    Mz=stblrnd(alphaz,betaz,h^(1/alphaz),0,node,Nz0);
    Bh=sqrt(h)*randn(3*node,Nz0);
    for i=1:node
        Gi=0;
        for j=1:node
            Gi=Gi+Asys(i,j)*(zxf0s(j,:)-zxf0s(i,:));
        end
        zxf1s(i,:)=zxf0s(i,:)+(-zyf0s(i,:)-zzf0s(i,:)+Gi)*h+zxf0s(i,:).*Bh(i,:)+(1+0.5*zxf0s(i,:).^2).*Mx(i,:);
        zyf1s(i,:)=zyf0s(i,:)+(zxf0s(i,:)+0.35*zyf0s(i,:))*h+sqrt(1+zyf0s(i,:).^2).*Bh(node+i,:)+1./(1+0.5*zyf0s(i,:).^2).*My(i,:);
        zzf1s(i,:)=zzf0s(i,:)+(0.2+zyf0s(i,:).*(zxf0s(i,:)-5.7))*h+sqrt(1+0.5*zyf0s(i,:).^2+(zzf0s(i,:)).^2).*Bh(2*node+i,:)+(1+sin(zzf0s(i,:)).^2).*Mz(i,:);
    end
    zxf=[zxf,zxf1s-zxf0s];
    zyf=[zyf,zyf1s-zyf0s];
    zzf=[zzf,zzf1s-zzf0s];
    for i=1:node
        I=abs(zxf1s(i,:))>am;
        zxf1s(i,I)=rand(1,length(zxf1s(i,I)))*2*zm-zm;
        I=abs(zyf1s(i,:))>am;
        zyf1s(i,I)=rand(1,length(zyf1s(i,I)))*2*zm-zm;
        I=abs(zzf1s(i,:))>am;
        zzf1s(i,I)=rand(1,length(zzf1s(i,I)))*2*zm-zm;
    end
    zxf0s=zxf1s;
    zyf0s=zyf1s;
    zzf0s=zzf1s;
end
zs0=[zxf0;zyf0;zzf0];
zs=[zxf;zyf;zzf];
clear zxf0 zyf0 zzf0 zxf zyf zzf

%%% Identify alpha beta sigma
alphatotal=[];
betatotal=[];
ztotal=[];
sigmatotal=[];
for is=1:3*node
    zsi0=zs0(is,:);
    zsi=zs(is,:);
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
            I=(zsi0>=z0(i)-dz)&(zsi0<z0(i)+dz);
            zxfi0=zsi0(I);
            zxfi=zsi(I);
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
    alphatotal=[alphatotal,alphax0];
    betatotal=[betatotal,betax0];
    ztotal=[ztotal,zx0'];
    sigmatotal=[sigmatotal,sigmax0];
end
clear zsi zsi0

%%% Identify sigmax and sigmay using SSR sparsing
Ncoef=4;
nz=length(ztotal(:,1));
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
delta0=1.5;
NA=Ncoef+1;
Xsigma0=zeros(3*node,NA);
Xsigma=zeros(3*node,NA);
deleteordersigma=zeros(3*node,NA);
deltaSSRsigma=zeros(3*node,NA);
possigma=zeros(3*node,NA);
expsigma=zeros(3*node,NA);
LXsigma=zeros(3*node,1);
for is=1:3*node
    A=zeros(nz,NA);
    A(:,1)=1;
    for i=1:Ncoef
        A(:,i+1)=ztotal(:,is).^i;
    end
    Xs0=(A'*A)\(A'*sigmatotal(:,is));
    % SSR
    A3=A;
    A2=A;
    B=sigmatotal(:,is);
    deleteordersigmax=zeros(1,NA);
    deltaSSRsigmax=zeros(1,NA);
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
    Xs=(A'*A)\(A'*sigmatotal(:,is));
    exponx=0:Ncoef;
    exponx(Ix)=[];
    Xsigma0(is,:)=Xs0';
    Xsigma(is,1:length(Xs))=Xs';
    deleteordersigma(is,:)=deleteordersigmax;
    deltaSSRsigma(is,:)=deltaSSRsigmax;
    possigma(is,1:length(possx))=possx;
    expsigma(is,1:length(exponx))=exponx;
    LXsigma(is)=length(Xs);
end
n=14;
fxs=fsigma(ztotal(:,n),expsigma(n,1:LXsigma(n)),Xsigma(n,1:LXsigma(n)));
figure;
plot(ztotal(:,n),sigmatotal(:,n),'*');
hold on
% plot(ztotal(:,n),1+0.5*ztotal(:,n).^2,'o');
% plot(ztotal(:,n),1./(1+0.5*ztotal(:,n).^2),'o');
plot(ztotal(:,n),1+sin(ztotal(:,n)).^2,'o');
hold on
plot(ztotal(:,n),fxs,'r.');
hold off
% % neural network
% zx0=ztotal(:,n);
% sigmax0=sigmatotal(:,n);
% path = sprintf('zx0.mat');
% save(path,'zx0');
% path = sprintf('sigmax0.mat');
% save(path,'sigmax0');
% testdata=linspace(-2,2,1000);
% path = sprintf('testdata.mat');
% save(path,'testdata');
load('sxtest15.mat');
figure;
plot(testdata,sxtest);


%%% Identify drift and diffusion by hard-thresholding sparse regression
epsilong=1;
Ncoef=2;
I1=ones(1,Nz,'logical');
for i=1:3*node
    I1=I1&(abs(zs0(i,:))<zm);
end
zs01=zs0(:,I1);
zs1=zs(:,I1);
Nz1=length(zs1(1,:));
I=ones(1,Nz1,'logical');
for i=1:3*node
    I=I&(abs(zs1(i,:))<epsilong);
end
zxinitial=zs01(:,I);
x=zs1(:,I);
n=length(zxinitial(1,:));
sx=[];
for i=1:3*node
    sx=[sx,fsigma(zxinitial(i,:)',expsigma(i,1:LXsigma(i)),Xsigma(i,1:LXsigma(i)))];
end
A=zeros(n,1+3*node+3*node*(3*node+1)/2);
A(:,1)=1;
for i=1:3*node
    A(:,i+1)=zxinitial(i,:)';
end
for i=1:3*node
    A(:,i+1+3*node)=(zxinitial(i,:)').^2;
end
pairs = nchoosek(1:3*node, 2);
for i=1:3*node*(3*node-1)/2
    A(:,i+1+6*node)=(zxinitial(pairs(i,1),:)').*(zxinitial(pairs(i,2),:)');
end
Delta=0.1;
A2=A;
NA=1+3*node+3*node*(3*node+1)/2;
X=zeros(NA,3*node+3*node*(3*node+1)/2);
posdd=zeros(NA,3*node+3*node*(3*node+1)/2);
LX=zeros(1,3*node+3*node*(3*node+1)/2);
for i=1:3*node
    kalphax=alphatotal(i)*(1-alphatotal(i))/(cos(pi*alphatotal(i)/2)*gamma(2-alphatotal(i)));
    if alphatotal(i)==1
        B=n/Nz1*x(i,:)'/h;
    else
        B=n/Nz1*x(i,:)'/h-sx(:,i).^alphatotal(i)*betatotal(i)*epsilong^(1-alphatotal(i))*kalphax/(1-alphatotal(i));
    end

    A=A2;
    posx=1:NA;
    for k=1:NA
        X1=(A'*A)\(A'*B);
        I=abs(X1)<Delta;
        A(:,I)=[];
        posx(I)=[];
        if isempty(X1(I))
            break;
        end
    end
    LX(i)=length(posx);
    X(1:LX(i),i)=X1;
    posdd(1:LX(i),i)=posx';
end
for i=1:3*node
    kalphax=alphatotal(i)*(1-alphatotal(i))/(cos(pi*alphatotal(i)/2)*gamma(2-alphatotal(i)));
    B=n/Nz1*x(i,:).^2'/h-sx(:,i).^alphatotal(i)*epsilong^(2-alphatotal(i))*kalphax/(2-alphatotal(i));

    A=A2;
    posx=1:NA;
    for k=1:NA
        X1=(A'*A)\(A'*B);
        I=abs(X1)<Delta;
        A(:,I)=[];
        posx(I)=[];
        if isempty(X1(I))
            break;
        end
    end
    LX(i+3*node)=length(posx);
    X(1:LX(i+3*node),i+3*node)=X1;
    posdd(1:LX(i+3*node),i+3*node)=posx';
end
for i=1:3*node*(3*node-1)/2
    B=n/Nz1*(x(pairs(i,1),:).*x(pairs(i,2),:))'/h;

    A=A2;
    posx=1:NA;
    for k=1:NA
        X1=(A'*A)\(A'*B);
        I=abs(X1)<Delta;
        A(:,I)=[];
        posx(I)=[];
        if isempty(X1(I))
            break;
        end
    end
    LX(i+6*node)=length(posx);
    X(1:LX(i+6*node),i+6*node)=X1;
    posdd(1:LX(i+6*node),i+6*node)=posx';
end
