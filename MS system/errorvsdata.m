% % errors for sparse regression and neural network
% testdata=linspace(-2,2,1000);
% sigmaxtrue=1-testdata+testdata.^2;
% sigmaytrue=1./(1+0.5*testdata.^2);
% srx=fsigma(testdata,exponx,Xsigma);
% sry=fsigma(testdata,expony,Ysigma);
% errorxnn=sqrt(mean((sxtest'-sigmaxtrue).^2))/max(abs(sigmaxtrue));
% errorynn=sqrt(mean((sytest'-sigmaytrue).^2))/max(abs(sigmaytrue));
% errorxsr=sqrt(mean((srx-sigmaxtrue).^2))/max(abs(sigmaxtrue));
% errorysr=sqrt(mean((sry-sigmaytrue).^2))/max(abs(sigmaytrue));

N=1000;
x=linspace(4.2,9,N);
M=10.^x;
ealpha=M.^(-0.5);
ealpha=ealpha*20;

% error vs data volume
Ndata=[3e4,1e5,3e5,1e6,3e6,1e7,3e7,1e8,3e8,1e9];
erroralpha1=[0.2428,0.1263,0.0602,0.0326,0.0336,0.0138,0.0024,0.0052,0.0079,0.0056];
erroralpha2=[0.6238,0.2424,0.0342,0.0357,0.0290,0.0182,0.0222,0.0022,0.0095,0.0059];
errorbeta1=[0.1490,0.0085,0.0610,0.0089,0.0046,0.0086,0.0048,0.0042,0.0027,0.0020];
errorbeta2=[0.0984,0.0593,0.0839,0.0292,0.0024,0.0215,0.0063,0.0027,0.0032,0.0005];
errorsigma1=[26.8932,1.0694,0.2685,0.0655,0.0535,0.0292,0.0166,0.0072,0.0111,0.0045];
errorsigma2=[0.5744,0.2884,0.1524,0.0810,0.0362,0.0307,0.0226,0.0085,0.0108,0.0091];
figure;
plot(Ndata,erroralpha1);
hold on
plot(Ndata,erroralpha2);
hold on
plot(M,ealpha);
hold off

figure;
plot(Ndata,errorbeta1);
hold on
plot(Ndata,errorbeta2);
hold on
plot(M,ealpha);
hold off

figure;
plot(Ndata,errorsigma1);
hold on
plot(Ndata,errorsigma2);
hold on
plot(M,ealpha);
hold off

% error vs data volume with Mh fixed
Ndata2=[1e6,3e6,1e7,3e7,1e8,3e8,1e9];
errora1=[0.9554,0.1713,0.0138,0.0097,0.0105,0.0061,0.0091];
errora2=[0.2278,0.0495,0.0182,0.0116,0.0321,0.0185,0.0242];
errorb1=[0.7245,0.1923,0.0086,0.0012,0.0026,0.0090,0.0034];
errorb2=[0.1218,0.0032,0.0215,0.0132,0.0045,0.0071,0.0127];
errors1=[0.3282,0.2540,0.0292,0.0403,0.0443,0.0240,0.0304];
errors2=[0.4874,0.0557,0.0307,0.0221,0.0230,0.0189,0.0278];
figure;
plot(Ndata2,errora1);
hold on
plot(Ndata2,errora2);
hold off

figure;
plot(Ndata2,errorb1);
hold on
plot(Ndata2,errorb2);
hold off

figure;
plot(Ndata2,errors1);
hold on
plot(Ndata2,errors2);
hold off

% error vs h^-1 with M fixed
hinv=[1/0.01,1/0.003,1/0.001,1/0.0003,1/0.0001,1/0.00003];
ea1=[0.9638,0.1179,0.0138,0.0293,0.0284,0.0681];
ea2=[0.2190,0.0060,0.0182,0.1010,0.0607,0.3760];
eb1=[0.7222,0.1307,0.0086,0.0098,0.0264,0.0523];
eb2=[0.1249,0.0050,0.0215,0.0496,0.0625,0.0728];
es1=[0.3265,0.2009,0.0292,0.0784,0.1888,0.1415];
es2=[0.4706,0.0248,0.0307,0.0480,0.0827,0.1916];
figure;
plot(hinv,ea1);
hold on
plot(hinv,ea2);
hold off

figure;
plot(hinv,eb1);
hold on
plot(hinv,eb2);
hold off

figure;
plot(hinv,es1);
hold on
plot(hinv,es2);
hold off

N=1000;
x=linspace(0.1,10,N);
fx=1./x+sqrt(x);
figure;
plot(x,fx);
