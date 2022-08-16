close all;
clear;
clc;

rng(2,'philox');
theta1=randi([0,179],18,1);
I1=zeros(217);
I1(19:199,:)=double(imread('slice_50.png'));
% figure;
% imshow(I1,[]);

R1=radon(I1,theta1);
I1a=iradon(R1,theta1);
figure;
imshow(I1a,[]);

yc1=reshape(R1,[],1);
m=309;
mc1=309*18;
nc1=217*217;
A=CSHelperCode(m,nc1,theta1);
At=A';

[xc1,status]=l1_ls(A,At,mc1,nc1,yc1,0.05);
xc1=reshape(xc1,217,217);
xc1=idct2(xc1);
figure;
imshow(xc1,[]);

theta2=randi([0,179],18,1);
I2=zeros(217,217);
I2(19:199,:)=double(imread('slice_51.png'));
% figure;
% imshow(I2,[]);

R2=radon(I2,theta2);
yc2=reshape(R2,[],1);
yc2=[yc1;yc2];
mc2=2*309*18;
nc2=2*217*217;
Ac2=CoupledCSHelperCode(m,nc1,2,{theta1,theta2});
Atc2=Ac2';

[xc2,statusc2]=l1_ls(Ac2,Atc2,mc2,nc2,yc2,0.5);
xc21=reshape(xc2(1:47089),217,217);
xc22=reshape(xc2(47090:94178),217,217);
xc22=xc21+xc22;
xc21=idct2(xc21);
xc22=idct2(xc22);
figure;
imshow(xc21,[]);
figure;
imshow(xc22,[]);

theta3=randi([0,179],18,1);
I3=zeros(217,217);
I3(19:199,:)=double(imread('slice_52.png'));
% figure;
% imshow(I3,[]);

R3a=radon(I3,theta3);
yc3=reshape(R3a,[],1);
yc3=[yc2;yc3];
mc3=3*309*18;
nc3=3*217*217;
Ac3=CoupledCSHelperCode(m,nc1,3,{theta1,theta2,theta3});
Atc3=Ac3';

[xc3,statusc3]=l1_ls(Ac3,Atc3,mc3,nc3,yc3,5);
xc31=reshape(xc3(1:47089),217,217);
xc32=reshape(xc3(47090:94178),217,217);
xc33=reshape(xc3(94179:141267),217,217);
xc32=xc31+xc32;
xc33=xc32+xc33;
xc31=idct2(xc31);
xc32=idct2(xc32);
xc33=idct2(xc33);
figure;
imshow(xc31,[]);
figure;
imshow(xc32,[]);
figure;
imshow(xc33,[]);
