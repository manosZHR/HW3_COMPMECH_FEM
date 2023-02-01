clear;close all;clc

nex = 8;
ney = 8;
lambda = 0;
alpha = 1;
omega = 1;   

nnx = 2*nex+1;
nny = 2*ney+1;
np=nnx*nny;

uinit=zeros(np,1);
tic;
[xpt ypt unew res_ev] = hw3_final(nex,ney,lambda,alpha,omega,uinit);
toc;

n2=nny;
l2=1;
k2=1;
for j = 1:nnx
    x(1:nny,k2)=xpt(l2:n2,1);
    y(1:nny,k2)=ypt(l2:n2,1);
    n2=n2+nny;
    k2=k2+1;
    l2=l2+nny;
    
 end
z=zeros(nnx,nny);

figure(1)
mesh(x,y,z,'edgecolor','k')
title('2D-Mesh')
xlabel('x')
ylabel('y')


[xi, yi] = meshgrid(linspace(min(xpt),max(xpt),length(xpt)),linspace(min(ypt),max(ypt),length(ypt)));
zi = griddata(xpt,ypt,unew,xi,yi);

figure(2)
contour(xi,yi,zi)
h=colorbar;
colormap jet
ylabel(h,'concentration','FontSize',14)
xlabel('x')
ylabel('y')

figure(3)
surf(xi,yi,zi,'edgecolor','none')
h=colorbar;
colormap jet
ylabel(h,'concentration','FontSize',14)
xlabel('x')
ylabel('y')

figure(4)
semilogy(1:length(res_ev),res_ev,'k')
xlabel('# iterations')
ylabel('Residual')

