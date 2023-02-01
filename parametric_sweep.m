clear;close all, clc
nex = 8;
ney = 8;

nnx = 2*nex+1;
nny = 2*ney+1;
np=nnx*nny;

alpha = 1;
omega = 1; 
lambda_origin = 0;
lambda_last = 10;
dl=0.5;

uinit=zeros(np,1);
uold=uinit;

maxiter = (lambda_last-lambda_origin)/dl;

normu=zeros(maxiter,1);

lambda=zeros(maxiter,1);
lambda(1)=lambda_origin;

tic;
iter=1;
while true %parametric sweep loop
   
   fprintf('Iteration: %d',iter)
   fprintf('/%d',maxiter)
   fprintf('\n')

   [xpt ypt unew res_ev] = hw3_final(nex,ney,lambda(iter),alpha,omega,uold);
   normu(iter) = norm(unew,2);
   uold=unew;

   tot_iter(iter)=length(res_ev); %number of iterations to converge per outside loop (parametric sweep loop)

   lambda(iter+1) = lambda(iter)+dl;
   iter=iter+1;

   if iter==maxiter+1
       break;end
end
toc;

figure(1)
plot(lambda(1:end-1),normu')
xlabel('lambda'),ylabel('||u||'),title('Parametric sweep (lambda)')

figure(2)
plot(tot_iter)
xlabel('Parametric sweep iteration'),ylabel('# iterations until convergense'),title('Parametric sweep (lambda)')


[xi, yi] = meshgrid(linspace(min(xpt),max(xpt),length(xpt)),linspace(min(ypt),max(ypt),length(ypt)));
zi = griddata(xpt,ypt,unew,xi,yi);

figure(3)
contour(xi,yi,zi)
h=colorbar;
colormap jet
ylabel(h,'concentration','FontSize',14)
xlabel('x')
ylabel('y')
title('Contour plot for λ= ',lambda(iter),'Fontweight','bold','Fontsize',12)



