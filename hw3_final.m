function [xpt,ypt,unew,res_ev]=hw3_final(nex,ney,lambda,alpha,omega,u0)
% 2-dimensional linear problem. 
% Quadratic finite element basis functions.
%
% MATLAB version of the original FORTRAN code 1d_FEM_quad.f90
%
% Syntax
%  [x y u]=FEM_2d_quad(nex,ney)
%
% Parent function: FEM_2d_quad
% Nested functions: xydiscr, nodnumb, xycoord, abfind, tsfun
%
% This code uses nested functions.
% The primary difference between nested functions and other types
% of functions is that they can access and modify variables
% that are defined in their parent functions.
% http://www.mathworks.com/help/matlab/matlab_prog/nested-functions.html

%Global variables
nnx = 2*nex+1;
nny = 2*ney+1;
ne  = nex*ney;
np  = nnx*nny;
r1 = []; J=[];
xorigin = []; yorigin = []; xlast = []; ylast = [];
deltax = []; deltay = []; xpt = []; ypt = [];
nop = []; ntop = []; nlat = [];
ncod = []; bc = [];
w = []; gp =[];
phi = []; phic = []; phie = [];
unew=[]; uold=[]; du=[]; uinit=[]; iter=[]; maxiter=[]; tol_1=[]; tol_2=[]; r=[]; res_ev=[]; err=[];
res=[]; TF=[];

%fprintf('2-D problem. Biquadratic basis functions \n')

xydiscr();
nodnumb();
xycoord();

%fprintf('nex=%d, ney=%d, ne=%d, np=%d\n',nex,ney,ne,np)


% prepare for essential boundary conditions
ncod = zeros(np,1);
bc   = zeros(np,1);

%Left Boundary
ncod(1:nny) = 1;
bc(1:nny)   = 0;

%Right Boundary
ncod(np-nny+1:np) = 1;
bc(np-nny+1:np)   = 0;

%Top Boundary
ncod(nny:nny:np) = 1;
bc(1:nny)   = 0;

%Bottom Boundary
ncod(1:nny:np-nny+1) = 1;
bc(1:nny:np-nny+1)   = 0;

% initialization
r1 = zeros(np,1);
J = zeros(np,np);


% solve the system of equations by Gauss elimination
du=zeros(np,1);
uinit=u0;%zeros(np,1);
uold=uinit;

     
r = zeros(np,1);
res_ev=[];

tol_1=1e-5; %tolerance for residual
tol_2=1e-8; %tolerance for unew-uold

maxiter=3000;

iter=0;
    while true 
        
        iter=iter+1;
        % matrix assembly
        for nell=1:ne
            abfind(nell,uold)
        end

        % impose essential boundary conditions
        for i=1:np
            if(ncod(i)==1)
                r1(i)= -( uold(i)-bc(i) );%bc(i);
                J(i,1:np)=0.;
                J(i,i)=1.;
            end
        end
        
        du = sparse(J)\r1;
        unew = uold + omega*du;
        
        TF = isnan(unew); %Flag for NaN values
        
        err  = abs(norm(unew-uold));
        res = err / (norm(unew,2));
        res_ev(iter)=res; %norm(du);
        
        uold = unew;

        if norm(du,2) <tol_1 %norm(r1,2) <tol_1
            fprintf('Iteration stopped with tol_1 \n')
            fprintf('Total number of iterations: %d \n',iter)
            break;end
        if res < tol_2
            fprintf('Iteration stopped with tol_2 \n')
            fprintf('Total number of iterations: %d \n',iter)
            break;end
        if iter>=maxiter
            fprintf('Maximum iterations reached \n')
            fprintf('Total number of iterations: %d \n',iter)
            break;end
        if ismember(1,TF)==1
            break;end

    
    end
   
    res_drop=-log10(res_ev(end)/res_ev(1));
    fprintf('Residual drop factor: %.2f \n',res_drop) % how many orders the residual drops from the first iteration 
    
    function xydiscr()
        xorigin = 0.;
        yorigin = 0.;
        xlast   = 1.;
        ylast   = 1.;
        deltax = (xlast-xorigin)/nex;
        deltay = (ylast-yorigin)/ney;
    end

    function nodnumb()
    % *** nodal numbering
        nel=0;
        for i=1:nex
            for j=1:ney
                nel=nel+1;
                for k=1:3
                    l=3*k-2;
                    nop(nel,l)=nny*(2*i+k-3)+2*j-1;
                    nop(nel,l+1)=nop(nel,l)+1;
                    nop(nel,l+2)=nop(nel,l)+2;
                end
            end
        end
    end

    function xycoord()
    % *** (x,y) coordinates of each node
        xpt=zeros(np,1);
        ypt=zeros(np,1);
        xpt(1)=xorigin;
        ypt(1)=yorigin;
        for i=1:nnx
            nnode=(i-1)*nny+1;
            xpt(nnode)=xpt(1)+(i-1)*deltax/2.;
            ypt(nnode)=ypt(1);
            for j=2:nny
                xpt(nnode+j-1)=xpt(nnode);
                ypt(nnode+j-1)=ypt(nnode)+(j-1)*deltay/2.;
            end
        end
    end

    
    function tsfun(x,y)
        l1  =@(c)  2.*c^2-3.*c+1.;
        l2  =@(c) -4.*c^2+4.*c;
        l3  =@(c)  2.*c^2-c;
        dl1 =@(c)  4.*c-3.;
        dl2 =@(c) -8.*c+4.;
        dl3 =@(c)  4.*c-1.;
        
        phi(1)=l1(x)*l1(y);
        phi(2)=l1(x)*l2(y);
        phi(3)=l1(x)*l3(y);
        phi(4)=l2(x)*l1(y);
        phi(5)=l2(x)*l2(y);
        phi(6)=l2(x)*l3(y);
        phi(7)=l3(x)*l1(y);
        phi(8)=l3(x)*l2(y);
        phi(9)=l3(x)*l3(y);
        phic(1)=l1(y)*dl1(x);
        phic(2)=l2(y)*dl1(x);
        phic(3)=l3(y)*dl1(x);
        phic(4)=l1(y)*dl2(x);
        phic(5)=l2(y)*dl2(x);
        phic(6)=l3(y)*dl2(x);
        phic(7)=l1(y)*dl3(x);
        phic(8)=l2(y)*dl3(x);
        phic(9)=l3(y)*dl3(x);
        phie(1)=l1(x)*dl1(y);
        phie(2)=l1(x)*dl2(y);
        phie(3)=l1(x)*dl3(y);
        phie(4)=l2(x)*dl1(y);
        phie(5)=l2(x)*dl2(y);
        phie(6)=l2(x)*dl3(y);
        phie(7)=l3(x)*dl1(y);
        phie(8)=l3(x)*dl2(y);
        phie(9)=l3(x)*dl3(y);
    end
    
    function abfind(nell,u)
        w  = [0.27777777777778, 0.444444444444, 0.27777777777778];
        gp = [0.1127016654    , 0.5           , 0.8872983346    ];
        
        ngl = nop(nell,1:9);
        
        for j=1:3 %LOOP j
            for k=1:3 %LOOPk
                tsfun(gp(j),gp(k))
                % *** isoparametric transformation
                x1=0;x2=0;y1=0;y2=0;
                for n=1:9
                    x1=x1+xpt(ngl(n))*phic(n);
                    x2=x2+xpt(ngl(n))*phie(n);
                    y1=y1+ypt(ngl(n))*phic(n);
                    y2=y2+ypt(ngl(n))*phie(n);
                end
                dett=x1*y2-x2*y1;

                for i=1:9
                    tphx(i)=(y2*phic(i)-y1*phie(i))/dett;
                    tphy(i)=(x1*phie(i)-x2*phic(i))/dett;
                end
                
         temp1=u(ngl).*phi';
         temp2=u(ngl).*tphx';
         temp3=u(ngl).*tphy';
         temp1=sum(temp1);
         temp2=sum(temp2);
         temp3=sum(temp3);
         
                
                % *** Jij
                for l=1:9
                    for m=1:9
                       
                     J(ngl(l),ngl(m)) = J(ngl(l),ngl(m)) ...
                     -w(j)*w(k)*dett*(tphx(l)*tphx(m)+tphy(l)*tphy(m))...
                     +w(j)*w(k)*dett*( phi(l)*phi(m)/(1+alpha*temp1)^2 ...
                     * lambda*exp( temp1/(1+alpha*temp1) ) ) ;   
                                      
                    end
                    
                    r1(ngl(l)) = r1(ngl(l)) - w(j)*w(k)*dett*( lambda* phi(l)*exp( temp1/(1+alpha*temp1) ) )...
                        +w(j)*w(k)*dett*( tphx(l)*temp2+tphy(l)*temp3 );
                                     
                end
                
                
            end %LOOP k
        end %LOOP j
        
        
        
    end
end
