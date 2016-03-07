function[u1]=JFNKmethod(F,u0,tol,itmax,precond,alpha,sigma0,sigma1)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 19th, 2014
%    Last update: May 21st, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

N = size(u0,1);

if precond
    epsdiff = sqrt(1+eps(sqrt(sum(u0.^2))))/sqrt(sum(u0.^2));
    M = (feval(F,u0+epsdiff)-feval(u0-epsdiff))/(2*epsdiff);
    invM = luinv(M);
else
    invM = eye(N);
end

F0 = feval(F,u0);
errnewton = F0;
k = 0;

while errnewton>=tol && k<=itmax  
    du0 = zeros(N,1);
    %------ GMRES(m) -----
    err = 1;
    restart = 0;
    while err>=tol && restart<=maxrestart
        v = zeros(N,m+1);
        h = zeros(m+1,m);
        if restart == 0
            r0 = -F0;
        else
            epsdiff = sqrt(1+eps(sqrt(sum(u0.^2))))/sqrt(sum(du0.^2));
            r0 = -F0 -(feval(F,u0+epsdiff*du0)-feval(u0-epsdiff*du0))/(2*epsdiff);
        end
        beta = sqrt(sum(r0.^2));
        v(:,1) = r0/beta;
        j = 1;
        innerloop = 1;
        while innerloop && j<=m
            zj = invM*v(:,j);
            epsdiff = sqrt(1+eps(sqrt(sum(u0.^2))))/sqrt(sum(zj.^2));
            wj = (feval(F,u0+epsdiff*zj)-feval(u0-epsdiff*zj))/(2*epsdiff);
            for i=1:j
                h(i,j) = wj'*v(:,i);
                wj = wj -h(i,j)*v(:,i);
            end
            h(j+1,j) = sqrt(sum(wj.^2));
            if h(j+1,j)==0
                innerloop = 0;
            else
                v(:,j+1) = wj/h(j+1,j);
                if j==m
                    innerloop = 0;
                else
                    j = j + 1;
                end
            end
        end
        if j~=m
            v = v(:,1:j);
            h = h(1:j+1,1:j);
        end
        r = zeros(j+1,1);
        r(1,1) = beta;
        for k=1:j
            P = eye(j+1);
            s = h(k+1,k)/(sqrt(h(k+1,k)^2+h(k,k)^2));
            c = h(k,k)/(sqrt(h(k+1,k)^2+h(k,k)^2));
            P(k,k) = c;
            P(k+1,k+1) = c;
            P(k,k+1) = s;
            P(k+1,k) = -s;
            h = P*h;
            r = P*r;
        end
        y = backsolve(h(1:j,:),r(1:j,:));
        duk = du0 + v(:,1:j)*y;
        err = abs(r(j+1,1));
        du0 = duk;
        restart = restart + 1;
    end
    lambda = 1;
    u1trial = u0 + lambda*duk;
    while sqrt(sum(feval(F,u1trial).^2))>=(1-alpha*lambda)*sqrt(sum(feval(F,u0).^2))
        sigma = sigma0 + rand(1)*(sigma1-sigma0);
        lambda = sigma*lambda;
        u1trial = u0 + lambda*duk;
    end
    u1 = u1trial;
    u0 = u1;
    F0 = feval(F,u0);
    errnewton = F0;
    k = k + 1;
end

return