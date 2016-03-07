function[lattice]=transfiniteinterpolation2D(N,compdomain,xi_min,xi_max,eta_min,eta_max,Ndim1,Ndim2,interpolanttype,e1,e2,e3,e4,c1,c2,c3,c4)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 15th, 2014
%    Last update: July 10th, 2014
%
%          Input: meshed computational domain compdomain
%                 numerical flag to choose the interpolant type
%                 e1 = e1(xi,eta_min) = e1(xi) = (x(xi),y(xi))
%                 e2 = e2(xi_min,eta) = e2(eta) = (x(eta),y(eta))
%                 e3 = e3(xi,eta_max) = e3(xi) = (x(xi),y(xi))
%                 e4 = e4(xi_max,eta) = e4(eta) = (x(eta),y(eta))
%                 c1 = c1(xi_min,eta_min) = (x1,y1)
%                 c2 = c2(xi_max,eta_min) = (x2,y2)
%                 c3 = c3(xi_max,eta_max) = (x3,y3)
%                 c4 = c4(xi_min,eta_max) = (x4,y4)
%                 e1 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e2 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e3 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e4 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 c1 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c2 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c3 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c4 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%
%    Interpolant: 1 --> Lagrange
%                 2 --> Hermite
%
%         Output: mesh in the physical domain

%%

mesh = zeros(N,2);

switch interpolanttype
    case 1 %--------------------------- Lagrange
        for r=1:Ndim2
            for p=1:Ndim1
                xi  = compdomain(p + (r-1)*Ndim1,1);
                eta = compdomain(p + (r-1)*Ndim1,2);
                mesh(p + (r-1)*Ndim1,1) = Lagrangeinterps1D([xi_min;xi_max],[e2(r,1);e4(r,1)],xi) + Lagrangeinterps1D([eta_min;eta_max],[e1(p,1);e3(p,1)],eta) - Lagrangeinterps2D([xi_min;xi_max],[eta_min;eta_max],[c1(1,1) c4(1,1);c2(1,1) c3(1,1)],xi,eta);
                mesh(p + (r-1)*Ndim1,2) = Lagrangeinterps1D([xi_min;xi_max],[e2(r,2);e4(r,2)],xi) + Lagrangeinterps1D([eta_min;eta_max],[e1(p,2);e3(p,2)],eta) - Lagrangeinterps2D([xi_min;xi_max],[eta_min;eta_max],[c1(1,2) c4(1,2);c2(1,2) c3(1,2)],xi,eta); 
            end
        end
    case 2 %--------------------------- Hermite
        for r=1:Ndim2
            for p=1:Ndim1
                xi  = compdomain(p + (r-1)*Ndim1,1);
                eta = compdomain(p + (r-1)*Ndim1,2);
                mesh(p + (r-1)*Ndim1,1) = Hermiteinterps1D([xi_min;xi_max],[e2(r,1);e4(r,1)],[e2(r,3);e4(r,3)],xi) + Hermiteinterps1D([eta_min;eta_max],[e1(p,1);e3(p,1)],[e2(r,4);e4(r,4)],eta) - Hermiteinterps2D([xi_min;xi_max],[eta_min;eta_max],[c1(1,1) c4(1,1);c2(1,1) c3(1,1)],[c1(1,3) c4(1,3);c2(1,3) c3(1,3)],[c1(1,4) c4(1,4);c2(1,4) c3(1,4)],[c1(1,7) c4(1,7);c2(1,7) c3(1,7)],xi,eta);
                mesh(p + (r-1)*Ndim1,2) = Hermiteinterps1D([xi_min;xi_max],[e2(r,2);e4(r,2)],[e2(r,5);e4(r,5)],xi) + Hermiteinterps1D([eta_min;eta_max],[e1(p,2);e3(p,2)],[e2(r,6);e4(r,6)],eta) - Hermiteinterps2D([xi_min;xi_max],[eta_min;eta_max],[c1(1,2) c4(1,2);c2(1,2) c3(1,2)],[c1(1,5) c4(1,5);c2(1,5) c3(1,5)],[c1(1,6) c4(1,6);c2(1,6) c3(1,6)],[c1(1,8) c4(1,8);c2(1,8) c3(1,8)],xi,eta); 
            end
        end
end

[indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(Ndim1,Ndim2);

mesh(indicesC1,1) = c1(1,1);
mesh(indicesC1,2) = c1(1,2);

mesh(indicesC2,1) = c2(1,1);
mesh(indicesC2,2) = c2(1,2);

mesh(indicesC3,1) = c3(1,1);
mesh(indicesC3,2) = c3(1,2);

mesh(indicesC4,1) = c4(1,1);
mesh(indicesC4,2) = c4(1,2);

mesh(indicesE1,1) = e1(2:end-1,1);
mesh(indicesE1,2) = e1(2:end-1,2);

mesh(indicesE2,1) = e4(2:end-1,1);
mesh(indicesE2,2) = e4(2:end-1,2);

mesh(indicesE3,1) = e3(2:end-1,1);
mesh(indicesE3,2) = e3(2:end-1,2);

mesh(indicesE4,1) = e2(2:end-1,1);
mesh(indicesE4,2) = e2(2:end-1,2);

lattice = [compdomain mesh mesh];

return