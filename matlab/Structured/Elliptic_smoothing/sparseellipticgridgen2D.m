function[lattice]=sparseellipticgridgen2D(Nx,N,lattice,deltaq,flagperiodicity,periodicity,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,spyflag)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 17th, 2014
%    Last update: August 7th, 2014
%
%          Input: meshed computational domain compdomain
%         Output: mesh in the physical domain

%%

boundaryindices = [];
periodicboundaryindices = [];
periodicrows = [];
periodiccolumns = [];
periodicvalues = [];

if ~flagperiodicity
    boundaryindices = [indicesE1;indicesE2;indicesE3;indicesE4;indicesC1;indicesC2;indicesC3;indicesC4];
elseif  any(periodicity==1) && ~any(periodicity==2) && ~any(periodicity==3) && ~any(periodicity==4)
    boundaryindices = [indicesE2;indicesE4;indicesC1;indicesC2;indicesC3;indicesC4];
    periodicboundaryindices = [indicesE1;indicesE3];
    periodicrows    = [indicesE1;indicesE1;  indicesE1;  indicesE1;   indicesE1;   indicesE1;     indicesE1;     indicesE1;     indicesE1;...
                       indicesE3;indicesE3;  indicesE3;  indicesE3;   indicesE3;   indicesE3;     indicesE3;     indicesE3;     indicesE3];
    periodiccolumns = [indicesE1;indicesE1-1;indicesE1+1;indicesE3;   indicesE1+Nx;indicesE3-1;   indicesE3+1;   indicesE1-1+Nx;indicesE1+1+Nx;...
                       indicesE3;indicesE3-1;indicesE3+1;indicesE3-Nx;indicesE1;   indicesE3-1-Nx;indicesE3+1-Nx;indicesE1-1;   indicesE1+1];
elseif ~any(periodicity==1) &&  any(periodicity==2) && ~any(periodicity==3) && ~any(periodicity==4)
    boundaryindices = [indicesE1;indicesE3;indicesC1;indicesC2;indicesC3;indicesC4];
    periodicboundaryindices = [indicesE2;indicesE4];
    periodicrows    = [indicesE2;indicesE2;  indicesE2;  indicesE2;   indicesE2;   indicesE2;     indicesE2;     indicesE2;     indicesE2;...
                       indicesE4;indicesE4;  indicesE4;  indicesE4;   indicesE4;   indicesE4;     indicesE4;     indicesE4;     indicesE4];
    periodiccolumns = [indicesE2;indicesE2-1;indicesE4;  indicesE2-Nx;indicesE2+Nx;indicesE2-1-Nx;indicesE4-Nx;  indicesE2-1+Nx;indicesE4+Nx;...
                       indicesE4;indicesE2;  indicesE4+1;indicesE4-Nx;indicesE4+Nx;indicesE2-Nx;  indicesE4+1-Nx;indicesE2+Nx;  indicesE4+1+Nx];
%elseif ~any(periodicity==1) && ~any(periodicity==2) &&  any(periodicity==3) && ~any(periodicity==4)
            
%elseif ~any(periodicity==1) && ~any(periodicity==2) && ~any(periodicity==3) &&  any(periodicity==4)
    
elseif  any(periodicity==1) &&  any(periodicity==2) && ~any(periodicity==3) && ~any(periodicity==4)
    boundaryindices = [indicesC1;indicesC2;indicesC3;indicesC4];
    periodicboundaryindices = [indicesE1;indicesE2;indicesE3;indicesE4];
elseif  any(periodicity==1) &&  any(periodicity==2) &&  any(periodicity==3) &&  any(periodicity==4)
    boundaryindices = [];
    periodicboundaryindices = [indicesE1;indicesE2;indicesE3;indicesE4;indicesC1;indicesC2;indicesC3;indicesC4];
end

effectivebulk = [];

for i=1:N
    if ~any(boundaryindices==i) && ~any(periodicboundaryindices==i)
        effectivebulk = [effectivebulk;i];
    end
end

it = 0;
err = 1;

if spyflag
    covariantbase = computecovariantbase2D(N,deltaq,lattice,firstdevneighbours);
    metriccoefficients = computemetriccoefficients2D(covariantbase); 
    clear covariantbase
    if ~flagperiodicity
        periodicvalues = [];
    elseif  any(periodicity==1) && ~any(periodicity==2) && ~any(periodicity==3) && ~any(periodicity==4)
        periodicvalues =  [-2*(metriccoefficients(indicesE1,2)./(deltaq(1).^2)+metriccoefficients(indicesE1,1)./(deltaq(2).^2));metriccoefficients(indicesE1,2)./(deltaq(1).^2);metriccoefficients(indicesE1,2)./(deltaq(1).^2);metriccoefficients(indicesE1,1)./(deltaq(2).^2);metriccoefficients(indicesE1,1)./(deltaq(2).^2);-0.5*metriccoefficients(indicesE1,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE1,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE1,3)./(deltaq(1).*deltaq(2));-0.5*metriccoefficients(indicesE1,3)./(deltaq(1).*deltaq(2));...
                           -2*(metriccoefficients(indicesE3,2)./(deltaq(1).^2)+metriccoefficients(indicesE3,1)./(deltaq(2).^2));metriccoefficients(indicesE3,2)./(deltaq(1).^2);metriccoefficients(indicesE3,2)./(deltaq(1).^2);metriccoefficients(indicesE3,1)./(deltaq(2).^2);metriccoefficients(indicesE3,1)./(deltaq(2).^2);-0.5*metriccoefficients(indicesE3,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE3,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE3,3)./(deltaq(1).*deltaq(2));-0.5*metriccoefficients(indicesE3,3)./(deltaq(1).*deltaq(2))];
    elseif ~any(periodicity==1) &&  any(periodicity==2) && ~any(periodicity==3) && ~any(periodicity==4)
        periodicvalues =  [-2*(metriccoefficients(indicesE2,2)./(deltaq(1).^2)+metriccoefficients(indicesE2,1)./(deltaq(2).^2));metriccoefficients(indicesE2,2)./(deltaq(1).^2);metriccoefficients(indicesE2,2)./(deltaq(1).^2);metriccoefficients(indicesE2,1)./(deltaq(2).^2);metriccoefficients(indicesE2,1)./(deltaq(2).^2);-0.5*metriccoefficients(indicesE2,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE2,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE2,3)./(deltaq(1).*deltaq(2));-0.5*metriccoefficients(indicesE2,3)./(deltaq(1).*deltaq(2));...
                           -2*(metriccoefficients(indicesE4,2)./(deltaq(1).^2)+metriccoefficients(indicesE4,1)./(deltaq(2).^2));metriccoefficients(indicesE4,2)./(deltaq(1).^2);metriccoefficients(indicesE4,2)./(deltaq(1).^2);metriccoefficients(indicesE4,1)./(deltaq(2).^2);metriccoefficients(indicesE4,1)./(deltaq(2).^2);-0.5*metriccoefficients(indicesE4,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE4,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE4,3)./(deltaq(1).*deltaq(2));-0.5*metriccoefficients(indicesE4,3)./(deltaq(1).*deltaq(2))];
    %elseif ~any(periodicity==1) && ~any(periodicity==2) &&  any(periodicity==3) && ~any(periodicity==4)
            
    %elseif ~any(periodicity==1) && ~any(periodicity==2) && ~any(periodicity==3) &&  any(periodicity==4)
    
    elseif  any(periodicity==1) &&  any(periodicity==2) && ~any(periodicity==3) && ~any(periodicity==4)
    
    elseif  any(periodicity==1) &&  any(periodicity==2) &&  any(periodicity==3) &&  any(periodicity==4)

    end
    A = sparse([boundaryindices;periodicrows;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk],[boundaryindices;periodiccolumns;effectivebulk;effectivebulk-1;effectivebulk+1;effectivebulk-Nx;effectivebulk+Nx;effectivebulk-1-Nx;effectivebulk+1-Nx;effectivebulk-1+Nx;effectivebulk+1+Nx],[ones(size(boundaryindices,1),1);periodicvalues;-2*(metriccoefficients(effectivebulk,2)./(deltaq(1).^2)+metriccoefficients(effectivebulk,1)./(deltaq(2).^2));metriccoefficients(effectivebulk,2)./(deltaq(1).^2);metriccoefficients(effectivebulk,2)./(deltaq(1).^2);metriccoefficients(effectivebulk,1)./(deltaq(2).^2);metriccoefficients(effectivebulk,1)./(deltaq(2).^2);-0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2));-0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2))],N,N); 
    periodicvalues = [];
    figure;
    spy(A)
    hold on
    grid on
    title('Structure of solving matrix for 2D elliptic grid generation')
    clear A metriccoefficients
end

while it<=itmax && err>=tol
    covariantbase = computecovariantbase2D(N,deltaq,lattice,firstdevneighbours);
    [metriccoefficients,J] = computemetriccoefficients2D(covariantbase);
    Pvec = sparse(P(lattice));
    Qvec = sparse(Q(lattice));
    bx = sparse([boundaryindices;periodicrows;effectivebulk],[ones(length(boundaryindices),1);ones(length(periodicrows),1);ones(length(effectivebulk),1)],[lattice(boundaryindices,3);-(J(periodicrows,:).^2).*(Pvec(periodicrows,:).*covariantbase(periodicrows,1)+Qvec(periodicrows,:).*covariantbase(periodicrows,3));-(J(effectivebulk,:).^2).*(Pvec(effectivebulk,:).*covariantbase(effectivebulk,1)+Qvec(effectivebulk,:).*covariantbase(effectivebulk,3))],N,1);
    by = sparse([boundaryindices;periodicrows;effectivebulk],[ones(length(boundaryindices),1);ones(length(periodicrows),1);ones(length(effectivebulk),1)],[lattice(boundaryindices,4);-(J(periodicrows,:).^2).*(Pvec(periodicrows,:).*covariantbase(periodicrows,2)+Qvec(periodicrows,:).*covariantbase(periodicrows,4));-(J(effectivebulk,:).^2).*(Pvec(effectivebulk,:).*covariantbase(effectivebulk,2)+Qvec(effectivebulk,:).*covariantbase(effectivebulk,4))],N,1);
    clear covariantbase Pvec Qvec J
    if ~flagperiodicity
        periodicvalues = [];
    elseif  any(periodicity==1) && ~any(periodicity==2) && ~any(periodicity==3) && ~any(periodicity==4)
        periodicvalues =  [-2*(metriccoefficients(indicesE1,2)./(deltaq(1).^2)+metriccoefficients(indicesE1,1)./(deltaq(2).^2));metriccoefficients(indicesE1,2)./(deltaq(1).^2);metriccoefficients(indicesE1,2)./(deltaq(1).^2);metriccoefficients(indicesE1,1)./(deltaq(2).^2);metriccoefficients(indicesE1,1)./(deltaq(2).^2);-0.5*metriccoefficients(indicesE1,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE1,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE1,3)./(deltaq(1).*deltaq(2));-0.5*metriccoefficients(indicesE1,3)./(deltaq(1).*deltaq(2));...
                           -2*(metriccoefficients(indicesE3,2)./(deltaq(1).^2)+metriccoefficients(indicesE3,1)./(deltaq(2).^2));metriccoefficients(indicesE3,2)./(deltaq(1).^2);metriccoefficients(indicesE3,2)./(deltaq(1).^2);metriccoefficients(indicesE3,1)./(deltaq(2).^2);metriccoefficients(indicesE3,1)./(deltaq(2).^2);-0.5*metriccoefficients(indicesE3,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE3,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE3,3)./(deltaq(1).*deltaq(2));-0.5*metriccoefficients(indicesE3,3)./(deltaq(1).*deltaq(2))];
    elseif ~any(periodicity==1) &&  any(periodicity==2) && ~any(periodicity==3) && ~any(periodicity==4)
        periodicvalues =  [-2*(metriccoefficients(indicesE2,2)./(deltaq(1).^2)+metriccoefficients(indicesE2,1)./(deltaq(2).^2));metriccoefficients(indicesE2,2)./(deltaq(1).^2);metriccoefficients(indicesE2,2)./(deltaq(1).^2);metriccoefficients(indicesE2,1)./(deltaq(2).^2);metriccoefficients(indicesE2,1)./(deltaq(2).^2);-0.5*metriccoefficients(indicesE2,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE2,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE2,3)./(deltaq(1).*deltaq(2));-0.5*metriccoefficients(indicesE2,3)./(deltaq(1).*deltaq(2));...
                           -2*(metriccoefficients(indicesE4,2)./(deltaq(1).^2)+metriccoefficients(indicesE4,1)./(deltaq(2).^2));metriccoefficients(indicesE4,2)./(deltaq(1).^2);metriccoefficients(indicesE4,2)./(deltaq(1).^2);metriccoefficients(indicesE4,1)./(deltaq(2).^2);metriccoefficients(indicesE4,1)./(deltaq(2).^2);-0.5*metriccoefficients(indicesE4,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE4,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(indicesE4,3)./(deltaq(1).*deltaq(2));-0.5*metriccoefficients(indicesE4,3)./(deltaq(1).*deltaq(2))];
    %elseif ~any(periodicity==1) && ~any(periodicity==2) &&  any(periodicity==3) && ~any(periodicity==4)
            
    %elseif ~any(periodicity==1) && ~any(periodicity==2) && ~any(periodicity==3) &&  any(periodicity==4)
    
    elseif  any(periodicity==1) &&  any(periodicity==2) && ~any(periodicity==3) && ~any(periodicity==4)
    
    elseif  any(periodicity==1) &&  any(periodicity==2) &&  any(periodicity==3) &&  any(periodicity==4)

    end
    A = sparse([boundaryindices;periodicrows;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk],[boundaryindices;periodiccolumns;effectivebulk;effectivebulk-1;effectivebulk+1;effectivebulk-Nx;effectivebulk+Nx;effectivebulk-1-Nx;effectivebulk+1-Nx;effectivebulk-1+Nx;effectivebulk+1+Nx],[ones(size(boundaryindices,1),1);periodicvalues;-2*(metriccoefficients(effectivebulk,2)./(deltaq(1).^2)+metriccoefficients(effectivebulk,1)./(deltaq(2).^2));metriccoefficients(effectivebulk,2)./(deltaq(1).^2);metriccoefficients(effectivebulk,2)./(deltaq(1).^2);metriccoefficients(effectivebulk,1)./(deltaq(2).^2);metriccoefficients(effectivebulk,1)./(deltaq(2).^2);-0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2));-0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2))],N,N); 
    periodicvalues = [];
    lattice(:,5) = full(A\bx);
    lattice(:,6) = full(A\by);
    err = sqrt(sum([bx-A*lattice(:,5);by-A*lattice(:,6)].^2));
    it = it + 1;
    clear A bx by metriccoefficients
end

return