function[indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4]=getindices2D(Nx,Ny)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 10th, 2014
%    Last update: July 11th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

% ---> compute indices

% bulk fluid

indicesbulk = zeros((Nx-2)*(Ny-2),1);


for j=1:Ny-2
    for i=1:Nx-2
        indicesbulk(i+(j-1)*(Nx-2),1) = (i+1)+j*Nx;
    end
end


indicesinternalbulk = zeros((Nx-4)*(Ny-4),1);


for j=1:Ny-4
    for i=1:Nx-4
        indicesinternalbulk(i+(j-1)*(Nx-4),1) = (i+2)+(j+1)*Nx;
    end
end

% edges

indicesE1 = zeros(Nx-2,1); % x

for i=1:Nx-2
    j = 1;
    indicesE1(i,1) = (i+1)+(j-1)*Nx;
end

indicesE2 = zeros(Ny-2,1); % y

for j=1:Ny-2
    i = Nx;
    indicesE2(j,1) = i+j*Nx;
end

indicesE3 = zeros(Nx-2,1); % x

for i=1:Nx-2
    j = Ny;
    indicesE3(i,1) = (i+1)+(j-1)*Nx;
end

indicesE4 = zeros(Ny-2,1); % y

for j=1:Ny-2
    i = 1;
    indicesE4(j,1) = i+j*Nx;
end

indicesinternalE1 = zeros(Nx-4,1); % x

for i=1:Nx-4
    j = 2;
    indicesinternalE1(i,1) = (i+2)+(j-1)*Nx;
end

indicesinternalE2 = zeros(Ny-4,1); % y

for j=1:Ny-4
    i = Nx-1;
    indicesinternalE2(j,1) = i+(j+1)*Nx;
end

indicesinternalE3 = zeros(Nx-4,1); % x

for i=1:Nx-4
    j = Ny-1;
    indicesinternalE3(i,1) = (i+2)+(j-1)*Nx;
end

indicesinternalE4 = zeros(Ny-4,1); % y

for j=1:Ny-4
    i = 2;
    indicesinternalE4(j,1) = i+(j+1)*Nx;
end

indicesexternalE1 = indicesE1(indicesE1(indicesE1~=2)~=Nx-1);

indicesexternalE2 = indicesE2(indicesE2(indicesE2~=2*Nx)~=Nx+(Ny-2)*Nx);

indicesexternalE3 = indicesE3(indicesE3(indicesE3~=Nx*Ny-1)~=2+(Ny-1)*Nx);

indicesexternalE4 = indicesE4(indicesE4(indicesE4~=1+Nx)~=1+(Ny-2)*Nx);

% corners

indicesC1 = 1;

indicesC2 = Nx;

indicesC3 = Nx + Nx*(Ny-1);

indicesC4 = 1 + Nx*(Ny-1);

indicesinternalC1 = [2 + Nx + Nx*Ny;2;1+Nx];

indicesinternalC2 = [(Nx-1) + Nx + Nx*Ny;Nx-1;2*Nx];

indicesinternalC3 = [Nx*Ny - Nx - 1 + Nx*Ny;Nx+(Ny-2)*Nx;Nx*Ny-1];

indicesinternalC4 = [2 + Nx*(Ny-2) + Nx*Ny;1+(Ny-2)*Nx;2+(Ny-1)*Nx];

return