function[]=meshGenCurvHollowSec(section)
%%
%==============================================================================
% Copyright (c) 2016 Universit� de Lorraine & Lule� tekniska universitet
% Author: Luca Di Stasio <luca.distasio@gmail.com>
%                        <luca.distasio@ingpec.eu>
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%
% Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
% Neither the name of the Universit� de Lorraine or Lule� tekniska universitet
% nor the names of its contributors may be used to endorse or promote products
% derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%==============================================================================
%
%  DESCRIPTION
%
%  A function to generate a structured mesh of a hollow (one hole) section with a generic
%  curvilinear boundary
%
%  Output:
%
%%

% assign properties to local variables
c1 = section.c1;
c2 = section.c2;
c3 = section.c3;
c4 = section.c4;
e1 = section.e1;
e2 = section.e2;
e3 = section.e3;
e4 = section.e4;

N1 = length(section.e1)+1;
N2 = length(section.e2)+1;

itmax = section.itmax;
tol = section.tol;

% function call
% lattice=generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny)
mesh = generatelattice2D(c1(1),N1+1,(c2(1)-c1(1))/N1,c1(2),N2+1,(c4(2)-c1(2))/N2,c1(1),c2(1),N1+1,c1(2),c4(2),N2+1);

%function call
%lattice=transfiniteinterpolation2D(N,compdomain,xi_min,xi_max,eta_min,eta_max,Ndim1,Ndim2,interpolanttype,e1,e2,e3,e4,c1,c2,c3,c4)
mesh = transfiniteinterpolation2D((N1+1)*(N2+1),mesh(:,1:2),c1(1),c2(1),c1(2),c4(2),N1+1,N2+1,1,e1,e2,e3,e4,c1,c2,c3,c4);

[indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(N1+1,N2+1);
[temp1,temp2,temp3,firstdevneighbours] = build_neighbourhoods2D((N1+1)*(N2+1),N1+1,0,0,0,0,0,indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);
%function call
%lattice=sparseellipticgridgen2D(Nx,N,lattice,deltaq,flagperiodicity,periodicity,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,spyflag)
mesh = sparseellipticgridgen2D(N1+1,(N1+1)*(N2+1),mesh,[(c2(1)-c1(1))/N1,(c4(2)-c1(2))/N2],0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0);

section.mesh = mesh;

return
