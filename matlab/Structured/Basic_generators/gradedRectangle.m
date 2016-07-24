function[mesh] = gradedRectangle(x0,y0,lx,ly,Nx,Ny)

%%
%==============================================================================
% Copyright (c) 2016 Université de Lorraine & Luleå tekniska universitet
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
% Neither the name of the Université de Lorraine or Luleå tekniska universitet
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
%  A function to create a regular mesh of quadrilaterals in a rectangular
%  geometry with n zones with different element size
%
%  Input: x0 - scalar - x-coordinate of center
%         y0 - scalar - y-coordinate of center
%         lx - [M x 1] vector - Side length in each one of the M mesh regions in x-direction
%         ly - [N x 1] vector - Side length in each one of the N mesh regions in y-direction
%         Nx - [M x 1] vector - Number of ELEMENTS in each one of the M mesh regions in x-direction
%         Ny - [N x 1] vector - Number of ELEMENTS in each one of the N mesh regions in y-direction
%
%  Output: mesh - (sum(Nx)+1)*(sum(Ny)+1) x 2 matrix - mesh nodes ordered through helical
%          indexing
%
%%

M = size(Nx,1);
N = size(Ny,1);

deltax = lx./Nx;
deltay = ly./Ny;

nodesx = zeros(M+1,1);
nodesy = zeros(N+1,1);

nodesx(1) = x0-0.5*sum(lx);
nodesy(1) = y0-0.5*sum(ly);

for m=2:M+1
   nodesx(m) = nodesx(m-1) + lx(m-1);  
end

for n=2:N+1
   nodesy(n) = nodesy(n-1) + ly(n-1);  
end

xs = [];
ys = [];

for m=1:M
    addVec = (nodesx(m):deltax(m):nodesx(m+1))';
    xs = [xs;...
         addVec(1:end-1)];
end
xs = [xs;...
      x0+0.5*sum(lx)];
     
for n=1:N
    addVec = (nodesy(n):deltay(n):nodesy(n+1))';
    ys = [ys;...
         addVec(1:end-1)]; 
end
ys = [ys;...
      y0+0.5*sum(ly)];

mesh = zeros((sum(Nx)+1)*(sum(Ny)+1),2);

for j=1:(sum(Ny)+1)
    mesh((j-1)*(sum(Nx)+1)+1:j*(sum(Nx)+1),1) = xs;
    mesh((j-1)*(sum(Nx)+1)+1:j*(sum(Nx)+1),2) = ys(j)*ones((sum(Nx)+1),1);
end

return