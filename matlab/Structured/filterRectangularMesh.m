function[nodes,elements,edges,...
         nodesSWcorner,nodesSEcorner,nodesNEcorner,nodesNWcorner,nodesSOUTHside,nodesEASTside,nodesNORTHside,nodesWESTside,...
         edgesSOUTHside,edgesEASTside,edgesNORTHside,edgesWESTside,...
         elementsSWcorner,elementsSEcorner,elementsNEcorner,elementsNWcorner,elementsSOUTHside,elementsEASTside,elementsNORTHside,elementsWESTside]=filterRectangularMesh(logfullfile,latexFolder,...
                                                                                                                                                                          elType,elOrder,Nx,Ny,NxEquiv,NyEquiv,...
                                                                                                                                                                          baseNodes,isCircular,circularity)
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
%  A function to mesh a simply connected 2D rectangular geometry with elements of
%  shape and order of choice
%
%  Available:
%  - 1st and 2nd order quads
%  - 1st and 2nd order tris
%
%  Input: x0 - scalar - x-coordinate of center
%         y0 - scalar - y-coordinate of center
%         lx - [M x 1] vector - Side length in each one of the M mesh regions in x-direction
%         ly - [N x 1] vector - Side length in each one of the N mesh regions in y-direction
%         Nx - [M x 1] vector - Number of ELEMENTS in each one of the M mesh regions in x-direction
%         Ny - [N x 1] vector - Number of ELEMENTS in each one of the N mesh regions in y-direction
%
%         For
%            --> 1st order quadrilaterals
%                    NxEquiv = Nx
%                    NyEquiv = Ny
%            --> 2nd order quadrilaterals
%                    NxEquiv = 2*Nx
%                    NyEquiv = 2*Ny
%            --> 1st order triangles
%                    NxEquiv = Nx
%                    NyEquiv = Ny
%            --> 2nd order triangles
%                    NxEquiv = 2*Nx
%                    NyEquiv = 2*Ny
%%


writeToLogFile(logfullfile,'In function: filterRectangularMesh\n')
writeToLogFile(logfullfile,'\nStarting timer\n')
start = tic;
% filter nodes if elements are not linear quadrilaterals
writeToLogFile(logfullfile,'Filtering nodes if elements are not linear quadrilaterals and building connectivity of edges and elements ...\n')
try
  if strcomp(elType,'quads') || strcomp(elType,'quad') || strcomp(elType,'quadrilaterals') || strcomp(elType,'quadrilateral')
    if strcomp(elType,'first') || strcomp(elType,'First') || strcomp(elType,'1st') || strcomp(elType,'1')
      if isCircular
        if strcomp(circularity,'EW') || strcomp(circularity,'E-W') || strcomp(circularity,'EastWest') || strcomp(circularity,'East-West')
          filteredNodes = zeros(sum(Nx)*(sum(Ny)+1),2);
          edges = zeros((2*sum(Nx))*sum(Ny)+sum(Nx),2);
          elements = zeros(sum(Nx)*sum(Ny),4);
          for j=1:sum(Ny)
            filteredNodes((j-1)*sum(Nx)+1:j*sum(Nx),1:2) = baseNodes((j-1)*(sum(Nx)+1)+1:(j-1)*(sum(Nx)+1)+sum(Nx),1:2);
            elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
            elements((j-1)*sum(Nx)+1:j*sum(Nx)-1,2) = ((j-1)*sum(Nx)+2:j*sum(Nx))';
            elements((j-1)*sum(Nx)+1:j*sum(Nx)-1,3) = (j*sum(Nx)+2:(j+1)*sum(Nx))';
            elements(j*sum(Nx),2) = (j-1)*sum(Nx)+1;
            elements(j*sum(Nx),3) = j*sum(Nx)+1;
            elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
            edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
            edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
            edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
            edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),2) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
          end
          j = sum(Ny)+1;
          filteredNodes((j-1)*sum(Nx)+1:j*sum(Nx)) = baseNodes((j-1)*(sum(Nx)+1)+1:(j-1)*(sum(Nx)+1)+sum(Nx),1:2);
          edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),1) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
          edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),2) = (j*sum(Nx)+2:(j+1)*sum(Nx)+1)';
          nodes = filteredNodes;
          nodesSWcorner = -1;
          nodesSEcorner = -1;
          nodesNEcorner = -1;
          nodesNWcorner = -1;
          nodesSOUTHside = (1:sum(Nx))'; %except corners
          nodesEASTside = []; %except corners
          nodesNORTHside = sum(Ny)*sum(Nx)+(1:sum(Nx))'; %except corners
          nodesWESTside = []; %except corners
          edgesSOUTHside = (1:sum(Nx))';
          edgesEASTside = [];
          edgesNORTHside = ((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx))';
          edgesWESTside = [];
          elementsSWcorner = -1;
          elementsSEcorner = -1;
          elementsNEcorner = -1;
          elementsNWcorner = -1;
          elementsSOUTHside = (1:sum(Nx))'; %except corners
          elementsEASTside = []; %except corners
          elementsNORTHside = sum(Ny)*sum(Nx)+(1:sum(Nx))'; %except corners
          elementsWESTside = []; %except corners
        else % North-South circularity
          filteredNodes = zeros((sum(Nx)+1)*sum(Ny),2);
          edges = zeros((2*sum(Nx)+1)*sum(Ny),2);
          elements = zeros(sum(Nx)*sum(Ny),4);
          for j=1:sum(Ny)-1
            filteredNodes((j-1)*(sum(Nx)+1)+1:j*(sum(Nx)+1),1:2) = baseNodes((j-1)*(sum(Nx)+1)+1:j*(sum(Nx)+1),1:2);
            elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),3) = (j*sum(Nx)+2:(j+1)*sum(Nx)+1)';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
            edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
            edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
            edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
            edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),2) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
            edges(2*j*sum(Nx)+1,1) = j*sum(Nx)+1;
            edges(2*j*sum(Nx)+1,2) = (j+1)*sum(Nx)+1;
          end
          j = sum(Ny);
          filteredNodes((j-1)*(sum(Nx)+1)+1:j*(sum(Nx)+1),1:2) = baseNodes((j-1)*(sum(Nx)+1)+1:j*(sum(Nx)+1),1:2);
          elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),3) = (2:sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (1:sum(Nx))';
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),2) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
          edges(2*j*sum(Nx)+1,1) = j*sum(Nx)+1;
          edges(2*j*sum(Nx)+1,2) = (j+1)*sum(Nx)+1;
          nodes = filteredNodes;
          nodesSWcorner = -1;
          nodesSEcorner = -1;
          nodesNEcorner = -1;
          nodesNWcorner = -1;
          nodesSOUTHside = []; %except corners
          nodesEASTside = (sum(Nx)+1:sum(Nx)+1:(sum(Ny)-1)*(sum(Nx)+1)+sum(Nx)+1)'; %except corners
          nodesNORTHside = []; %except corners
          nodesWESTside = (1:sum(Nx)+1:(sum(Ny)-1)*(sum(Nx)+1)+1)'; %except corners
          edgesSOUTHside = [];
          edgesEASTside = (2*sum(Nx)+1:2*sum(Nx)+1:(sum(Ny)-1)*(2*sum(Nx)+1)+2*sum(Nx)+1)';
          edgesNORTHside = [];
          edgesWESTside = (sum(Nx)+1:2*sum(Nx)+1:(sum(Ny)-1)*(2*sum(Nx)+1)+sum(Nx)+1)';
          elementsSWcorner = -1;
          elementsSEcorner = -1;
          elementsNEcorner = -1;
          elementsNWcorner = -1;
          elementsSOUTHside = []; %except corners
          elementsEASTside = (sum(Nx):sum(Nx):(sum(Ny)-1)*sum(Nx)+sum(Nx))'; %except corners
          elementsNORTHside = []; %except corners
          elementsWESTside = (1:sum(Nx):(sum(Ny)-1)*sum(Nx)+1)'; %except corners
        end
      else
        edges = zeros((2*sum(Nx)+1)*sum(Ny)+sum(Nx),2);
        elements = zeros(sum(Nx)*sum(Ny),4);
        nodes = baseNodes;
        for j=1:sum(Ny)
          elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),3) = (j*sum(Nx)+2:(j+1)*sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),2) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
          edges(2*j*sum(Nx)+1,1) = j*sum(Nx)+1;
          edges(2*j*sum(Nx)+1,2) = (j+1)*sum(Nx)+1;
        end
        j = sum(Ny)+1;
        edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),1) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
        edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),2) = (j*sum(Nx)+2:(j+1)*sum(Nx)+1)';
        nodesSWcorner = 1;
        nodesSEcorner = sum(Nx)+1;
        nodesNEcorner = (sum(Nx)+1)*(sum(Ny)+1);
        nodesNWcorner = (sum(Nx)+1)*sum(Ny)+1;
        nodesSOUTHside = (2:sum(Nx))'; %except corners
        nodesEASTside = (2*(sum(Nx)+1):sum(Nx)+1:sum(Ny)*(sum(Nx)+1))'; %except corners
        nodesNORTHside = ((sum(Nx)+1)*sum(Ny)+2:(sum(Nx)+1)*(sum(Ny)+1)-1)'; %except corners
        nodesWESTside = (sum(Nx)+1+1:sum(Nx)+1:(sum(Ny)-1)*(sum(Nx)+1)+1)'; %except corners
        edgesSOUTHside = (1:sum(Nx))';
        edgesEASTside = (2*sum(Nx)+1:2*sum(Nx)+1:sum(Ny)*(2*sum(Nx)+1))';
        edgesNORTHside = ((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx))';
        edgesWESTside = (sum(Nx)+1:2*sum(Nx)+1:(sum(Ny)-1)*(2*sum(Nx)+1)+sum(Nx)+1)';
        elementsSWcorner = 1;
        elementsSEcorner = sum(Nx);
        elementsNEcorner = sum(Nx)*sum(Ny);
        elementsNWcorner = sum(Nx)*(sum(Ny)-1)+1;
        elementsSOUTHside = (2:sum(Nx)-1)'; %except corners
        elementsEASTside = (2*sum(Nx):sum(Nx):sum(Nx)*(sum(Ny)-1))'; %except corners
        elementsNORTHside = (sum(Nx)*(sum(Ny)-1)+2:sum(Nx)*sum(Ny)-1)'; %except corners
        elementsWESTside = (sum(Nx)+1:sum(Nx):(sum(Ny)-2)*sum(Nx)+1)'; %except corners
      end
    elseif strcomp(elType,'second') || strcomp(elType,'Second') || strcomp(elType,'2nd') || strcomp(elType,'2')
      if isCircular
        if strcomp(circularity,'EW') || strcomp(circularity,'E-W') || strcomp(circularity,'EastWest') || strcomp(circularity,'East-West')

        else % North-South circularity

        end
      else
        filteredNodes = zeros((3*sum(Nx)+2)*sum(Ny)+2*sum(Nx)+1,2);
        edges = zeros((2*sum(Nx)+1)*sum(Ny)+sum(Nx),3);
        elements = zeros(sum(Nx)*sum(Ny),8);
        for j=1:sum(Ny)
          filteredNodes((j-1)*(3*sum(Nx)+2)+1:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1,1:2) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+1:2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv)+1,1:2);
          filteredNodes((j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+1:j*(3*sum(Nx)+2),1:2) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv)+1+1:2:2*j*(sum(NxEquiv)+1),1:2);
          elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),2) = ((j-1)*(3*sum(Nx)+2)+3:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),3) = (j*(3*sum(Nx)+2)+3:2:j*(3*sum(Nx)+2)+2*sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (j*(3*sum(Nx)+2)+1:2:j*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),5) = ((j-1)*(3*sum(Nx)+2)+2:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx))';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),6) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+2:2:j*(3*sum(Nx)+2);
          elements((j-1)*sum(Nx)+1:j*sum(Nx),7) = (j*(3*sum(Nx)+2)+2:2:j*(3*sum(Nx)+2)+2*sum(Nx))';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),8) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+1:2:j*(3*sum(Nx)+1);
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),2) = ((j-1)*(3*sum(Nx)+2)+2:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx))';
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),3) = ((j-1)*(3*sum(Nx)+2)+3:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1)';
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),2) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+1:2:j*(3*sum(Nx)+1);
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),3) = (j*(3*sum(Nx)+2)+1:2:j*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          edges(2*j*sum(Nx)+1,1) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1;
          edges(2*j*sum(Nx)+1,2) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx);
          edges(2*j*sum(Nx)+1,3) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1;
        end
        j = sum(Ny)+1;
        filteredNodes((j-1)*(3*sum(Nx)+2)+1:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1,1:2) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+1:2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv)+1,1:2);
        edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
        edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),2) = ((j-1)*(3*sum(Nx)+2)+2:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx))';
        edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),3) = ((j-1)*(3*sum(Nx)+2)+3:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1)';
        nodes = filteredNodes;
      end
    end
  elseif strcomp(elType,'tris') || strcomp(elType,'tri') || strcomp(elType,'triangles') || strcomp(elType,'triangle')
    if strcomp(elType,'first') || strcomp(elType,'First') || strcomp(elType,'1st') || strcomp(elType,'1')
      NxTri = Nx/2;
      NyTri = Ny/2;
      filteredNodes = zeros((2*sum(NxTri)+1)*sum(NyTri)+sum(NxTri)+1,2);
      edges = zeros((6*sum(NxTri)+1)*sum(NyTri)+sum(NxTri),2);
      elements = zeros(sum(Nx)*sum(Ny),3);
      for j=1:sum(NyTri)
        filteredNodes((j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1,1:2) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+1:2:(j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv)+1,1:2);
        filteredNodes((j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1),1:2) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv)+1+2:2:j*2*(sum(NxEquiv)+1)-1,1:2);
        elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
        elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
        elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),3) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
        elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
        elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),2) = j*(2*sum(NxTri)+1)+2:j*(2*sum(NxTri)+1)+sum(NxTri)+1;
        elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),3) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
        elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
        elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),2) = j*(2*sum(NxTri)+1)+2:j*(2*sum(NxTri)+1)+sum(NxTri)+1;
        elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),3) = j*(2*sum(NxTri)+1)+1:j*(2*sum(NxTri)+1)+sum(NxTri);
        elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
        elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
        elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),3) = j*(2*sum(NxTri)+1)+1:j*(2*sum(NxTri)+1)+sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
        edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
        edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
        edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
        edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
        edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),2) = j*(2*sum(NxTri)+1)+2:j*(2*sum(NxTri)+1)+sum(NxTri)+1;
        edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),1) = j*(2*sum(NxTri)+1)+1:j*(2*sum(NxTri)+1)+sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
        edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),2) = j*(2*sum(NxTri)+1)+1:j*(2*sum(NxTri)+1)+sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
        edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,2) = j*(2*sum(NxTri)+1)+sum(NxTri)+1;
      end
      j = sum(NyTri)+1;
      filteredNodes((j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1,1:2) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+1:2:(j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv)+1,1:2);
      edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
      edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
      nodes = filteredNodes;
    elseif strcomp(elType,'second') || strcomp(elType,'Second') || strcomp(elType,'2nd') || strcomp(elType,'2')
      NxTri = Nx/2;
      NyTri = Ny/2;
      filteredNodes = zeros((8*sum(NxTri)+2)*sum(NyTri)+2*sum(NxTri)+1,2);
      edges = zeros((6*sum(NxTri)+1)*sum(NyTri)+sum(NxTri),3);
      elements = zeros(sum(Nx)*sum(Ny),6);
      for j=1:sum(NyTri)
        filteredNodes((j-1)*(8*sum(NxTri)+2)+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1,1:2) = baseNodes((j-1)*4*sum(NxEquiv)+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1,1:2);
        filteredNodes((j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri),1:2) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+2:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv),1:2);
        filteredNodes((j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1,1:2) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1,1:2);
        filteredNodes((j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2*sum(NxTri),1:2) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1+2:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv),1:2);
        elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
        elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
        elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
        elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri);
        elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),5) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri);
        elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)-1;
        elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
        elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),2) = j*(8*sum(NxTri)+2)+3:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
        elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
        elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1;
        elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),5) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri);
        elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri);
        elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
        elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),2) = j*(8*sum(NxTri)+2)+3:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
        elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),3) = j*(8*sum(NxTri)+2)+1:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
        elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri);
        elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),5) = j*(8*sum(NxTri)+2)+2:2:j*(8*sum(NxTri)+2)+2*sum(NxTri);
        elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri)-1;
        elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
        elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
        elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),3) = j*(8*sum(NxTri)+2)+1:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
        elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)-1;
        elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),5) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri)-1;
        elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)-1;
        edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
        edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
        edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
        edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)-1;
        edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
        edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),3) = j*(8*sum(NxTri)+2)+3:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
        edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),1) = j*(8*sum(NxTri)+2)+1:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
        edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri)-1;
        edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
        edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)-1;
        edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),3) = j*(8*sum(NxTri)+2)+1:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
        edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
        edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1;
        edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,3) = j*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
      end
      j = sum(NyTri)+1;
      filteredNodes((j-1)*(8*sum(NxTri)+2)+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1,1:2) = baseNodes((j-1)*4*sum(NxEquiv)+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1,1:2);
      edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
      edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri);
      edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
      nodes = filteredNodes;
    end
  end
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier,'\n'])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end
writeToLogFile(logfullfile,['... done.','\n'])

elapsed = toc(start);
writeToLogFile(logfullfile,'Timer stopped.\n')
writeToLogFile(logfullfile,['\nELAPSED WALLCLOCK TIME: ', num2str(elapsed),' [s]\n\n'])
writeToLogFile(logfullfile,'Exiting function: filterRectangularMesh\n')

return
