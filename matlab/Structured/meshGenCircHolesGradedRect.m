function[nodes,elements,edges]=meshGenCircHolesGradedRect(logfullfile,elType,elOrder,x0,y0,lx,ly,Nx,Ny,holes)
%%
%==============================================================================
% Copyright (c) 2016 - 2017 Université de Lorraine & Luleå tekniska universitet
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
% Neither the name of the Université de Lorraine & Luleå tekniska universitet
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
%  A function to mesh a non-simply connected, i.e. with rectangular holes,
%  rectangular geometry with rectangular elements
%
%  Input: x0 - scalar - x-coordinate of center
%         y0 - scalar - y-coordinate of center
%         lx - [M x 1] vector - Side length in each one of the M mesh regions in x-direction
%         ly - [N x 1] vector - Side length in each one of the N mesh regions in y-direction
%         Nx - [M x 1] vector - Number of ELEMENTS in each one of the M mesh regions in x-direction
%         Ny - [N x 1] vector - Number of ELEMENTS in each one of the N mesh regions in y-direction
%         holes - [H x 3] matrix - H is the number of holes; for each hole the following data must pe provided:
%                                  xC - scalar - x-coordinate of hole's center
%                                  yC - scalar - y-coordinate of hole's center
%                                  R - scalar - Radius
%%

writeToLogFile(logfullfile,'In function: meshGenCircHolesGradedRect\n')
writeToLogFile(logfullfile,'\nStarting timer\n')
start = tic;

writeToLogFile(logfullfile,['Creating base rectangular mesh ...','\n'])
try
  if strcomp(elType,'quads') || strcomp(elType,'quad') || strcomp(elType,'quadrilaterals') || strcomp(elType,'quadrilateral')
    if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
      writeToLogFile(logfullfile,['    Type ad order of elements : First order quadrilaterals','\n'])
      NxEquiv = Nx;
      NyEquiv = Ny;
    elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
      writeToLogFile(logfullfile,['    Type ad order of elements : Second order quadrilaterals','\n'])
      NxEquiv = 2*Nx;
      NyEquiv = 2*Ny;
    end
  elseif strcomp(elType,'tris') || strcomp(elType,'tri') || strcomp(elType,'triangles') || strcomp(elType,'triangle')
    if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
      writeToLogFile(logfullfile,['    Type ad order of elements : First order triangles','\n'])
      NxEquiv = Nx;
      NyEquiv = Ny;
    elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
      writeToLogFile(logfullfile,['    Type ad order of elements : Second order triangles','\n'])
      NxEquiv = 2*Nx;
      NyEquiv = 2*Ny;
    end
  end
  writeToLogFile(logfullfile,['    Calling function ', 'gradedRectangle',' ...\n']);
  baseMesh = gradedRectangle(logfullfile,x0,y0,lx,ly,NxEquiv,NyEquiv);
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['... done.','\n'])
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end

writeToLogFile(logfullfile,['Create holes and deform internal boundaries ...','\n'])
try
  baseMesh = [baseMesh (1:length(baseMesh))' zeros(length(baseMesh),1) ones(length(baseMesh),1) 2*ones(length(baseMesh),1) zeros(length(baseMesh),1)]; % x y <node index in full mesh> <node index in reduced mesh> isInBulk isInHole/isOnBoundary/isInBulk <internal boundary topology>
  defMesh = baseMesh;
  DeltaX = sum(lx)/sum(Nx);
  DeltaY = sum(ly)/sum(Ny);
  x1 = x0 - 0.5*sum(lx) + 10*DeltaX;
  x2 = x0 + 0.5*sum(lx) - 10*DeltaX;
  y1 = y0 - 0.5*sum(ly) + 10*DeltaY;
  y2 = y0 + 0.5*sum(ly) - 10*DeltaY;
  for k=1:length(holes)
    xC = holes(k,1);
    yC = holes(k,2);
    R = holes(k,3);
    if x1>xC+R || x2<xC+R || x1>xC-R || x2<xC-R
      writeToLogFile(logfullfile,['Hole number ',num2str(k),' is at least partially outside the domain. ','\n'])
      writeToLogFile(logfullfile,['Skipping it. ','\n'])
      continue
    end
    A.center = [xC yC];
    A.radius = R;
    for h=1:k-1
      B.center = [holes(h,1) holes(h,2)];
      B.radius = holes(h,3);
      if checkCircleIntersections(A,B)>0
        writeToLogFile(logfullfile,['Hole number ',num2str(k),' intersects a preceding one. ','\n'])
        writeToLogFile(logfullfile,['Skipping it. ','\n'])
        continue
      end
    end
    L = R/sqrt(2);
    xH1 = x0 - L;
    xH2 = x0 + L;
    yH1 = y0 - L;
    yH2 = y0 + L;
    for j=2:sum(Ny)
      for i=2:sum(Nx)
        xP = baseMesh((j-1)*(sum(Nx)+1)+i,1);    % current point
        yP = baseMesh((j-1)*(sum(Nx)+1)+i,2);
        xRi = baseMesh((j-1)*(sum(Nx)+1)+i+1,1); % neighbour at the right
        yRi = baseMesh((j-1)*(sum(Nx)+1)+i+1,2);
        xLe = baseMesh((j-1)*(sum(Nx)+1)+i-1,1); % neighbour at the left
        yLe = baseMesh((j-1)*(sum(Nx)+1)+i-1,2);
        xUp = baseMesh(j*(sum(Nx)+1)+i,1);       % upper neighbour
        yUp = baseMesh(j*(sum(Nx)+1)+i,2);
        xLo = baseMesh((j-2)*(sum(Nx)+1)+i,1);   % lower neighbour
        yLo = baseMesh((j-2)*(sum(Nx)+1)+i,2);
        xLoLe = baseMesh((j-2)*(sum(Nx)+1)+i-1,1); % lower left neighbour
        yLoLe = baseMesh((j-2)*(sum(Nx)+1)+i-1,2);
        xLoRi = baseMesh((j-2)*(sum(Nx)+1)+i+1,1); % lower right neighbour
        yLoRi = baseMesh((j-2)*(sum(Nx)+1)+i+1,2);
        xUpRi = baseMesh(j*(sum(Nx)+1)+i+1,1); % upper right neighbour
        yUpRi = baseMesh(j*(sum(Nx)+1)+i+1,2);
        xUpLe = baseMesh(j*(sum(Nx)+1)+i-1,1); % upper left neighbour
        yUpLe = baseMesh(j*(sum(Nx)+1)+i-1,2);
        if xH1<=xP && xH2>=xP && (yH1==yP || yH2==yP)     % either on the south or north hole's boundary line
          baseMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
          defMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
          defMesh((j-1)*(sum(Nx)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
          defMesh((j-1)*(sum(Nx)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
          if yH1==yP
            baseMesh((j-1)*(sum(Nx)+1)+i,5) = 7;
            defMesh((j-1)*(sum(Nx)+1)+i,5) = 7;
          else
            baseMesh((j-1)*(sum(Nx)+1)+i,5) = 5;
            defMesh((j-1)*(sum(Nx)+1)+i,5) = 5;
          end
        elseif (xH1==xP || xH2==xP) && yH1<=yP && yH2>=yP % either on the left or right hole's boundary line
          baseMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
          defMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
          defMesh((j-1)*(sum(Nx)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
          defMesh((j-1)*(sum(Nx)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
          if xH1==xP
            baseMesh((j-1)*(sum(Nx)+1)+i,5) = 6;
            defMesh((j-1)*(sum(Nx)+1)+i,5) = 6;
          else
            baseMesh((j-1)*(sum(Nx)+1)+i,5) = 8;
            defMesh((j-1)*(sum(Nx)+1)+i,5) = 8;
          end
        elseif xH1<xP && xH2>xP && yH1<yP && yH2>yP     % inside hole
          baseMesh((j-1)*(sum(Nx)+1)+i,3) = 0;
          defMesh((j-1)*(sum(Nx)+1)+i,3) = 0;
          baseMesh((j-1)*(sum(Nx)+1)+i,4) = 0;
          defMesh((j-1)*(sum(Nx)+1)+i,4) = 0;
        else                                            % outside hole
          if xH1<xLoLe && xH2>xLoLe && yH1<yLoLe && yH2>yLoLe       % lower left neighbour inside hole
            baseMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            baseMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(Nx)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(Nx)+1)+i,5) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,5) = 1;
          elseif xH1<xLoRi && xH2>xLoRi && yH1<yLoRi && yH2>yLoRi   % lower right neighbour inside hole
            baseMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            baseMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(Nx)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(Nx)+1)+i,5) = 2;
            defMesh((j-1)*(sum(Nx)+1)+i,5) = 2;
          elseif xH1<xUpRi && xH2>xUpRi && yH1<yUpRi && yH2>yUpRi   % upper right neighbour inside hole
            baseMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            baseMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(Nx)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(Nx)+1)+i,5) = 3;
            defMesh((j-1)*(sum(Nx)+1)+i,5) = 3;
          elseif xH1<xUpLe && xH2>xUpLe && yH1<yUpLe && yH2>yUpLe   % upper left neighbour inside hole
            baseMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            baseMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(Nx)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(Nx)+1)+i,5) = 4;
            defMesh((j-1)*(sum(Nx)+1)+i,5) = 4;
          elseif xH1<xRi && xH2>xRi && yH1<yRi && yH2>yRi   % right neighbour inside hole
            baseMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            baseMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(Nx)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(Nx)+1)+i,5) = 6;
            defMesh((j-1)*(sum(Nx)+1)+i,5) = 6;
          elseif xH1<xLe && xH2>xLe && yH1<yLe && yH2>yLe   % left neighbour inside hole
            baseMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            baseMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(Nx)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(Nx)+1)+i,5) = 8;
            defMesh((j-1)*(sum(Nx)+1)+i,5) = 8;
          elseif xH1<xUp && xH2>xUp && yH1<yUp && yH2>yUp   % upper neighbour inside hole
            baseMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            baseMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(Nx)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(Nx)+1)+i,5) = 7;
            defMesh((j-1)*(sum(Nx)+1)+i,5) = 7;
          elseif xH1<xLo && xH2>xLo && yH1<yLo && yH2>yLo   % lower neighbour inside hole
            baseMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            baseMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,4) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(Nx)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(Nx)+1)+i,2)-yC),(baseMesh((j-1)*(sum(Nx)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(Nx)+1)+i,5) = 5;
            defMesh((j-1)*(sum(Nx)+1)+i,5) = 5;
          else                                              % bulk
            baseMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            defMesh((j-1)*(sum(Nx)+1)+i,3) = 1;
            baseMesh((j-1)*(sum(Nx)+1)+i,4) = 2;
            defMesh((j-1)*(sum(Nx)+1)+i,4) = 2;
          end
        end
      end
    end
  end
  writeToLogFile(logfullfile,['... done.','\n'])
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end

writeToLogFile(logfullfile,['Deforming the mesh ...','\n'])
try
  writeToLogFile(logfullfile,['    Calling function ', 'getindices2D',' ...\n']);
  [indicesbulk,indicesinternalbulk,...
   indicesE1,indicesE2,indicesE3,indicesE4,...
   indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,...
   indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,...
   indicesC1,indicesC2,indicesC3,indicesC4,...
   indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(sum(Nx)+1,sum(Ny)+1);
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['    Calling function ', 'build_neighbourhoods2D',' ...\n']);
  [temp1,temp2,temp3,firstdevneighbours] = build_neighbourhoods2D((sum(Nx)+1)*(sum(Ny)+1),sum(Nx)+1,...
                                                                  0,0,...
                                                                  1,0,0,...
                                                                  indicesbulk,indicesinternalbulk,...
                                                                  indicesE1,indicesE2,indicesE3,indicesE4,...
                                                                  indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,...
                                                                  indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,...
                                                                  indicesC1,indicesC2,indicesC3,indicesC4,...
                                                                  indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['    Calling function ', 'sparseellipticgridgen2Dinternalboundaries',' ...\n']);
  sparseellipticgridgen2Dinternalboundaries(sum(Nx)+1,(sum(Nx)+1)*(sum(Ny)+1),lattice,deltaq,boundary,flagperiodicity,periodicity,...
                                            indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,...
                                            indicesC1,indicesC2,indicesC3,indicesC4,...
                                            firstdevneighbours,itmax,tol,0)

  writeToLogFile(logfullfile,['... done.','\n'])
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end

elapsed = toc(start);
writeToLogFile(logfullfile,'Timer stopped.\n')
writeToLogFile(logfullfile,['\nELAPSED WALLCLOCK TIME: ', num2str(elapsed),' [s]\n\n'])
writeToLogFile(logfullfile,'Exiting function: meshGenCircHolesGradedRect\n')

return
