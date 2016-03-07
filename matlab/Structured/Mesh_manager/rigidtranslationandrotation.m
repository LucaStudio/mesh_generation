function[mesh]=rigidtranslationandrotation(D,oldmesh,r0,T)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: March 27th, 2014
%    Last update: April 1st, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

mesh = oldmesh;

switch D
    case 2
        if oldmesh.edgeflag && oldmesh.faceflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = r0(1) + T(1,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i)];
                mesh.nodes.meshcoordinates(2,i) = r0(2) + T(2,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i)];
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = r0(1) + T(1,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i)];
                mesh.edges.centroid(2,i) = r0(2) + T(2,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i)];
            end
            for i=1:oldmesh.totF
                mesh.faces.centroid(1,i) = r0(1) + T(1,:)*[oldmesh.faces.centroid(1,i);oldmesh.faces.centroid(2,i)];
                mesh.faces.centroid(2,i) = r0(2) + T(2,:)*[oldmesh.faces.centroid(1,i);oldmesh.faces.centroid(2,i)];
            end
        elseif oldmesh.edgeflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = r0(1) + T(1,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i)];
                mesh.nodes.meshcoordinates(2,i) = r0(2) + T(2,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i)];
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = r0(1) + T(1,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i)];
                mesh.edges.centroid(2,i) = r0(2) + T(2,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i)];
            end
        else
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = r0(1) + T(1,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i)];
                mesh.nodes.meshcoordinates(2,i) = r0(2) + T(2,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i)];
            end
        end
    case 3
        if oldmesh.edgeflag && oldmesh.faceflag && oldmesh.cellflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = r0(1) + T(1,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i);oldmesh.nodes.meshcoordinates(3,i)];
                mesh.nodes.meshcoordinates(2,i) = r0(2) + T(2,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i);oldmesh.nodes.meshcoordinates(3,i)];
                mesh.nodes.meshcoordinates(3,i) = r0(3) + T(3,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i);oldmesh.nodes.meshcoordinates(3,i)];
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = r0(1) + T(1,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i);oldmesh.edges.centroid(3,i)];
                mesh.edges.centroid(2,i) = r0(2) + T(2,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i);oldmesh.edges.centroid(3,i)];
                mesh.edges.centroid(3,i) = r0(3) + T(3,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i);oldmesh.edges.centroid(3,i)];
            end
            for i=1:oldmesh.totF
                mesh.faces.centroid(1,i) = r0(1) + T(1,:)*[oldmesh.faces.centroid(1,i);oldmesh.faces.centroid(2,i);oldmesh.faces.centroid(3,i)];
                mesh.faces.centroid(2,i) = r0(2) + T(2,:)*[oldmesh.faces.centroid(1,i);oldmesh.faces.centroid(2,i);oldmesh.faces.centroid(3,i)];
                mesh.faces.centroid(3,i) = r0(3) + T(3,:)*[oldmesh.faces.centroid(1,i);oldmesh.faces.centroid(2,i);oldmesh.faces.centroid(3,i)];
            end
            for i=1:oldmesh.totC
                mesh.cells.centroid(1,i) = r0(1) + T(1,:)*[oldmesh.cells.centroid(1,i);oldmesh.cells.centroid(2,i);oldmesh.cells.centroid(3,i)];
                mesh.cells.centroid(2,i) = r0(2) + T(2,:)*[oldmesh.cells.centroid(1,i);oldmesh.cells.centroid(2,i);oldmesh.cells.centroid(3,i)];
                mesh.cells.centroid(3,i) = r0(3) + T(3,:)*[oldmesh.cells.centroid(1,i);oldmesh.cells.centroid(2,i);oldmesh.cells.centroid(3,i)];
            end
        elseif oldmesh.edgeflag && oldmesh.faceflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = r0(1) + T(1,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i);oldmesh.nodes.meshcoordinates(3,i)];
                mesh.nodes.meshcoordinates(2,i) = r0(2) + T(2,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i);oldmesh.nodes.meshcoordinates(3,i)];
                mesh.nodes.meshcoordinates(3,i) = r0(3) + T(3,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i);oldmesh.nodes.meshcoordinates(3,i)];
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = r0(1) + T(1,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i);oldmesh.edges.centroid(3,i)];
                mesh.edges.centroid(2,i) = r0(2) + T(2,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i);oldmesh.edges.centroid(3,i)];
                mesh.edges.centroid(3,i) = r0(3) + T(3,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i);oldmesh.edges.centroid(3,i)];
            end
            for i=1:oldmesh.totF
                mesh.faces.centroid(1,i) = r0(1) + T(1,:)*[oldmesh.faces.centroid(1,i);oldmesh.faces.centroid(2,i);oldmesh.faces.centroid(3,i)];
                mesh.faces.centroid(2,i) = r0(2) + T(2,:)*[oldmesh.faces.centroid(1,i);oldmesh.faces.centroid(2,i);oldmesh.faces.centroid(3,i)];
                mesh.faces.centroid(3,i) = r0(3) + T(3,:)*[oldmesh.faces.centroid(1,i);oldmesh.faces.centroid(2,i);oldmesh.faces.centroid(3,i)];
            end
        elseif oldmesh.edgeflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = r0(1) + T(1,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i);oldmesh.nodes.meshcoordinates(3,i)];
                mesh.nodes.meshcoordinates(2,i) = r0(2) + T(2,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i);oldmesh.nodes.meshcoordinates(3,i)];
                mesh.nodes.meshcoordinates(3,i) = r0(3) + T(3,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i);oldmesh.nodes.meshcoordinates(3,i)];
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = r0(1) + T(1,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i);oldmesh.edges.centroid(3,i)];
                mesh.edges.centroid(2,i) = r0(2) + T(2,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i);oldmesh.edges.centroid(3,i)];
                mesh.edges.centroid(3,i) = r0(3) + T(3,:)*[oldmesh.edges.centroid(1,i);oldmesh.edges.centroid(2,i);oldmesh.edges.centroid(3,i)];
            end
        else
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = r0(1) + T(1,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i);oldmesh.nodes.meshcoordinates(3,i)];
                mesh.nodes.meshcoordinates(2,i) = r0(2) + T(2,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i);oldmesh.nodes.meshcoordinates(3,i)];
                mesh.nodes.meshcoordinates(3,i) = r0(3) + T(3,:)*[oldmesh.nodes.meshcoordinates(1,i);oldmesh.nodes.meshcoordinates(2,i);oldmesh.nodes.meshcoordinates(3,i)];
            end
        end
    otherwise
        disp('Dimension of space requested is currently not supported')
end

return