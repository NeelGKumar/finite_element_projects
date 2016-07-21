classdef QuadMesher < Brep2D
    % QuadMesher
    % Version 2013.01
    % Version 2013.02
    %   - modified the path to quadmesher
    %   - added shape argument in initializeMesher
    % Version 2013.03
    % - removed area computation and myArea (added to brep)
    % Version 2013.04
    % - bug in quadmesher line 146
    % - Added myBoundaryNodes
    % - mChanged setPoints to setMeshPoints
    % Version 2013.05
    % - Line 146: bug in distance to mesher
    properties(GetAccess = 'public', SetAccess = 'private')
        myShape;
        myModelScale;
        myMesh;
        myElemSize;
        myNodesPerEdge;
        myNodesPerElement;
        myNumNodes;
        myNumElemsDesired;
        myNumElems;
        myPseudoDensity; % for topology optimization
        myNumBoundarySegments;
        myBoundaryNodes;
    end
    methods
        function obj = QuadMesher(brep,nElements,shape) 
            % brep can be a brep structure or a file name
            obj = obj@Brep2D(brep); % call superclass 
            obj.myShape = 'Linear';% only linear shape is supported
            obj.myNumElemsDesired = nElements;
            obj = obj.initializeMesher(nElements,shape);
  
        end
        function obj = initializeMesher(obj,nElements,~)
            % shape argument neglected
            % Creates a mesh using an external quad mesher
            % Only linear shape functions are supported.
            obj.myNumBoundarySegments = size(obj.myBrep.segments,2);
            x = obj.myBrep.vertices(1,:);
            y = obj.myBrep.vertices(1,:);
            obj.myModelScale = max(max(x)-min(x),max(y)-min(y));
            obj.myElemSize = sqrt(obj.myArea/nElements);
            h = obj.myElemSize;
            quadMesher = which('quadmesher.exe') ;
            quadMeshDir = fileparts(quadMesher);
            if (isempty(quadMeshDir))
                disp('Could not find quadmesher.exe; please add to matlab path');
                return
            end
                
            saveDir = pwd;
            cd(quadMeshDir);
            
            % Now write the data 
            inputFileName = 'tempquadmesh.inp'; % don't change, this is used in quadmesher.par
            dos(['del ' inputFileName ' &']);
            fid = fopen(inputFileName','w');
            fprintf(fid,'1\n');
            fprintf(fid,'%d\n',obj.myNumBoundarySegments);
            fprintf(fid,'%d\n',obj.myNumBoundarySegments);
            brep = obj.myBrep;
            v = brep.vertices;
            for i = 1:obj.myNumBoundarySegments
                breptype = brep.segments(1,i);
                vs = brep.segments(2,i);
                ve = brep.segments(3,i);
                xs = v(1,vs);xe = v(1,ve);
                ys = v(2,vs);ye = v(2,ve);
                if (i == 1)
                    fprintf(fid,'%f %f %f\n',xs,ys,h);
                end
                if ((breptype == 1) || (breptype == -1)) % line or line-connector between loops
                    fprintf(fid,'%f %f %f\n',xe,ye,h);
                    fprintf(fid,'lin\n');
                elseif (breptype == 2) % arc 
                    fprintf(fid,'%f %f %f\n',xe,ye,h);
                    vc = abs(brep.segments(4,i));
                    xc = v(1,vc);yc = v(2,vc);
                    R = sqrt((xc-xs)^2+(yc-ys)^2);
                    fprintf(fid,'arc\n');
                    if (brep.segments(4,i) > 0)
                        fprintf(fid,'%f \n',-R);
                    else
                        fprintf(fid,'%f \n',R);
                    end
                end
            end
            fprintf(fid,'lin');
            fclose(fid);
            outputFileName =  'tempquadmesh.mesh'; %don't change, used in quadmesher.par
            system(['del ' outputFileName]);% remove previous copy
            %disp('Attempting to create a quad mesh ... ');
            cmd ='quadmesher.exe quadmesher.par /NOCONSOLE  /f >nul 2>&1 &';
            system(cmd);
            while (exist(outputFileName,'file')==0)
                % wait for file
            end
            %disp('Meshing complete ... ');
            eval('!taskkill /t /im cmd.exe /f >nul 2>&1');
            % now read the quadmesh file     
            fid = fopen(outputFileName,'r');
            fgetl(fid); % read and discard first line
            data = fscanf(fid,' %d \n',5);
            obj.myNodesPerElement = 4;
            obj.myNumNodes = data(1);
            obj.myNumElems = data(2);
            obj.myMesh.p = zeros(2,obj.myNumNodes);
            obj.myMesh.q = zeros(4,obj.myNumElems);
            data = fscanf(fid,' %f ',[4,obj.myNumNodes]);
            obj.myMesh.p(:,data(1,:)) = data(2:3,:);
            %Now read items as below
            %1 0 quad        1        2        4        3
            %2 0 quad        2        5        6        4

            C = textscan(fid,'%d%d%s%d%d%d%d',obj.myNumElems);
            quads = C{1};
            obj.myMesh.q(1,quads)  = C{4};
            obj.myMesh.q(2,quads)  = C{5};
            obj.myMesh.q(3,quads)  = C{6};
            obj.myMesh.q(4,quads)  = C{7};

            fclose(fid);
            cd(saveDir);
            % now figure out the number of edges
            % gather all pairs of nodes
            mesh = obj.myMesh;
            pairs = [mesh.q([1,2],:), mesh.q([2,3],:), mesh.q([3,4],:),mesh.q([4,1],:)];
            toFlip = pairs(1,:) > pairs(2,:);% make sure the first node is smaller
            pairs(:,toFlip) = flipud(pairs(:,toFlip));
            
            [~,IF] = unique(pairs','rows','first'); % find the unique pairs of nodes
            [~,IL] = unique(pairs','rows','last'); % find the unique pairs of nodes
            index = intersect(IF,IL);
            pairs = pairs(:,index);
            
            
            % the unique pairs correspond to edges
            % We need to find the boundary segment the edges fall on
            nMeshEdges = size(pairs,2);
            xStart = mesh.p(1,pairs(1,:));yStart = mesh.p(2,pairs(1,:));  
            xEnd = mesh.p(1,pairs(2,:));yEnd = mesh.p(2,pairs(2,:));  
            pts(1,:) = (xStart + xEnd)/2;  pts(2,:) = (yStart + yEnd)/2; 
            % find the distance to each segment
            dist = 1e12*ones(obj.myNumBoundarySegments,nMeshEdges);
            brep = obj.myBrep;
            for seg = 1:obj.myNumBoundarySegments
                breptype = brep.segments(1,seg);
                vs = brep.segments(2,seg);
                ve = brep.segments(3,seg);
                if (breptype == 1) % line
                    dist(seg,:) = obj.distOfPointsToLineSegment(pts,v(:,vs),v(:,ve));
                elseif (breptype == 2) % arc
                    disp('Distance to arc not implemented');
                elseif (breptype == 3) % arc
                    disp('Not implemented');
                end
            end
            obj.myNodesPerEdge = 2;
            [~,bndryIndex] = min(dist);
            obj.myMesh.e = zeros(5,nMeshEdges);
            obj.myMesh.e(1:2,:) = pairs;
            obj.myMesh.e(5,:) = bndryIndex;
            obj.myPseudoDensity = ones(obj.myNumElems,1);
            
            bndryNodes = obj.myMesh.e(1:2,:);
            bndryNodes = bndryNodes(:);
            obj.myBoundaryNodes = unique(bndryNodes);
        end
        
        function obj = resetBrepAndMesh(obj,brep)
            % Resets a brep (useful during shape optimization) and
            % recreates the mesh.
            obj = obj.resetBrep(brep);
            obj = obj.initializeMesher(obj.myNumElemsDesired,obj.myShape);
        end
        
        function obj = setPseudoDensity(obj,elems,value)
            % Pseudo-density is useful for turning off elements; useful
            % during topology optimization
            obj.myPseudoDensity(:) = 1;
            obj.myPseudoDensity(elems) = value;
        end
         
        function obj = setMeshPoints(obj,nodes,p)
            obj.myMesh.p(:,nodes) = p;
        end
        function plotMesh(obj)
            X = zeros(obj.myNumElems,5);
            Y = zeros(obj.myNumElems,5);
            for count = 1:obj.myNumElems
                if (obj.myPseudoDensity(count) == 0), continue;end
                nodes = obj.myMesh.q(:,count);
                X(count,:) = obj.myMesh.p(1,[nodes' nodes(1)]);
                Y(count,:) = obj.myMesh.p(2,[nodes' nodes(1)]);
            end
            plot( X',Y','b','LineWidth',0.8); hold on;% quad
            axis equal; axis on;view(2);
            obj.adjustFigScale();
            hold on;
        end
        function [N,gradN] = edgeShapeFunction(obj,xi)
            if (strcmp(obj.myShape,'Linear'))
                N = [(1-xi)/2;
                    (1+xi)/2];
                gradN = [-1/2;
                    1/2];
            elseif strcmp(obj.myShape,'Quadratic')
                N = [xi*(xi-1)/2;
                    (1-xi)*(1+xi);
                    xi*(xi+1)/2];
                gradN = [(2*xi-1)/2;
                    (-2*xi);
                    (2*xi+1)/2];
            end
        end
        function  [N,gradN] = QuadShapeFunction(obj,xi,eta)
            if (strcmp(obj.myShape,'Linear')) % linear 
                N = 0.25*[(1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];
                gradN = 0.25*[eta-1 1-eta eta+1 -eta-1; xi-1 -xi-1 xi+1 1-xi];
            elseif (strcmp(obj.myShape,'Quadratic')) % quadratic 
               disp('Not implemented');
            end
        end
        function [J] = Jacobian(obj,elem,xi,eta)
            nodes = obj.myMesh.q(1:obj.myNodesPerElement,elem)';
            xNodes = obj.myMesh.p(1,nodes);
            yNodes = obj.myMesh.p(2,nodes);
            [~,gradN] = QuadShapeFunction(obj,xi,eta);
            J = zeros(2,2);
            J(1,1) = gradN(1,:)*xNodes';
            J(1,2) = gradN(2,:)*xNodes';
            J(2,1) = gradN(1,:)*yNodes';
            J(2,2) = gradN(2,:)*yNodes';
        end
    end
     methods(Static)
        function [xi_GQ,eta_GQ, wt_GQ] = GaussQuad()
            xi_GQ = [-1/sqrt(3) 1/sqrt(3)  1/sqrt(3) -1/sqrt(3)];
            eta_GQ = [-1/sqrt(3) -1/sqrt(3)  1/sqrt(3) 1/sqrt(3)];
            wt_GQ = [1 1 1 1]; 
        end
        function [xi_GQ, wt_GQ] = GaussQLine(nPoints)
            % Gauss quadrature pts for a line (-1 to 1)
            if (nPoints <= 1)
                xi_GQ = 0.0;
                wt_GQ = 2;
            elseif (nPoints == 2)
                xi_GQ = [-0.577350269189626 0.577350269189626];
                wt_GQ = [1 1];
            elseif (nPoints == 3)
                xi_GQ = [-0.774596669241483 0 0.774596669241483];
                wt_GQ = [ 0.555555555555556 0.88888888888888889 0.5555555555555556];
            elseif (nPoints == 4)
                xi_GQ = [-0.8611363115 -0.3399810435 0.3399810435 0.8611363115];
                wt_GQ = [0.3478548451 0.6521451548 0.6521451548 0.3478548451];
            elseif (nPoints >= 5)
                xi_GQ = [-0.906179846 -0.53846931 0  0.53846931 0.906179846];
                wt_GQ = [0.236926885 0.4786286704 0.5688888888 0.4786286704 0.236926885];
            end
        end
    end
end
