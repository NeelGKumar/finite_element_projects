classdef plateFEMindlin < QuadMesher
    properties(GetAccess = 'public', SetAccess = 'private')
        myE; % Young's modulus
        myGyz; %shear stiffness Gyz
        myGxz; %shear stiffness Gxz
        myNu; % Poisson's ratio
        myD; % Elasticity matrix Dm
        myK; % stiffness matrix;
        myF; % force vector
        myU; % Solution vector
        myH; %element thickness
        myElementAreas; % area of each element
        myTotalArea; % total area
        myFixedDOF; % fixed degrees of freedom
        myGaussPts_4; % 4 Gauss pts for integrating quad elements
        myGaussWts_4; % Gauss wts
        myGaussPts_9; % 9 Gauss pts for integrating quad elements
        myGaussWts_9; % Gauss wts
        myTotalLinearSolves
        myDensity %material density
        myResults
        myDOFperNode %5 for 8 Serendipity nodes, 4 for center node
        myBCtype
        myBCvalue
        myNumDOF %globalDOF
        myBodyForce %[Fx Fy Fw Mx My] set all zero except Fw
    end
    methods
        function obj= plateFEMindlin(brepFileName,nElements,shape)
            obj = obj@QuadMesher(brepFileName,nElements,shape); % call superclass
            obj.myE = 2e11;
            obj.myGyz = 1e9; %shear modulus
            obj.myGxz = 1e9; %shear modulus
            obj.myNu = 0.33;
            obj.myH = 0.1; %element thickness
            obj.myDensity=7800; %steel
            obj.myBodyForce=[0;0;1;0;0];%[Fx Fy Fw Mx My]
            E = obj.myE;
            Gyz = obj.myGyz;
            Gxz = obj.myGxz;
            nu = obj.myNu;
            h = obj.myH;
            %************Kirchoff plate bending elasticity matrix Dk
            d = (E*h^3)/(12*(1-nu^2));
            Dk = [d nu*d 0;nu*d d 0;0 0 (1-nu)*d/2]; %3X3
            %************Plane stress elasticity matrix Ematrix 
            Ematrix = 1/E*[1 -nu 0;-nu 1 0;0 0 2*(1+nu)]; %3X3
            Dm = zeros(8,8);
            Dm(1:3,1:3) = Ematrix*obj.myH;
            Dm(4:6,4:6) = Dk;
            Dm(7,7) = Gyz*obj.myH/1.2;
            Dm(8,8) = Gxz*obj.myH/1.2;
            %Elasticity matrix Dm
            obj.myD = Dm; 
            %[u,v,w,RX,RY]serendipity 8nodes;[u,v,RX,RY] center node
            obj.myDOFperNode=5;
            obj.myBCtype = zeros(obj.myNumBoundarySegments,5);% U V W RX RY
            obj.myBCvalue = zeros(obj.myNumBoundarySegments,5);
            %Selective integration rules
            obj.myGaussPts_4 =...%4 integration points
                [ -0.577350269189626 -0.577350269189626;
                0.577350269189626 -0.577350269189626;
                0.577350269189626  0.577350269189626;
                -0.577350269189626  0.577350269189626];
            obj.myGaussWts_4 = [1;1;1;1];
            obj.myGaussPts_9 =...%9 integration points
                [ -sqrt(0.6) -sqrt(0.6);
                0.00 -sqrt(0.6);
                sqrt(0.6) -sqrt(0.6);
                -sqrt(0.6) 0.00;
                0.00 0.00;
                sqrt(0.6) 0.00;
                -sqrt(0.6) sqrt(0.6);
                0.00 sqrt(0.6);
                sqrt(0.6) sqrt(0.6)];
            obj.myGaussWts_9 = [25/81;40/81;25/81;40/81;64/81;40/81;25/81;40/81;25/81];
            %DOF
            obj.myFixedDOF = []; %all fixed DOF 
            globaldof = obj.myDOFperNode*obj.myNumNodes-1*obj.myNumElems;%44 DOFs per elem
            obj.myNumDOF=globaldof;
            obj.myF = sparse(globaldof,1);
        end
        function obj = fixUOfEdge(obj,boundaryEdges)%deflection U
            obj.myBCtype(boundaryEdges,1) = 1;
        end
        function obj = fixVOfEdge(obj,boundaryEdges)%deflection V
            obj.myBCtype(boundaryEdges,2) = 1;
        end
        function obj = fixWOfEdge(obj,boundaryEdges)%deflection W
            obj.myBCtype(boundaryEdges,3) = 1;
        end
        function obj = fixRXOfEdge(obj,boundaryEdges)%rotation in X
            obj.myBCtype(boundaryEdges,4) = 1;
        end
        function obj = fixRYOfEdge(obj,boundaryEdges)%rotation in Y
            obj.myBCtype(boundaryEdges,5) = 1;
        end
        function obj = fixEdge(obj,boundaryEdges)
            obj.myBCtype(boundaryEdges,1:5) = 1;
        end
        function obj = assembleBC(obj)
            nDOF = obj.myNumDOF; %real global DOF
            DummynDOF = obj.myDOFperNode*obj.myNumNodes;%dummyNumDOF=5*Nnode
            % Gather Dirichlet boundary conditions
            isDirichlet = zeros(DummynDOF,1);
            dirValue = zeros(DummynDOF,1);
            for geomEdge = 1:size(obj.myBrep.segments,2)%num of segments
                typeU = obj.myBCtype(geomEdge,1);
                typeV = obj.myBCtype(geomEdge,2);
                typeW = obj.myBCtype(geomEdge,3);
                typeRX = obj.myBCtype(geomEdge,4);
                typeRY = obj.myBCtype(geomEdge,5);
                valueU = obj.myBCvalue(geomEdge,1);
                valueV = obj.myBCvalue(geomEdge,2);
                valueW = obj.myBCvalue(geomEdge,3);
                valueRX = obj.myBCvalue(geomEdge,4);
                valueRY = obj.myBCvalue(geomEdge,5);
                if (typeU == 1)
                    boundarySegments = find(obj.myMesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.myMesh.e(1:obj.myNodesPerEdge,seg);
                        nodes = unique(nodes(:));
                        Udof = 5*nodes-4;%dummyDOF
                        isDirichlet(Udof) = 1;
                        dirValue(Udof) =  valueU;
                    end
                end
                if (typeV == 1)
                    boundarySegments = find(obj.myMesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.myMesh.e(1:obj.myNodesPerEdge,seg);
                        nodes = unique(nodes(:));
                        Vdof = 5*nodes-3;%dummyDOF
                        isDirichlet(Vdof) = 1;
                        dirValue(Vdof) =  valueV;
                    end
                end
                if (typeW == 1)
                    boundarySegments = find(obj.myMesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.myMesh.e(1:obj.myNodesPerEdge,seg);
                        nodes = unique(nodes(:));
                        Wdof = 5*nodes-2;%dummyDOF
                        isDirichlet(Wdof) = 1;
                        dirValue(Wdof) =  valueW;
                    end
                end
                if (typeRX == 1)
                    boundarySegments = find(obj.myMesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.myMesh.e(1:obj.myNodesPerEdge,seg);
                        nodes = unique(nodes(:));
                        RXdof = 5*nodes-1;%dummyDOF
                        isDirichlet(RXdof) = 1;
                        dirValue(RXdof) =  valueRX;
                    end
                end
                if (typeRY == 1)
                    boundarySegments = find(obj.myMesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.myMesh.e(1:obj.myNodesPerEdge,seg);
                        nodes = unique(nodes(:));
                        RYdof = 5*nodes;%dummyDOF
                        isDirichlet(RYdof) = 1;
                        dirValue(RYdof) =  valueRY;
                    end
                end
            end
            %remove W DOF of center node
            dummyCenterNode = zeros(obj.myNumElems,1);
            for elem=1:obj.myNumElems
                elementnodenumbers = obj.myMesh.q(:,elem);%9 nodes
                colum=find(abs(obj.myMesh.p(1,:)-sum(obj.myMesh.p(1,elementnodenumbers))/9)<0.1);%searching X coordinate 
                for count=1:size(colum)%searching Y coordinate 
                    if (abs(obj.myMesh.p(2,colum(count))-sum(obj.myMesh.p(2,elementnodenumbers))/9)<0.1)
                        dummyCenterNode(elem)=colum(count); continue; %center node #
                    end
                end
            end
            WdummyCenterNode=5*dummyCenterNode-2;
            isDirichlet(WdummyCenterNode)=[];%remove dummy W DOF
            %fixed DOF
            dirichletDOF = find(isDirichlet(:) == 1);
            obj.myFixedDOF = dirichletDOF;
            %body force
            N=zeros(5,44); elementDof = zeros(44,1);
            bodyForce=zeros(nDOF,1);
            for elem = 1:obj.myNumElems %element body force + assembly
                elementnodenumbers = obj.myMesh.q(:,elem);%9 nodes
                dummy = [5*elementnodenumbers-4 5*elementnodenumbers-3 5*elementnodenumbers-2 ...
                    5*elementnodenumbers-1 5*elementnodenumbers];
                dummy(9,3) = [];%assume center node is NUMBERED the last one<<<<<-----need check
                %reorder to:U1 V1 w1 thetax1 thetay1,U2 V2 w2 thetax2...
                elementDof(1:5) = dummy(1,:); elementDof(6:10) = dummy(2,:);
                elementDof(11:15) = dummy(3,:); elementDof(16:20) = dummy(4,:);
                elementDof(21:25) = dummy(5,:); elementDof(26:30) = dummy(6,:);
                elementDof(31:35) = dummy(7,:); elementDof(36:40) = dummy(8,:);
                elementDof(41:44) = dummy(9,:);
                nodeCoordinates = obj.myMesh.p(:,elementnodenumbers)';
                % cycle for Gauss point
                for q = 1:size(obj.myGaussWts_4,1) %4 Gaussian integration 
                    GaussPoint = obj.myGaussPts_4(q,:);
                    s = GaussPoint(1); %epsilon
                    t = GaussPoint(2); %eta
                    %Lagrange 9 nodes' shape functions
                    N9 = (1-s^2)*(1-t^2);
                    N5 = 0.5*(1-s^2)*(1-t)-0.5*N9; N6 = 0.5*(1+s)*(1-t^2)-0.5*N9;
                    N7 = 0.5*(1-s^2)*(1+t)-0.5*N9; N8 = 0.5*(1-s)*(1-t^2)-0.5*N9;
                    N1 = 0.25*(1-s)*(1-t)-0.5*(N5+N8)-0.25*N9;
                    N2 = 0.25*(1+s)*(1-t)-0.5*(N5+N6)-0.25*N9;
                    N3 = 0.25*(1+s)*(1+t)-0.5*(N6+N7)-0.25*N9;
                    N4 = 0.25*(1-s)*(1+t)-0.5*(N7+N8)-0.25*N9;
                    N_Lagrange = [N1;N2;N3;N4;N5;N6;N7;N8;N9];
                    %Serendipity 8 nodes' shape functions
                    N85 = 0.5*(1-s^2)*(1-t); N86 = 0.5*(1+s)*(1-t^2);
                    N87 = 0.5*(1-s^2)*(1+t); N88 = 0.5*(1-s)*(1-t^2);
                    N81 = 0.25*(1-s)*(1-t)-0.5*(N85+N88);
                    N82 = 0.25*(1+s)*(1-t)-0.5*(N85+N86);
                    N83 = 0.25*(1+s)*(1+t)-0.5*(N86+N87);
                    N84 = 0.25*(1-s)*(1+t)-0.5*(N87+N88);
                    N_Serendipity = [N81;N82;N83;N84;N85;N86;N87;N88];
                    %shape function MATRIX
                    N(1,1:9) = N_Lagrange'; N(2,10:18) = N_Lagrange';
                    N(3,19:26) = N_Serendipity';
                    N(4,27:35) = N_Lagrange'; N(5,36:44) = N_Lagrange';
                    %Jacob matrix (domain transformation)
                    [obj,Ns,Nt] = derivativeFun(obj);
                    Ns=double(Ns); Nt=double(Nt); %transform symble to double 
                    Xs = Ns'*nodeCoordinates(1,:); Ys = Ns'*nodeCoordinates(2,:);
                    Xt = Nt'*nodeCoordinates(1,:); Yt = Nt'*nodeCoordinates(2,:);
                    Jacob=[Xs Ys;Xt Yt];
                    %body force
                    bodyForce(elementDof,1) =  bodyForce(elementDof,1)+...
                        obj.myBodyForce'*N*det(Jacob)*obj.myGaussWts(q); %assembly
                end
            end
            obj.myF=bodyForce;%assemble force vector
        end
        function obj = runFE(obj)
            D = obj.myD; %Elasticity matrix Dm
            obj = assembleBC(obj); %fixedDOF & assembled force vector
            globaldof = obj.myNumDOF;
            elementDof = zeros(44,1);
            obj.myK = zeros(globaldof,globaldof); %global stiffness (K) matrix
            ElemKs = zeros(globaldof,globaldof); %shear stiffness matrix
            ElemKm = zeros(globaldof,globaldof); %membrane stiffness matrix
            ElemKb = zeros(globaldof,globaldof); %bending stiffness matrix
            %*********************************************************
            for elem = 1:obj.myNumElems %Ke + assembly
                elementnodenumbers = obj.myMesh.q(:,elem);%9 nodes
                dummy = [5*elementnodenumbers-4 5*elementnodenumbers-3 5*elementnodenumbers-2 ...
                    5*elementnodenumbers-1 5*elementnodenumbers];
                dummy(9,3) = [];%assume center node is the last node
                %change order to:U1 V1 w1 thetax1 thetay1,U2 V2 w2 thetax2...
                elementDof(1:5) = dummy(1,:); elementDof(6:10) = dummy(2,:);
                elementDof(11:15) = dummy(3,:); elementDof(16:20) = dummy(4,:);
                elementDof(21:25) = dummy(5,:); elementDof(26:30) = dummy(6,:);
                elementDof(31:35) = dummy(7,:); elementDof(36:40) = dummy(8,:);
                elementDof(41:44) = dummy(9,:);
                nodeCoordinates = obj.myMesh.p(:,elementnodenumbers)';
                % cycle for Gauss point
                for count=1:3 %Ks Km Kb
                    if (count==1) %Ks
                        GaussWts=obj.myGaussWts_4;
                        GaussPts=obj.myGaussPts_4;
                    else %Km Kb
                        GaussWts=obj.myGaussWts_9;
                        GaussPts=obj.myGaussPts_9;
                    end                    
                    for q = 1:size(GaussWts,1) %Gaussian integration
                        GaussPoint = GaussPts(q,:);
                        s = GaussPoint(1); %epsilon
                        t = GaussPoint(2); %eta
                        N9 = (1-s^2)*(1-t^2);
                        N5 = 0.5*(1-s^2)*(1-t)-0.5*N9; N6 = 0.5*(1+s)*(1-t^2)-0.5*N9;
                        N7 = 0.5*(1-s^2)*(1+t)-0.5*N9; N8 = 0.5*(1-s)*(1-t^2)-0.5*N9;
                        N1 = 0.25*(1-s)*(1-t)-0.5*(N5+N8)-0.25*N9;
                        N2 = 0.25*(1+s)*(1-t)-0.5*(N5+N6)-0.25*N9;
                        N3 = 0.25*(1+s)*(1+t)-0.5*(N6+N7)-0.25*N9;
                        N4 = 0.25*(1-s)*(1+t)-0.5*(N7+N8)-0.25*N9;
                        N_Lagrange = [N1;N2;N3;N4;N5;N6;N7;N8;N9];
                        N85 = 0.5*(1-s^2)*(1-t); N86 = 0.5*(1+s)*(1-t^2);
                        N87 = 0.5*(1-s^2)*(1+t); N88 = 0.5*(1-s)*(1-t^2);
                        N81 = 0.25*(1-s)*(1-t)-0.5*(N85+N88);
                        N82 = 0.25*(1+s)*(1-t)-0.5*(N85+N86);
                        N83 = 0.25*(1+s)*(1+t)-0.5*(N86+N87);
                        N84 = 0.25*(1-s)*(1+t)-0.5*(N87+N88);
                        N_Serendipity = [N81;N82;N83;N84;N85;N86;N87;N88];
                        %domain transformation shape function
                        N(1,1:9) = N_Lagrange'; N(2,10:18) = N_Lagrange';
                        N(3,19:26) = N_Serendipity';
                        N(4,27:35) = N_Lagrange'; N(5,36:44) = N_Lagrange';
                        %Jacob 
                        [obj,Ns,Nt,N8s,N8t] = derivativeFun();
                        Ns=double(Ns); Nt=double(Nt); %transform symble to double 
                        N8s=double(N8s);N8t=double(N8t);
                        Xs = Ns'*nodeCoordinates(1,:); Ys = Ns'*nodeCoordinates(2,:);
                        Xt = Nt'*nodeCoordinates(1,:); Yt = Nt'*nodeCoordinates(2,:);
                        Jacob=[Xs Ys;Xt Yt];
                        %B matrix
                        Nix = zeros(9,1); Niy = zeros(9,1);
                        N8ix = zeros(8,1); N8iy = zeros(8,1);
                        for i=1:8
                            Nis = Ns(i); Nit = Nt(i);
                            xyDerivativeN = Jacob\[Nis;Nit];
                            Nix(i) = xyDerivativeN(1); Niy(i) = xyDerivativeN(2); %Equ C.15 
                            N8is = N8s(i); N8it = N8t(i);
                            xyDerivativeN8 = Jacob\[N8is;N8it];
                            N8ix(i) = xyDerivativeN8(1); N8iy(i) = xyDerivativeN8(2);
                        end
                        Nis = Ns(9); Nit = Nt(9);
                        xyDerivativeN = Jacob\[Nis;Nit];
                        Nix(9) = xyDerivativeN(1); Niy(9) = xyDerivativeN(2);
                        B = obj.creatBmatrix(Nix,Niy,N8ix,N8iy,N_Lagrange); %whole B matrix
                        %decompose B into Bs Bm and Bb; assemble [K]
                        Bdummy = zeros(8,44);
                        if (count==1)
                            Bdummy(end-1:end,:) = B(end-1:end,:); %shear Bs
                            %assembly Ks matrix
                            ElemKs(elementDof,elementDof) =  ElemKs(elementDof,elementDof)+...
                                Bdummy'*D*Bdummy*GaussWts(q)*det(Jacob);
                        elseif (count==2)
                            Bdummy(1:3,:) = B(1:3,:); %in-plane membrane strain Bm
                            %assembly Ks matrix
                            ElemKm(elementDof,elementDof) =  ElemKm(elementDof,elementDof)+...
                                Bdummy'*D*Bdummy*GaussWts(q)*det(Jacob);
                        elseif (count==3)
                            Bdummy(4:6,:) = B(4:6,:); %in-plane bending Bb
                            %assembly Ks matrix
                            ElemKb(elementDof,elementDof) =  ElemKb(elementDof,elementDof)+...
                                Bdummy'*D*Bdummy*GaussWts(q)*det(Jacob);
                        end
                    end
                end
            end
            obj.myK = ElemKs+ElemKm+ElemKb;%compose [K]
            obj.myK = sparse(obj.myK);
            %********************************************************
            activeDof = setdiff((1:globaldof)',(obj.myFixedDOF));
            % Use iterative method instead of direct due to hanging nodes
            obj.myU = zeros(globaldof,1);
            obj.myU(activeDof) = obj.CG(obj.myK(activeDof,activeDof),obj.myF(activeDof),1e-10,2000,obj.myU(activeDof));
            obj.myTotalLinearSolves = obj.myTotalLinearSolves + 1;
            max(obj.myU)
        end
        function [obj,Ns,Nt,N8s,N8t] = derivativeFun(obj)
            syms s t N1 N2 N3 N4 N5 N6 N7 N8 N9;
            N9 = (1-s^2)*(1-t^2);
            N5 = 0.5*(1-s^2)*(1-t)-0.5*N9; N6 = 0.5*(1+s)*(1-t^2)-0.5*N9;
            N7 = 0.5*(1-s^2)*(1+t)-0.5*N9; N8 = 0.5*(1-s)*(1-t^2)-0.5*N9;
            N1 = 0.25*(1-s)*(1-t)-0.5*(N5+N8)-0.25*N9;
            N2 = 0.25*(1+s)*(1-t)-0.5*(N5+N6)-0.25*N9;
            N3 = 0.25*(1+s)*(1+t)-0.5*(N6+N7)-0.25*N9;
            N4 = 0.25*(1-s)*(1+t)-0.5*(N7+N8)-0.25*N9;
            Ns = [diff(N1,s);diff(N2,s);diff(N3,s);diff(N4,s);diff(N5,s);diff(N6,s);...
                diff(N7,s);diff(N8,s);diff(N9,s)];
            Nt = [diff(N1,t);diff(N2,t);diff(N3,t);diff(N4,t);diff(N5,t);diff(N6,t);...
                diff(N7,t);diff(N8,t);diff(N9,t)];
            N85 = 0.5*(1-s^2)*(1-t); N86 = 0.5*(1+s)*(1-t^2);
            N87 = 0.5*(1-s^2)*(1+t); N88 = 0.5*(1-s)*(1-t^2);
            N81 = 0.25*(1-s)*(1-t)-0.5*(N85+N88);
            N82 = 0.25*(1+s)*(1-t)-0.5*(N85+N86);
            N83 = 0.25*(1+s)*(1+t)-0.5*(N86+N87);
            N84 = 0.25*(1-s)*(1+t)-0.5*(N87+N88);
            N8s = [diff(N81,s);diff(N82,s);diff(N83,s);diff(N84,s);diff(N85,s)'...
                diff(N86,s);diff(N87,s);diff(N88,s)];
            N8t = [diff(N81,t);diff(N82,t);diff(N83,t);diff(N84,t);diff(N85,t)'...
                diff(N86,t);diff(N87,t);diff(N88,t)];
        end
        function B = obj.creatBmatrix(Nix,Niy,N8ix,N8iy,N_Lagrange)
            b1=[Nix(1) 0 0 0 0;0 Niy(1) 0 0 0;Niy(1) Nix(1) 0 0 0;
                0 0 0 0 -Nix(1);0 0 0 Niy(1) 0;0 0 0 Nix(1) -Niy(1);
                0 0 -N8iy(1) N_Lagrange(1) 0;0 0 -N8ix(1) 0 -N_Lagrange(1)];
            b2=[Nix(2) 0 0 0 0;0 Niy(2) 0 0 0;Niy(2) Nix(2) 0 0 0;
                0 0 0 0 -Nix(2);0 0 0 Niy(2) 0;0 0 0 Nix(2) -Niy(2);
                0 0 -N8iy(2) N_Lagrange(2) 0;0 0 -N8ix(2) 0 -N_Lagrange(2)];
            b3=[Nix(3) 0 0 0 0;0 Niy(3) 0 0 0;Niy(3) Nix(3) 0 0 0;
                0 0 0 0 -Nix(3);0 0 0 Niy(3) 0;0 0 0 Nix(3) -Niy(3);
                0 0 -N8iy(3) N_Lagrange(3) 0;0 0 -N8ix(3) 0 -N_Lagrange(3)];
            b4=[Nix(4) 0 0 0 0;0 Niy(4) 0 0 0;Niy(4) Nix(4) 0 0 0;
                0 0 0 0 -Nix(4);0 0 0 Niy(4) 0;0 0 0 Nix(4) -Niy(4);
                0 0 -N8iy(4) N_Lagrange(4) 0;0 0 -N8ix(4) 0 -N_Lagrange(4)];
            b5=[Nix(5) 0 0 0 0;0 Niy(5) 0 0 0;Niy(5) Nix(5) 0 0 0;
                0 0 0 0 -Nix(5);0 0 0 Niy(5) 0;0 0 0 Nix(5) -Niy(5);
                0 0 -N8iy(5) N_Lagrange(5) 0;0 0 -N8ix(5) 0 -N_Lagrange(5)];
            b6=[Nix(6) 0 0 0 0;0 Niy(6) 0 0 0;Niy(6) Nix(6) 0 0 0;
                0 0 0 0 -Nix(6);0 0 0 Niy(6) 0;0 0 0 Nix(6) -Niy(6);
                0 0 -N8iy(6) N_Lagrange(6) 0;0 0 -N8ix(6) 0 -N_Lagrange(6)];
            b7=[Nix(7) 0 0 0 0;0 Niy(7) 0 0 0;Niy(7) Nix(7) 0 0 0;
                0 0 0 0 -Nix(7);0 0 0 Niy(7) 0;0 0 0 Nix(7) -Niy(7);
                0 0 -N8iy(7) N_Lagrange(7) 0;0 0 -N8ix(7) 0 -N_Lagrange(7)];
            b8=[Nix(8) 0 0 0 0;0 Niy(8) 0 0 0;Niy(8) Nix(8) 0 0 0;
                0 0 0 0 -Nix(8);0 0 0 Niy(8) 0;0 0 0 Nix(8) -Niy(8);
                0 0 -N8iy(8) N_Lagrange(8) 0;0 0 -N8ix(8) 0 -N_Lagrange(8)];
            b9=[Nix(9) 0 0 0;0 Niy(9) 0 0;Niy(9) Nix(9) 0 0;0 0 0 -Nix(9);
                0 0 Niy(9) 0;0 0 Nix(9) -Niy(9);0 0 N_Lagrange(9) 0;0 0 0 -N_Lagrange(9)];
            B = [b1 b2 b3 b4 b5 b6 b7 b8 b9];
        end
    end
end
