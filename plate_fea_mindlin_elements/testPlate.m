clc; clear all; close all; format compact; format long
problem = 1;
addpath ('.\breps');
elasticityClass = @plateFEMindlin; % QuadElasticity or TriElasticity
if (problem == 1) % 
    fem = elasticityClass('Rectangle.brep',500,'Quadratic');
    fem = fem.fixEdge(6);fem = fem.fixEdge(5);fem = fem.fixEdge(4);
    fem = fem.fixEdge(3);fem = fem.fixEdge(2);fem = fem.fixEdge(1);
    
    %fem = fem.applyYForceOnEdge(3,-1); %no surface force 
elseif  (problem == 2) % 
    fem = elasticityClass('SquareHole.brep',2000,'Linear','PlaneStrain');
    fem = fem.fixEdge(9);
    fem = fem.applyXForceOnEdge(7,1);
elseif  (problem == 3) %     
    fem = elasticityClass('LBracket.brep',2000,'Quadratic','PlaneStress');
    fem = fem.fixEdge(7);
    fem = fem.applyYForceOnEdge(3,-1); 
elseif  (problem == 4) % 
    fem = elasticityClass('tube.brep',500,'Quadratic','AxiSymmetric');
    fem = fem.fixEdge(1);
    fem = fem.applyXForceOnEdge(4,1); 
elseif  (problem == 5) % 
    fem = elasticityClass('Square.brep',500,'Quadratic','PlaneStress');
    % The number '1' of fixed edge is defined in *.brep file by starting and ending points
    fem = fem.fixEdge(1); 
    fem = fem.applyYForceOnEdge(3,1); 
    % applyYForceOnEdge=Y-direction force; applyXForceOnEdge= X-direction
    % force; '3'=# of edge; '1'=value of force
elseif  (problem == 6) % 
    fem = elasticityClass('CantileverBeam.brep',1500,'Quadratic','PlaneStress');
    fem = fem.fixEdge(6);
    fem = fem.applyYForceOnEdge(3,-1); 
elseif  (problem == 7) % 
    fem = elasticityClass('MichelleBeam.brep',500,'Quadratic','PlaneStress');
    fem = fem.fixEdge(7);
    fem = fem.applyYForceOnEdge(3,-1); 
elseif  (problem == 8) % 
    fem = elasticityClass('Mast.brep',500,'Quadratic','PlaneStress');
    fem = fem.fixEdge(1);
    fem = fem.applyYForceOnEdge([3 9],-1);     
elseif  (problem == 9) % 
    fem = elasticityClass('RectanglePlateThreeHoles.brep',2000,'Linear','PlaneStrain');
    fem = fem.fixEdge(4);
    fem = fem.applyXForceOnEdge(2,1);
end

%fem.plotMesh();
fem = fem.runFE(); 

return;
figure; fem.plotDeformation();
figure; fem.plotStress();
