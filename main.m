% Clear workspace
clear
close(findall(0,'Type','figure'));
clc;
addpath('Geometria','Solver');
%% Domain
% Piramid
nHeds = 12;
nPts = 4;
nEdges = 6;
nEl = 4;
%Boundary condition per face
q1 = [2;1];
q2 = [2;2];
q3 = [2;3];
q4 = [2;4];
% Type of element
typeEl = 'T3'; % Const, T3, T6, T9
points(4) = Point;
points(1) = Point([-3;-1.5;0],1);
points(2) = Point([3;-1.5;0],2);
points(3) = Point([0;1.5;0],3);
points(4) = Point([0;0;4],4);
heds(6) = Hed;   % number of half edges
heds(1) = Hed([1,2],1,1,1,2);
heds(2) = Hed([2,3],2,2,1,3);
heds(3) = Hed([3,1],3,3,1,1);
heds(4) = Hed([4,1],4,4,2,5);
heds(5) = Hed([1,3],5,3,2,6);
heds(6) = Hed([3,4],6,6,2,4);
heds(7) = Hed([3,2],7,2,3,8);
heds(8) = Hed([2,4],8,5,3,9);
heds(9) = Hed([4,3],9,6,3,7);
heds(10) = Hed([2,1],10,1,4,11);
heds(11) = Hed([1,4],11,4,4,12);
heds(12) = Hed([4,2],12,5,4,10);
elements(4) = Face;
elements(1) = Face(1,heds,q1);
elements(2) = Face(4,heds,q2);
elements(3) = Face(7,heds,q3);
elements(4) = Face(10,heds,q4);
edges(6) = Edge;
edges(1) = Edge(heds(1),heds(10),1);
edges(2) = Edge(heds(2),heds(9),2);
edges(3) = Edge(heds(3),heds(5),3);
edges(4) = Edge(heds(4),heds(11),4);
edges(5) = Edge(heds(8),heds(12),5);
edges(6) = Edge(heds(6),heds(9),6);
nL = 4; %4 levels
solid = Solid(points,edges,heds,elements,nL,typeEl);
nv = 2;
hold on;
solid.plot(nv,true); %by level
nodes = solid.getPointsInLevel(nv);
solver = Solver(nodes);

keyboard;
