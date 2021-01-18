classdef Face < handle
    properties
        hedInc;
        pointsInc;
        points;
        heds = Hed.empty;
        adjacentFaces;
        id;
        typeEl;
        refinedPointsInc;
        refinedPoints;
        children;
        q; %boundary condition
        MEH;
        MEG;
    end
    methods
        function this = Face(hedInc,heds,q)
            if nargin > 0
                this.hedInc = hedInc;
                this.getIndcPoints(heds);
                this.id = heds(hedInc).elId;
                this.q = q;
            end
        end
        function getIndcPoints(this,heds)
            hedsEl(3) = Hed;
            pointsInc = zeros(1,3);
            hed = heds(this.hedInc);
            p1 = hed.inc(1);
            p2 = hed.inc(2);
            hedsEl(1) = hed;
            pointsInc(1) = p1;
            c = 1;
            while p2 ~= p1
                hed = heds(hed.heNext);
                c = c + 1;
                hedsEl(c) = hed;
                p2 = hed.inc(2);
                pointsInc(c) = hed.inc(1);
            end
            this.heds = hedsEl;
            this.pointsInc = pointsInc;
        end
        function findAdjacents(this,heds,edges,elements)
            elIds = [];
            id = this.id;
            for i = 1:length(this.heds)
                hed = this.heds(i);
                p2Target = hed.inc(2);
                hed = heds(edges(hed.edgeId).getTwin(hed.id).heNext);
                while true
                    p2 = hed.inc(2);
                    if p2 == p2Target
                        break;
                    end
                    newId = hed.elId;
                    if ~ismember(hed.elId,elIds)
                        elIds(end+1) = newId;
                    end
                    hed = heds(edges(hed.edgeId).getTwin(hed.id).heNext);
                end
            end
            this.adjacentFaces = elements(elIds);
        end
        function newPoints = refine(this,edges,points,typeEl)
            n = length(points);
            this.typeEl = typeEl;
            newPoints = Point.empty;
            switch this.typeEl
                case 'Const'
                    n = n + 1;
                    coord = (points(this.pointsInc(1)).coord + ...
                        points(this.pointsInc(2)).coord + ...
                        points(this.pointsInc(3)).coord) /3;
                    pMid = Point(coord,n);
                    newPoints(end+1) = pMid;
                    this.refinedPointsInc(end+1) = n;
                case 'T3'
                    %do nothing
                case 'T6'
                    for i = 1:length(this.heds)
                        hed = this.heds(i);
                        edge = edges(hed.edgeId);
                        if ~edge.isSplitedToRefine
                            p1 = points(hed.inc(1));
                            p2 = points(hed.inc(2));
                            coord = (p1.coord + p2.coord)/2;
                            n = n+1;
                            pMid = Point(coord,n);
                            newPoints(end+1) = pMid;
                            this.refinedPointsInc(end+1) = n;
                            edge.pRefined = n;
                            edge.isSplitedToRefine = true;
                        else
                            this.refinedPointsInc(end+1) = edge.pRefined;
                        end
                    end
            end
        end
        function initialize(this,points,heds,edges,elements)
            this.points = points(this.pointsInc);
            this.refinedPoints = points(this.refinedPointsInc);
            this.passBoundaryCondition();
            this.findAdjacents(heds,edges,elements);
        end
        function passBoundaryCondition(this)
            for i = 1:length(this.pointsInc)
                pt = this.points(i);
                if isempty(pt.q)
                    pt.q = this.q;
                end
            end
            for i = 1:length(this.refinedPointsInc)
                pt = this.refinedPoints(i);
                if isempty(pt.q)
                    pt.q = this.q;
                end
            end
        end
        function points = getPoints(this)
            points = Point.empty;
            points = this.points;
            if strcmp(this.typeEl,'Const') || strcmp(this.typeEl,'T3')
                return
            else
                points(end+1:end  + length(this.refinedPoints)) = ...
                    this.refinedPoints;
                return
            end
        end
        %% Methods used in BEM
        function shapeF = getShapeF(this,xi,eta)
            switch this.typeEl
                case 'Const'
                    shapeF = zeros(3,1);
                    shapeF(1) = xi;
                    shapeF(2) = eta;
                    shapeF(3) = 1 - xi - eta;
                case 'T3'
                    shapeF = zeros(3,1);
                    shapeF(1) = xi;
                    shapeF(2) = eta;
                    shapeF(3) = 1 - xi - eta;
                case 'T6'
                    shapeF = zeros(6,1);
                    zeta = 1 - xi - eta;
                    shapeF(1) = (2*xi-1)*xi;
                    shapeF(2) = (2*eta - 1)*eta;
                    shapeF(3) = (2*zeta - 1)*zeta;
                    shapeF(4) = 4*xi*eta;
                    shapeF(5) = 4*eta*zeta;
                    shapeF(6) = 4*zeta*xi;
            end
        end
        function dShapeF = getDShapeF(this,xi,eta)
            dShapeF = zeros(3,2);
            switch this.typeEl
                case 'Const'
                    dShapeF(:,1) = this.points(1).coord - this.points(3).coord;
                    dShapeF(:,2) = this.points(2).coord - this.points(3).coord;
                case 'T3'
                    dShapeF(:,1) = this.points(1).coord - this.points(3).coord;
                    dShapeF(:,2) = this.points(2).coord - this.points(3).coord;
                case 'T6'
                    zeta = 1 - xi - eta;
                    dShapeF(:,1) = 4*eta*this.refinedPoints(1).coord + ...
                        2*this.points(1).coord*xi + ...
                        this.points(1).coord*((-1)+2*xi)+...
                        +4*this.refinedPoints(3).coord*zeta;
                    dShapeF(:,2) = 2*eta*this.points(2).coord + ...
                        ((-1)+2*eta)*this.points(2).coord + ...
                        4*this.refinedPoints(1).coord*xi + ...
                        4*this.refinedPoints(2).coord*zeta;
            end
        end
        function jacobian = getJacobian(this,xi,eta)
            der = this.getDShapeF(xi,eta);
            jacobian = norm(cross(der(:,1),der(:,2)));
        end
        function normal = getNormal(this,xi,eta)
            der = this.getDShapeF(xi,eta);
            normal = -(1/(norm(der(:,1))*norm(der(:,1))))* ...
                cross(der(:,1),der(:,2));
        end
        function shapePt = getShapePt(this,id,xi,eta)
            switch this.typeEl
                case 'Const'
                    shapePt = 1;
                case 'T3'
                    [~,loc] = ismember(id,this.pointsInc);
                    shapeF = this.getShapeF(xi,eta);
                    shapePt = shapeF(loc);
                case 'T6'
                    ids = [this.pointsInc,this.refinedPointsInc];
                    [~,loc] = ismember(id,ids);
                    shapeF = this.getShapeF(xi,eta);
                    shapePt = shapeF(loc);
            end
        end
        %% Methods used in FMM
        % Get expasion coordinate point
        function yc = getYc(this)
            yc = zeros(3,1);
            n  = length(this.pointsInc);
            for i = 1:n
                yc = yc + this.points(i).coord;
            end
            yc = yc/n;
        end
        function points = getFieldPoints(this)
            switch this.typeEl
                case 'Const'
                    points = this.refinedPoints;
                case 'T3'
                    points = [this.points,this.refinedPoints];
                case 'T6'
                    points = [this.points,this.refinedPoints];
            end
        end
        function points = geElementPoints(this)
            switch this.typeEl
                case 'Const'
                    points = this.points;
                case 'T3'
                    points = [this.points,this.refinedPoints];
                case 'T6'
                    points = [this.points,this.refinedPoints];
            end
        end
        function plot(this,color)
            if nargin == 1
                color = [0,0,1];
            end
            p = zeros(length(this.pointsInc),3);
            for i = 1:length(this.pointsInc)
                p(i,:) = this.points(i).coord;
            end
            line([p(1,1),p(2,1)],[p(1,2),p(2,2)],...
                [p(1,3),p(2,3)],'Color', color);
            line([p(2,1),p(3,1)],[p(2,2),p(3,2)],...
                [p(2,3),p(3,3)],'Color',color);
            line([p(3,1),p(1,1)],[p(3,2),p(1,2)],...
                [p(3,3),p(1,3)],'Color',color);
            for i = 1:length(this.pointsInc)
                pi = p(i,:);
                plot3(pi(1),pi(2),pi(3),'o','Color','red');
            end
            for i = 1:length(this.refinedPointsInc)
                pi = this.refinedPoints(i).coord;
                plot3(pi(1),pi(2),pi(3),'o','Color','cyan');
            end
            
        end
        function plotAdjacent(this,points,elements)
            n = length(this.adjacentFaces);
            for i = 1:n
                el = elements(this.adjacentFaces(i));
                el.plot(points,[1,0,0]);
            end
        end
    end
end