classdef Solid < handle
    properties
        points;
        edges;
        heds;
        elementsMother;
        elements = Face.empty;
        nL;
        nHeds;
        nPts;
        nEdges;
        nEl;
        elementsLevel;
        typeEl;
    end
    methods
        function this = Solid(points,edges,heds,elements,nL,typeEl)
            this.points = points;
            this.edges = edges;
            this.heds = heds;
            this.elementsMother = elements;
            this.nL = nL;
            this.nHeds = length(this.heds);
            this.nPts = length(this.points);
            this.nEdges = length(this.edges);
            this.nEl = length(elements);
            this.elementsLevel = cell(nL,1);
            this.typeEl = typeEl;
            this.start();
            this.initializeEl();
        end
        function start(this)
            this.elementsLevel{1} = this.elementsMother;
            for i = 2:this.nL
                Elements_m = this.elementsLevel{i-1};
                nElNv = length(Elements_m);
                newElements(4*nElNv) = Face;
                count = 1;
                for j = 1:nElNv
                    el = Elements_m(j);
                    localHeds = Hed.empty;
                    localPoints = Point.empty;
                    for k = 1:length(el.heds)
                        hed_m = el.heds(k);
                        edge_m = this.edges(hed_m.edgeId);
                        if ~edge_m.isSplited
                            coord = (this.points(hed_m.inc(1)).coord + this.points(hed_m.inc(2)).coord )/2;
                            this.nPts = this.nPts + 1;
                            pt = Point(coord,this.nPts);
                            this.points(end+1) = pt;
                            localPoints(end+1) = pt;
                            for l = 1:2 % number of subdivisions
                                this.nEdges = this.nEdges + 1;
                                this.nHeds = this.nHeds + 1;
                                inc1 = [hed_m.inc(l),this.nPts];
                                inc2 = [this.nPts,hed_m.inc(l)];
                                incs = [inc1;inc2];
                                newHed = Hed(incs(l,:),this.nHeds,this.nEdges);
                                this.nHeds = this.nHeds + 1;
                                newHedTwin = Hed(incs(mod(l,2)+1,:),this.nHeds,this.nEdges);
                                newEdge = Edge(newHed,newHedTwin,this.nEdges);
                                localHeds(end+1) = newHed;
                                this.heds(end+1) = newHed;
                                this.heds(end+1) = newHedTwin;
                                this.edges(end+1) = newEdge;
                                edge_m.childrenId(l) = this.nEdges;
                            end
                            edge_m.isSplited = true;
                            edge_m.mid = pt;
                        else
                            localPoints(end+1) = edge_m.mid;
                            for l = 2:-1:1 % the this.heds reverse the order
                                edge_new = this.edges(edge_m.childrenId(l));
                                localHeds(end+1) = edge_new.hed2;
                            end
                        end
                    end
                    internalHeds(3) = Hed;
                    for k = 1:3
                        this.nEdges = this.nEdges + 1;
                        this.nHeds = this.nHeds + 1;
                        p2 = localPoints(mod(k,3)+1).id;
                        p1 = localPoints(k).id;
                        newHed = Hed([p1,p2],this.nHeds,this.nEdges);
                        this.nHeds = this.nHeds + 1;
                        newHedTwin = Hed([p2,p1],this.nHeds,this.nEdges);
                        internalHeds(k) = newHed;
                        this.heds(end+1) = newHed;
                        this.heds(end+1) = newHedTwin;
                        newEdge = Edge(newHed,newHedTwin,this.nEdges);
                        this.edges(end+1) = newEdge;
                    end
                    %reorder
                    internalHeds = internalHeds([3,1,2]) ;
                    % Update hNext for border triangles
                    c1 = 1;
                    c2 = 5;
                    for l = 1:2:5
                        this.nEl = this.nEl + 1;
                        firstId = localHeds(l).id;
                        idMidle = this.edges(internalHeds(c1).edgeId).hed2.id;
                        localHeds(l).heNext = idMidle;
                        localHeds(l).elId = this.nEl;
                        idBefore = mod(6 + l - 2,6)+1;
                        this.heds(idMidle).heNext = localHeds(idBefore).id;
                        this.heds(idMidle).elId = this.nEl;
                        localHeds(mod(c2,6)+1).heNext = firstId;
                        localHeds(mod(c2,6)+1).elId = this.nEl;
                        c1 = c1 + 1;
                        c2 = c2 + 2;
                    end
                    % Update hNext for internal triangle
                    this.nEl = this.nEl + 1;
                    for l = 1:3
                        internalHeds(l).heNext = internalHeds(mod(l,3)+1).id;
                        internalHeds(l).elId = this.nEl;
                    end
                    %Create new Elements
                    elementsChildren(4) = Face;
                    elementsChildren(1) = Face(localHeds(1).id,this.heds,el.q);
                    elementsChildren(2) = Face(localHeds(3).id,this.heds,el.q);
                    elementsChildren(3) = Face(localHeds(5).id,this.heds,el.q);
                    elementsChildren(4) = Face(internalHeds(1).id,this.heds,el.q);
                    newElements(count:count+3) = elementsChildren;
                    el.children = elementsChildren;
                    count = count + 4;
                end
                % Create interpolation functions
                for k = 1:length(newElements)
                    el = newElements(k);
                    newPoints = el.refine(this.edges,this.points,this.typeEl);
                    % As new this.points was created there is a need to updtate this.nPts
                    this.points(end+1:end+length(newPoints)) = newPoints;
                    this.nPts = length(this.points);
                end
                this.elementsLevel{i} = newElements;
            end
            % vectorize elements
            this.elements(this.nEl) = Face;
            count = 1;
            for i = 1:length(this.elementsLevel)
                els = this.elementsLevel{i};
                this.elements(count:count+length(els)-1) = els;
                count = count+length(els);
            end
        end
        function initializeEl(this)
            n = length(this.elements);
            nMother = length(this.elementsMother);
            %The mother elements doesnt need to initialize
            for i = nMother+1:n
                el = this.elements(i);
                el.initialize(this.points,this.heds,this.edges,this.elements);
            end
        end
        function points = getPointsInLevel(this,nv)
            elements = this.elementsLevel{nv};
            ids = [];
            points = Point.empty;
            for i = 1:length(elements)
                el = elements(i);
                for j = 1:length(el.refinedPoints)
                    pi = el.refinedPoints(j);
                    if ~ismember(pi.id,ids)
                        pi.elements(end+1) = el;
                        points(end+1) = pi;
                        ids(end+1) = pi.id;
                    else
                        pi.elements(end+1) = el;
                    end
                end
                if strcmp(this.typeEl,'Const')
                    continue;
                end
                for j = 1:length(el.points)
                    pi = el.points(j);
                    if ~ismember(pi.id,ids)
                        pi.elements(end+1) = el;
                        points(end+1) = pi;
                        ids(end+1) = pi.id;
                    else
                        pi.elements(end+1) = el;
                    end
                end
            end
        end
        function plot(this,nv)
            if nargin == 1
                flagRefined = false;
            end
            axis equal;
            ax = gca;
            set(gca,'XColor', 'none','YColor','none','ZColor','none')
            view(30,70);
            elementsNv  = this.elementsLevel{nv};
            for i = 1:length(elementsNv)
                el = elementsNv(i);
                el.plot;
            end
        end
    end
end

