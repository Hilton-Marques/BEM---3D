classdef BEM3D < handle
    properties
        Nodes;
        nNodes;
        nT; %Neuman boundary condition
        nD; %Dirichlet boundary condition
        elements;
        ng = 3;
    end
    methods(Static)
        function [X,Y,Wx,Wy] = triQuad(N)
            v = [0 0; 0 1; 1 0];
            n=1:N;  nnk=2*n+1; A=[1/3 repmat(1,1,N)./(nnk.*(nnk+2))];
            n=2:N; nnk=nnk(n); B1=2/9; nk=n+1; nnk2=nnk.*nnk;
            B=4*(n.*nk).^2./(nnk2.*nnk2-nnk2); ab=[A' [2; B1; B']]; s=sqrt(ab(2:N,2));
            [V,X]=eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
            [X,I]=sort(diag(X)); x=(X+1)/2; wx=ab(1,2)*V(1,I)'.^2/4;
            N=N-1; N1=N+1; N2=N+2;  y=cos((2*(N:-1:0)'+1)*pi/(2*N+2));
            L=zeros(N1,N2);  y0=2;  iter=0;
            while max(abs(y-y0))>eps
                L(:,1)=1;    L(:,2)=y;
                for k=2:N1
                    L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
                end
                Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
                y0=y;    y=y0-L(:,N2)./Lp;  iter=iter+1;
            end
            cd=[ 1, 0, 0; -1, 0, 1; 0, 1,-1]*v;
            t=(1+y)/2;  Wx=abs(det(cd(2:3,:)))*wx;  Wy=1./((1-y.^2).*Lp.^2)*(N2/N1)^2;
            [tt,xx]=meshgrid(t,x); yy=tt.*xx;
            X=cd(1,1)+cd(2,1)*xx+cd(3,1)*yy;    Y=cd(1,2)+cd(2,2)*xx+cd(3,2)*yy;
        end
        function points = relabel(points)
            for i = 1:length(points)
                points(i).idTemp = i;
            end
        end
    end
    methods
        function this = BEM3D(nodes,elements)
            this.Nodes = nodes;
            this.nNodes = length(nodes);
            this.relabel(this.Nodes);
            % count number of neuman boundary condition
            this.nT = 0;
            this.elements = elements;
            for i = 1:this.nNodes
                node = this.Nodes(i);
                this.nT = this.nT + length(node.elements);
            end
            this.nD = this.nNodes;
            this.Calculate();
        end
        function Calculate(this)
            indexU = [];
            indexQ = [];
            bQ = [];
            bU = [];
            HH = this.getH(this.Nodes);
            GG = this.getG(this.Nodes);
            t = zeros(this.nT,1);
            count = 1;
            for i = 1:length(this.Nodes)
                p = this.Nodes(i);
                for j = 1:length(p.elements)
                    el = p.elements(j);
                    if el.q(1) == 2 %Neumann condition
                        t(count) = el.q(2);
                    end
                    count = count + 1;
                end
            end
            d = HH\(GG*t)
            %             bU = zeros(this.nNodes,1);
            %             bQ = zeros(this.nNodes,1);
            %             HHtemp = zeros(this.nNodes,this.nNodes);
            %             GGtemp = zeros(this.nNodes,this.nNodes);
            %             A = zeros(this.nNodes,this.nNodes);
            %             for i = 1:this.nNodes
            %                 element = this.Elements(i);
            %                 if (element.q(1) == 1)
            %                     indexU(end+1) = i;
            %                     bU(i) = element.q(2);
            %                 else
            %                     indexQ(end+1) = i;
            %                     bQ(i) = element.q(2);
            %                 end
            %             end
            %             if(~isempty(indexQ))
            %                 HHtemp(indexQ,indexQ) = HH(indexQ,indexQ);
            %             end
            %             if(~isempty(indexU))
            %                 GGtemp(indexU,indexU) = GG(indexU,indexU);
            %             end
            %             b = HH*bU - GG*bQ;
            %             fileID = fopen('u.txt','w');
            %             fprintf(fileID,'%12.9f\n',b);
            %             fclose(fileID);
            %             %             A = GGtemp - HHtemp ;
            %             %             u = A\b
            %             %             nd = zeros(this.nEl,1);
            %             %             nq = zeros(this.nEl,1);
            %             %             nd(indexU,1) = bU(indexU);
            %             %             nd(indexQ,1) = u(indexQ);
            %             %             nq(indexQ,1) = bQ(indexQ);
            %             %             nq(indexU,1) = u(indexU);
            %             %             HH = this.FormHH(this.intNodes);
            %             %             GG = this.FormGG(this.intNodes);
            %             %             uint = GG*nq - HH*nd
        end
        tic
        function HH = getH(this,sources)
            nSource = length(sources);
            HH = zeros(nSource,this.nD);
            HHH = zeros(nSource,this.nD);
            [X,Y,Wx,Wy] = this.triQuad(this.ng);
            tic
            for i = 1:length(this.elements)
                element = this.elements(i);
                for j = 1:length(sources)
                    source = sources(j);
                    points = element.getPoints();
                    for m = 1:length(points)
                        field = points(m);
                        r = (source.coord - field.coord);
                        d = dot(r,r);
                        if (d > 0)
                            for l = 1:this.ng
                                for k = 1:this.ng
                                    xksi = zeros(3,1);
                                    shapeF = element.getShapeF(X(k,l),Y(k,l));
                                    for p = 1:length(points)
                                        xksi = xksi + shapeF(p)*points(p).coord;
                                    end
                                    r = xksi - source.coord;
                                    jacobian = element.getJacobian(X(k,l),Y(k,l));
                                    n = element.getNormal(X(k,l),Y(k,l));
                                    HH(source.idTemp,field.idTemp) = HH(source.idTemp,field.idTemp) + ...
                                        (jacobian/0.5)*...
                                        (-dot(r,n)/(4*pi*norm(r)^3))*shapeF(m)*...
                                        Wx(k)*Wy(l);
                                end
                            end
                        end
                    end
                end
            end
            toc
            tic
            for i = 1: nSource
                source = sources(i);
                for j = 1:this.nNodes
                    field = this.Nodes(j);
                    r = (source.coord - field.coord);
                    d = dot(r,r);
                    if (d > 0)
                        nEl = length(field.elements);
                        for m = 1:nEl
                            element = field.elements(m);
                            points = element.getPoints();
                            for l = 1:this.ng
                                for k = 1:this.ng
                                    xksi = zeros(3,1);
                                    shapeF = element.getShapeF(X(k,l),Y(k,l));
                                    for p = 1:length(points)
                                        xksi = xksi + shapeF(p)*points(p).coord;
                                    end
                                    r = xksi - source.coord;
                                    jacobian = element.getJacobian(X(k,l),Y(k,l));
                                    n = element.getNormal(X(k,l),Y(k,l));
                                    shapePt = element.getShapePt(field.id,X(k,l),Y(k,l));
                                    HHH(i,j) = HHH(i,j) + (jacobian/0.5)*...
                                        (-dot(r,n)/(4*pi*norm(r)^3))*shapePt*...
                                        Wx(k)*Wy(l);
                                end
                            end
                        end
                    else
                        HHH(i,j) = 0;
                    end
                end
            end
            toc
        end
        function GG = getG(this,sources)
            nSource = length(sources);
            GG = zeros(nSource,this.nT);
            GGG = zeros(nSource,this.nT);
            [X,Y,Wx,Wy] = this.triQuad(this.ng);
            count = 1;
            tic
            for i = 1:length(this.elements)
                element = this.elements(i);
                points = element.getPoints();
                for j = 1:length(points)
                    for m = 1:length(sources)
                        source = sources(m);
                        for l = 1:this.ng
                            for k = 1:this.ng
                                xksi = zeros(3,1);
                                shapeF = element.getShapeF(X(k,l),Y(k,l));
                                for p = 1:length(points)
                                    xksi = xksi + shapeF(p)*points(p).coord;
                                end
                                r = xksi - source.coord;
                                jacobian = element.getJacobian(X(k,l),Y(k,l));
                                n = element.getNormal(X(k,l),Y(k,l));
                                GG(source.idTemp,count) = GG(source.idTemp,count) + ...
                                    (jacobian/0.5)*...
                                    (1/(4*pi*norm(r)))*shapeF(j) *...
                                    Wx(k)*Wy(l);
                            end
                        end
                    end
                    count = count + 1;
                end
            end
            toc
            tic
            for i = 1: nSource
                source = sources(i);
                count = 1; % counter to deal with many normals at a node
                for j = 1:this.nNodes
                    field = this.Nodes(j);
                    nEl = length(field.elements);
                    for m = 1:nEl
                        element = field.elements(m);
                        points = element.getPoints();
                        for l = 1:this.ng
                            for k = 1:this.ng
                                xksi = zeros(3,1);
                                shapeF = element.getShapeF(X(k,l),Y(k,l));
                                for p = 1:length(points)
                                    xksi = xksi + shapeF(p)*points(p).coord;
                                end
                                r = xksi - source.coord;
                                jacobian = element.getJacobian(X(k,l),Y(k,l));
                                shapePt = element.getShapePt(field.id,X(k,l),Y(k,l));
                                GGG(i,count) = GGG(i,count) + (jacobian/0.5)*...
                                    (1/(4*pi*norm(r)))*shapePt *...
                                    Wx(k)*Wy(l);
                            end
                        end
                    end
                    count = count + 1;
                end
            end
            toc
        end
    end
    
end
