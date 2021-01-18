classdef FMM < handle
    properties
        nLMax
        elementsLevel
        nTotalPoints
        n % truncation terms
        solver
    end
    methods
        function this = FMM(solid,n,ng)
            this.nLMax = solid.nL;
            this.elementsLevel = solid.elementsLevel;
            this.nTotalPoints = solid.nPts;
            this.n = n;
            this.solver = Solver(ng,n,this.nTotalPoints);
            
            %% Upward Pass
            for i = this.nLMax:-1:1
                elements = this.elementsLevel{i};
                if i == this.nLMax
                    for j = 1:length(elements)
                        element = elements(j);
                        % Calculate multipole expansion Eq. 3.53 e 3.54
                        this.solver.calculateME(element)
                    end
                else
                    for j = 1:length(elements)
                        element = elements(j);
                        % Eq. 3.55
                        this.solver.calculateMMT(element)
                    end
                end
            end
            %% Downward Pass
            keyboard;
        end
    end
end