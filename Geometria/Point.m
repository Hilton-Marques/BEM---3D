classdef Point < handle
    properties
        coord
        id
        elements = Face.empty;
    end
    methods
        function this = Point(coord,id)
            if nargin == 2
                this.coord = coord;
                this.id = id;
            end
        end
        function out = getSolidAngle(this)
            out = 0.25;
        end
    end
end