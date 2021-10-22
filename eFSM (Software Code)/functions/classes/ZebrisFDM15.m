classdef ZebrisFDM15
    properties
        FrameRowSize = 176
        FrameColSize = 64
        Frames
        Frequency
        Data
        Plate
    end
    methods
        function Pressure = getPressure(obj)
            Pressure = zeros(size([obj.Data],1)/[obj.FrameRowSize], 1);
            i = 1; j = 1;
            while i < size([obj.Data], 1)
                Pressure(j, 1) = sum(sum(obj.Data(i:i + obj.FrameRowSize - 1, :), 1), 2);
                j = j + 1;
                i = i + 176;
            end
        end
        function F = getForce(obj)
            F = obj.getPressure.*(0.8469^2);
        end
        function M = getMass(obj)
            M = obj.getForce./9.81;
        end
        function A = currentMatrix(obj, i)
            A = obj.Data((i-1)*obj.FrameRowSize + 1:(i)*obj.FrameRowSize,1:obj.FrameColSize);
        end
        function [] = plotcurrentMatrix(obj, i)
            surf(obj.currentMatrix(i)); view(0,90); shading interp;
            title('Zebris FDM 1.5 (Pressure Plate)');
            axis equal
        end
        function [] = showWholeAnimation(obj, range)
            if nargin == 1
                range = 1:obj.Frames;
            end
            for i = range%1:obj.Frames
                obj.plotcurrentMatrix(i);
                axis equal
                drawnow;
            end
        end
    end
end