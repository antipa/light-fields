backwardRayTracing;

gridX = 4;
gridT = 4;
xRange = [0 100];
tRange = [-4 4];
stepX = (max(xRange) - min(xRange)) ./ gridX;
stepT = (max(tRange) - min(tRange)) ./ gridT;
xValues = [];
for l = xRange(1):stepX:xRange(2)
    for q = 1:gridX
        xValues = [xValues l];
    end
end
yValues = [];
for l = tRange(2):-stepT:tRange(1)
    for q = 1:gridT
        yValues = [yValues l];
    end
end
    figure(1);
    hold on
    grid on
    aMatrix = zeros(sensorSizeX, gridX .* gridT);
    for i = 1:sensorSizeX
        for j = 1: gridX * gridT
            xmin = xValues(j);
            xmax = xValues(j+gridX);
            ymax = yValues(j);
            ymin = yValues(j+gridT);
            list = find(output{2,1} == i);
            scatter(output{2,2}(list),output{2,3}(list));
            count = 0;
            for k = 1:length(list)
                if output{2,2}(list(k)) > xmin && output{2,2}(list(k)) < xmax && output{2,3}(list(k)) > ymin && output{2,3}(list(k)) < ymax
                    count = count + 1;
                end
            end
            aMatrix(i,j) = count ./ raysPerPixel;
        end
    end
    hold off
