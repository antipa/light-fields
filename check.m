backwardRayTracing;

gridX = 3;
gridT = 3;
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
yValues = [tRange(2):-stepT:tRange(1)];

figure(1);
hold on
grid on
aMatrix = zeros(sensorSizeX, gridX .* gridT);
for i = 1:sensorSizeX
    list = find(output{2,1} == i);
    scatter(output{2,2}(list),output{2,3}(list));
    for j = 1: gridX * gridT
        xmin = xValues(j);
        xmax = xValues(j+gridX);
        ymax = yValues(mod(j-1,length(yValues)-1) + 1);
        ymin = yValues(mod(j-1,length(yValues)-1) + 2);
        count = 0;
        a = output{2,2}(list) > xmin & output{2,2}(list) <= xmax & ...
            output{2,3}(list) > ymin & output{2,3}(list) <= ymax;
        aMatrix(i,j) = sum(a) ./ raysPerPixel;
    end
end
hold off
