backwardRayTracing;

gridX = 5;
gridT = 5;
xRange = [0 1000];
tRange = [4 -4];
stepX = (max(xRange) - min(xRange)) ./ gridX;
stepT = (tRange(2) - tRange(1)) ./ gridT;
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
%     scatter(output{2,2}(list),output{2,3}(list));
    for j = 1: gridX * gridT
        xmin = xValues(j);
        xmax = xValues(j+gridX);
        ymin = yValues(mod(j-1,length(yValues)-1) + 1);
        ymax = yValues(mod(j-1,length(yValues)-1) + 2);
        count = 0;
        a = output{2,2}(list) > xmin & output{2,2}(list) <= xmax & ...
            output{2,3}(list) > ymin & output{2,3}(list) <= ymax;
        b = output{2,2}(list) > xRange(1) & output{2,2}(list) <= xRange(2) & ...
            output{2,3}(list) > tRange(2) & output{2,3}(list) <= tRange(1);
        if sum(b)
        aMatrix(i,j) = sum(a) ./ sum(b);
        end
    end
end
hold off
