backwardRayTracing;

figure(1);
hold on
for i = 1:sensorSizeX
    list = find(output{2,1} == i);
    scatter(output{2,2}(list),output{2,3}(list));
end
hold off