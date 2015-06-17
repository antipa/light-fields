backwardRayTracing;

figure(1);
hold on
grid on
aMatrix = zeros(sensorSizeX,9);
for i = 1:sensorSizeX
    for j = 1:9
        list = find(output{2,1} == i);
        scatter(output{2,2}(list),output{2,3}(list));
        count = 0;
        for k = 1:length(list)
            switch j
                case 1
                    if output{2,2}(list(k)) > 0 && output{2,2}(list(k)) < 33.33 && output{2,3}(list(k)) > 1.333 && output{2,3}(list(k)) < 4
                        count = count + 1;
                    end
                case 2
                    if output{2,2}(list(k)) > 0 && output{2,2}(list(k)) < 33.33 && output{2,3}(list(k)) > -1.333 && output{2,3}(list(k)) < 1.33
                        count = count + 1;
                    end
                case 3
                    if output{2,2}(list(k)) > 0 && output{2,2}(list(k)) < 33.33 && output{2,3}(list(k)) > -4 && output{2,3}(list(k)) < -1.33
                        count = count + 1;
                    end
                case 4
                    if output{2,2}(list(k)) > 33.33 && output{2,2}(list(k)) < 66.66 && output{2,3}(list(k)) > 1.333 && output{2,3}(list(k)) < 4
                        count = count + 1;
                    end
                case 5
                    if output{2,2}(list(k)) > 33.33 && output{2,2}(list(k)) < 66.66 && output{2,3}(list(k)) > -1.333 && output{2,3}(list(k)) < 1.33
                        count = count + 1;
                    end
                case 6
                    if output{2,2}(list(k)) > 33.33 && output{2,2}(list(k)) < 66.66 && output{2,3}(list(k)) > -4 && output{2,3}(list(k)) < -1.33
                        count = count + 1;
                    end
                case 7
                    if output{2,2}(list(k)) > 66.66 && output{2,2}(list(k)) < 99.99 && output{2,3}(list(k)) > 1.333 && output{2,3}(list(k)) < 4
                        count = count + 1;
                    end
                case 8
                    if output{2,2}(list(k)) > 66.66 && output{2,2}(list(k)) < 99.99 && output{2,3}(list(k)) > -1.333 && output{2,3}(list(k)) < 1.33
                        count = count + 1;
                    end
                case 9
                    if output{2,2}(list(k)) > 66.66 && output{2,2}(list(k)) < 99.99 && output{2,3}(list(k)) > -4 && output{2,3}(list(k)) < -1.33
                        count = count + 1;
                    end
            end
        end
        aMatrix(i,j) = count ./ raysPerPixel;
    end
end
hold off
