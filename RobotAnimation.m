close all
clear
clc
%% Load Data
data = load('Mipsi.1');
% Set up the data points as Matrices of [X Y Z]
pointAO = data(:, 2:4);
pointBA = data(:, 5:7);
pointCB = data(:, 8:10);
pointDC = data(:, 11:13);
pointED = data(:, 14:16);
pointFE = data(:, 17:19);
pointQ = data(:, 20:22);
pointR = data(:, 23:25);

%% Set up 3D graph
figure
plot3([],[],[])


bodyALine = line([pointAO(1,1), pointBA(1,1)], [pointAO(1,2), pointBA(1,2)], [pointAO(1,3), pointBA(1,3)],'Color','r','LineWidth', 3.0);

bodyBLine = line([pointBA(1,1), pointCB(1,1)], [pointBA(1,2), pointCB(1,2)], [pointBA(1,3), pointCB(1,3)],'Color','g','LineWidth', 3.0);

bodyCLine = line([pointCB(1,1), pointDC(1,1)], [pointCB(1,2), pointDC(1,2)], [pointCB(1,3), pointDC(1,3)],'Color','b','LineWidth', 3.0);

bodyDLine = line([pointDC(1,1), pointED(1,1)], [pointDC(1,2), pointED(1,2)], [pointDC(1,3), pointED(1,3)],'Color','b','LineWidth', 3.0);

bodyELine = line([pointED(1,1), pointFE(1,1)], [pointED(1,2), pointFE(1,2)], [pointED(1,3), pointFE(1,3)],'Color','b','LineWidth', 3.0);

bodyFLine = line([pointFE(1,1), pointQ(1,1)] , [pointFE(1,2), pointQ(1,2)] , [pointFE(1,3), pointQ(1,3)] ,'Color','b','LineWidth', 3.0);


pathQ = animatedline;
addpoints(pathQ, pointQ(1,1), pointQ(1,2), pointQ(1,3)); 
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')
axisHandle = gca(); % gca = Get Current Axis
set(axisHandle, 'XLim', [0 3], 'YLim', [0 2], 'ZLim', [-2 2])
set(axisHandle, 'XLimMode', 'manual', 'YLimMode', 'manual', 'XLimMode', 'manual')
axis equal
grid on
view(3)
camup([0 1 0])
for iter=2:length(pointAB)
    % Update the XData, YData and ZData for each line
    set(bodyALine, 'XData', [pointAO(iter,1) , pointBA(iter,1)])
    set(bodyALine, 'YData', [pointAO(iter,2) , pointBA(iter,2)])
    set(bodyALine, 'ZData', [pointAO(iter,3) , pointBA(iter,3)])

    set(bodyBLine, 'XData', [pointBA(iter,1), pointCB(iter,1)])
    set(bodyBLine, 'YData', [pointBA(iter,2), pointCB(iter,2)])
    set(bodyBLine, 'ZData', [pointBA(iter,3), pointCB(iter,3)])

    set(bodyCLine, 'XData', [pointCB(iter,1), pointDC(iter,1)])
    set(bodyCLine, 'YData', [pointCB(iter,2), pointDC(iter,2)])
    set(bodyCLine, 'ZData', [pointCB(iter,3), pointDC(iter,3)])

    set(bodyDLine, 'XData', [pointDC(iter,1), pointED(iter,1)])
    set(bodyDLine, 'YData', [pointDC(iter,2), pointED(iter,2)])
    set(bodyDLine, 'ZData', [pointDC(iter,3), pointED(iter,3)])

    set(bodyELine, 'XData', [pointED(iter,1), pointFE(iter,1)])
    set(bodyELine, 'YData', [pointED(iter,2), pointFE(iter,2)])
    set(bodyELine, 'ZData', [pointED(iter,3), pointFE(iter,3)])

    set(bodyFLine, 'XData', [pointFE(iter,1), pointQ(iter,1)])
    set(bodyFLine, 'YData', [pointFE(iter,2), pointQ(iter,2)])
    set(bodyFLine, 'ZData', [pointFE(iter,3), pointQ(iter,3)])

    
    % Add Points to Animated Line for Q
    addpoints(pathQ, pointQ(iter,1), pointQ(iter,2), pointQ(iter,3)); 
    
    % Tell Matlab to draw
    drawnow
    % Pause for Rate 
    % If there are 800 Frames and we want it to be 8 seconds real time, we
    % draw 100 frames a second.
%     pause(1/100)
    
end%     pause(1/100)
    