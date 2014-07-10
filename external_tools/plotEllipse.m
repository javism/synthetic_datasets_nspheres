% Author: Jacob
% Source: http://stackoverflow.com/questions/2153768/draw-ellipse-and-ellipsoid-in-matlab/2155162#2155162
% No licence available

function h=plotEllipse(C,S,color)

    % range to plot over
    %------------------------------------
    N = 100;
    theta = 0:1/N:2*pi+1/N;

    % Parametric equation of the ellipse
    %----------------------------------------
    state(1,:) = S(1)*cos(theta); 
    state(2,:) = S(2)*sin(theta);

    % Coordinate transform (since your ellipse is axis aligned)
    %----------------------------------------
    X = state;
    X(1,:) = X(1,:) + C(1);
    X(2,:) = X(2,:) + C(2);

    % Plot
    %----------------------------------------
    h = plot(X(1,:),X(2,:),'color',color);
    hold on;
    %plot(C(1),C(2),'ro');
    axis equal;
    grid;

end
