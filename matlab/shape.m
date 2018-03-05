function c = shape(sp, tp, thickness, mid, plt)
% Generate shape corners given starting point, target point, thickness and
% midpoint(0 if quadrilateral)
if ~mid
    coeff = polyfit([sp(1),tp(1)],[sp(2),tp(2)],1);
    upper = coeff + [0, 0.5*thickness];
    lower = coeff - [0, 0.5*thickness];
    
    b = sp(2) + sp(1)/upper(1);
    perp = [-1/upper(1), b];
    b2 = 0.5*thickness/sin(atan(upper(1)));
    l2 = perp + [0, -b2];
    
    b = tp(2) + tp(1)/upper(1);
    perp = [-1/upper(1), b];
    u2 = perp + [0, b2];
    
    c = zeros(4,2);
    x = (upper(2)-l2(2))/(l2(1)-upper(1));
    y = upper(1)*x + upper(2);
    c(1,:) = [x,y];
    
    x = (upper(2)-u2(2))/(u2(1)-upper(1));
    y = upper(1)*x + upper(2);
    c(2,:) = [x,y];
    
    x = (lower(2)-u2(2))/(u2(1)-lower(1));
    y = lower(1)*x + lower(2);
    c(3,:) = [x,y];
    
    x = (lower(2)-l2(2))/(l2(1)-lower(1));
    y = lower(1)*x + lower(2);
    c(4,:) = [x,y];
    
    if plt
        figure;
        plot(c(:,1),c(:,2),'b');
        hold on;
        plot(sp(1),sp(2),'rx');
        hold on;
        plot(tp(1),tp(2),'rx');
    end
else
    c = 0;
end

end

