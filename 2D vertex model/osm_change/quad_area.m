function A = quad_area(P1, P2, P3, P4)
    % P1, P2, P3, P4 The vertices of the quadrilateral are arranged in order
    x = [P1(1); P2(1); P3(1); P4(1); P1(1)];
    y = [P1(2); P2(2); P3(2); P4(2); P1(2)];
    A = 0.5 * abs(sum(x(1:end-1).*y(2:end) - y(1:end-1).*x(2:end)));
end