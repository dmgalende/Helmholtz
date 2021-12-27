function [s1,s2,s3] = interpolationSplines(tauL,tauR,x,xL,xR)

a = -1;
b = 1;
tx = tmap(x,a,b);
s1 = spline1(tx);
s2 = spline2(tx);
s3 = spline3(tx);

    function f = beta(u)
        f = 0.25*u.^3 - 0.75*u + 0.5;
    end

    function f = spline1(u)
        f = zeros(size(u));
        pos1 = u <= -1;
        pos2 = u >= -tauL;
        pos3 = ~pos1 & ~pos2;
        f(pos3) = beta((1-tauL)^-1 * (2*u(pos3)+tauL+1));
        f(pos1) = 1;
        f(pos2) = 0;
    end

    function f = spline3(u)
        f = zeros(size(u));
        pos1 = u <= tauR;
        pos2 = u >= 1;
        pos3 = ~pos1 & ~pos2;
        f(pos3) = beta((1-tauR)^-1 * (-2*u(pos3)+tauR+1));
        f(pos1) = 0;
        f(pos2) = 1;
    end

    function f = spline2(u)
        f = zeros(size(u));
        pos1 = u <= -1;
        pos2 = ~pos1 & (u <= -tauL);
        pos3_ = pos1 | pos2;
        pos3 = ~pos3_ & (u <= tauR);
        pos4_ = pos3_ | pos3;
        pos4 = ~pos4_ & (u <= 1);
        pos5 = ~(pos4_ | pos4);
        f(pos2) = beta(-(1-tauL)^-1 * (2*u(pos2)+tauL+1));
        f(pos4) = beta((1-tauR)^-1 * (2*u(pos4)-tauR-1));
        f(pos1) = 0;
        f(pos3) = 1;
        f(pos5) = 0;
    end

    function f = tmap(u,a,b)
        f = (xL-xR)^-1 * ((a-b)*u + b*xL - a*xR);
    end

end

