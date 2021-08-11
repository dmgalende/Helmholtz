function [u1,u2]= nefemInput2DCorrectPeriodicParameters(u1,u2,aNurbs)

if aNurbs.periodic
    TOL = 1e-5;
    du1Ini = abs(u1-aNurbs.iniParam);
    du2Ini = abs(u2-aNurbs.iniParam);
    du1End = abs(u1-aNurbs.endParam);
    du2End = abs(u2-aNurbs.endParam);

    if du1Ini<TOL
        if du2End<du2Ini
            u1 = aNurbs.endParam;
        end
    end

    if du2Ini<TOL
        if du1End<du1Ini
            u2 = aNurbs.endParam;
        end
    end
end