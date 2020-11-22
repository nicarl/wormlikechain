%% Worm like chain
%% Description
% This function calculates the amplitude form factor  of worm like chains according to Pedersen & Schurtenberger with the corrections according to Rey et al (Chen, W. R., Butler, P. D., & Magid, L. J. (2006). Incorporating intermicellar interactions in the fitting of SANS data from cationic wormlike micelles. Langmuir, 22(15), 6539–6548. https://doi.org/10.1021/la0530440).
%
%
%% CODE
function [val] = kWormlikeChain_sasview(q,a)
%% Parameters
L = a(1);              % contour length
b = a(2);              % Kuhn length
scal = a(3);           % Scaling
%% Calculation
val = Sk_WR(q,L,b)*scal;

end
%% Support functions
function val = AlphaSquare(x)
    val = power(1.0+(x/3.12).^2+(x/8.67).^3, 0.176/3.0);
end

function val = Rgsquarezero(L,b)
    r = b/L;
    val = (L*b/6.0) * (1.0 + r*(-1.5 + r*(1.5 + r*0.75*expm1(-2.0/r))));
end

function val = Rgsquareshort(L,b)
    val = AlphaSquare(L/b) * Rgsquarezero(L,b);
end

function val = Rgsquare(L,b)
    val = AlphaSquare(L/b)*L*b/6.0;
end

function val = sech_WR(x)
    val = 1/cosh(x);
end

function val = a1long(L,b,p1,p2,q0)
    if( L/b > 10.0)
        C = 3.06/power((L/b),0.44);
     else 
        C = 1.0;
    end

    C1 = 1.22;
    C2 = 0.4288;
    C3 = -1.651;
    C4 = 1.523;
    C5 = 0.1477;
    miu = 0.585;

    Rg2 = Rgsquare(L,b);
    Rg22 = Rg2*Rg2;
    Rg = sqrt(Rg2);
    Rgb = Rg*q0/b;

    b2 = b*b;
    b3 = b*b*b;
    b4 = b3*b;
    q02 = q0*q0;
    q03 = q0*q0*q0;
    q04 = q03*q0;
    q05 = q04*q0;

    Rg02 = Rg2*q02;

    t1 = (b*C*((4.0/15.0 - exp((-(Rg02/b2))) *((11.0/15.0 +...
        (7.0*b2)/(15.0*Rg02))) +(7.0*b2)/(15.0*Rg02))));

    t2 = (2.0*b4*(((-1.0) + exp((-(Rg02/b2))) +Rg02/b2))*...
        ((1.0 + 0.5*(((-1.0) -tanh((-C4 + Rgb/C5)))))));

    t3 = ((C3*power(Rgb,((-3.0)/miu)) +C2*power(Rgb,((-2.0)/miu)) +...
        C1*power(Rgb,((-1.0)/miu))));

    t4 = ((1.0 + tanh(((-C4) + Rgb)/C5)));

    t5 = (1.0/(b*p1*power(q0,((-1.0) - p1 - p2)) -b*p2*power(q0,((-1.0) ...
        - p1 - p2))));

    t6 = (b*C*(((-((14.0*b3)/(15.0*q03*Rg2))) +(14.0*...
        b3*exp((-(Rg02/b2))))/(15.0*q03*Rg2) ...
        +(2.0*exp((-(Rg02/b2)))*q0*((11.0/15.0 ...
        +(7.0*b2)/(15.0*Rg02)))*Rg2)/b)));

    t7 = (Rg*((C3*power(((Rgb)),((-3.0)/miu)) +...
        C2*power(((Rgb)),((-2.0)/miu)) +...
        C1*power(((Rgb)),((-1.0)/miu))))*power(sech_WR(((-C4) ...
        +Rgb)/C5),2));

    t8 = (b4*Rg*(((-1.0) + exp((-(Rg02/b2)))...
        +Rg02/b2))*power(sech_WR(((-C4) + Rgb)/C5),2));

    t9 = (2.0*b4*(((2.0*q0*Rg2)/b -(2.0*exp((-(Rg02/b2)))...
        *q0*Rg2)/b))*((1.0 + 0.5*(((-1.0) -tanh(((-C4) + Rgb)/C5))))));

    t10 = (8.0*b4*b*(((-1.0) + exp((-(Rg02/b2))) + Rg02/b2))...
        *((1.0 + 0.5*(((-1.0) - tanh(((-C4) +Rgb)/C5))))));

    t11 = (((-((3.0*C3*Rg*power(((Rgb)),((-1.0) -3.0/miu)))/miu)) ...
        - (2.0*C2*Rg*power(((Rgb)),((-1.0) -2.0/miu)))/miu ...
        - (C1*Rg*power(((Rgb)),((-1.0) -1.0/miu)))/miu));

    t12 = ((1.0 + tanh(((-C4) + Rgb)/C5)));

    t13 = (b*C*((4.0/15.0 - exp((-(Rg02/b2)))...
        *((11.0/15.0 +(7.0*b2)/(15.0*q02* Rg2))) +(7.0*b2)/(15.0*Rg02))));

    t14 = (2.0*b4*(((-1.0) + exp((-(Rg02/b2))) +...
        Rg02/b2))*((1.0 + 0.5*(((-1.0) - tanh(((-C4) +Rgb)/C5))))));

    t15 = ((C3*power(((Rgb)),((-3.0)/miu)) +C2*power(((Rgb)),((-2.0)/miu))...
        +C1*power(((Rgb)),((-1.0)/miu))));


    val = (power(q0,p1)*(((-((b*pi)/(L*q0))) +t1/L +t2/(q04*Rg22) ...
        +0.5*t3*t4)) + (t5*((power(q0,(p1 - p2))*(((-power(q0,(-p1)))...
        *(((b2*pi)/(L*q02) +t6/L +t7/(2.0*C5) -t8/(C5*q04*Rg22) ...
        + t9/(q04*Rg22) -t10/(q05*Rg22) + 0.5*t11*t12)) ...
        -b*p1*power(q0,((-1.0) - p1))*(((-((b*pi)/(L*q0))) + ...
        t13/L +t14/(q04*Rg22) + 0.5*t15*((1.0 + tanh(((-C4) +Rgb)/C5...
        )))))))))));

end

function val = a2long(L,b,p1,p2,q0)
    if L/b > 10.0
        C = 3.06/power((L/b),0.44);
    else 
        C = 1.0;
    end

    C1 = 1.22;
    C2 = 0.4288;
    C3 = -1.651;
    C4 = 1.523;
    C5 = 0.1477;
    miu = 0.585;

    Rg2 = Rgsquare(L,b);
    Rg22 = Rg2*Rg2;
    b2 = b*b;
    b3 = b*b*b;
    b4 = b3*b;
    q02 = q0*q0;
    q03 = q0*q0*q0;
    q04 = q03*q0;
    q05 = q04*q0;
    Rg = sqrt(Rg2);
    Rgb = Rg*q0/b;
    Rg02 = Rg2*q02;
    
    t1 = (1.0/(b* p1*power(q0,((-1.0) - p1 - p2)) -b*p2...
        *power(q0,((-1.0) - p1 - p2)) ));

    t2 = (b*C*(((-1.0*((14.0*b3)/(15.0*q03*Rg2))) ...
        +(14.0*b3*exp((-(Rg02/b2))))/(15.0*q03*Rg2)...
        +(2.0*exp((-(Rg02/b2)))*q0*((11.0/15.0 +(7*b2)...
        /(15.0*Rg02)))*Rg2)/b)))/L;

    t3 = (Rg*((C3*power(((Rgb)),((-3.0)/miu)) +C2*power(((Rgb)),((-2.0)/miu))...
        +C1*power(((Rgb)),((-1.0)/miu))))* power(sech_WR(((-C4) ...
        +Rgb)/C5),2.0))/(2.0*C5);

    t4 = (b4*Rg*(((-1.0) + exp((-(Rg02/b2))) +Rg02/b2))...
        *power(sech_WR(((-C4) +Rgb)/C5),2))/(C5*q04*Rg22);

    t5 = (2.0*b4*(((2.0*q0*Rg2)/b -(2.0*exp((-(Rg02/b2)))...
        *q0*Rg2)/b))*((1.0 + 0.5*(((-1.0) - tanh(((-C4) +Rgb)/C5))))))...
        /(q04*Rg22);

    t6 = (8.0*b4*b*(((-1.0) + exp((-(Rg02/b2))) +Rg02/b2))...
        *((1.0 + 0.5*(((-1) - tanh(((-C4) +Rgb)/C5))))))/(q05*Rg22);

    t7 = (((-((3.0*C3*Rg*power(((Rgb)),((-1.0) -3.0/miu)))/miu)) ...
        - (2.0*C2*Rg*power(((Rgb)),((-1.0) - 2.0/miu)))/miu -...
        (C1*Rg*power(((Rgb)),((-1.0) - 1.0/miu)))/miu));

    t8 = ((1.0 + tanh(((-C4) + Rgb)/C5)));

    t9 = (b*C*((4.0/15.0 - exp((-(Rg02/b2)))*((11.0/15.0 ...
        +(7.0*b2)/(15*Rg02))) + (7.0*b2)/(15.0*Rg02))))/L;

    t10 = (2.0*b4*(((-1) + exp((-(Rg02/b2))) +Rg02/b2))...
        *((1.0 + 0.5*(((-1) - tanh(((-C4) +Rgb)/C5))))))/(q04*Rg22);

    val = ((-1.0*(t1* ((-power(q0,-p1)*(((b2*pi)/(L*q02) +t2 + t3 ...
        - t4 + t5 - t6 + 0.5*t7*t8)) - b*p1*power(q0,((-1.0) - p1))...
        *(((-((b*pi)/(L*q0))) + t9 + t10 + 0.5*((C3*power(((Rgb)),...
        ((-3.0)/miu)) + C2*power(((Rgb)),((-2.0)/miu)) +C1*power(((Rgb)),...
        ((-1.0)/miu))))*((1.0 + tanh(((-C4) + Rgb)/C5))))))))));
end

function val = ashort(L,b,p1short,p2short,factor,pdiff,q0)
    Rg2_sh = Rgsquareshort(L,b);
    Rg2_sh2 = Rg2_sh*Rg2_sh;
    b3 = b*b*b;
    t1 = ((q0*q0*Rg2_sh)/(b*b));
    Et1 = exp(t1);
    Emt1 = 1.0/Et1;
    q02 = q0*q0;
    q0p = power(q0,(-4.0 + p1short));

    val = ((factor/(L*(pdiff)*Rg2_sh2)*((b*Emt1*q0p*((8.0*b3*L - ...
        8.0*b3*Et1*L - 2.0*b3*L*p2short +2.0*b3*Et1*L*p2short + ...
        4.0*b*L*q02*Rg2_sh + 4.0*b*Et1*L*q02*Rg2_sh -2.0*b*Et1*L*...
        p2short*q02*Rg2_sh - Et1*pi*q02*q0*Rg2_sh2 +Et1*p2short*...
        pi*q02*q0*Rg2_sh2))))));
end

function val = a1short(L,b,p1short,p2short,q0)
    factor = 1.0;
    val = ashort(L, b, p1short, p2short, factor, p1short - p2short, q0);
end

function val = a2short(L,b,p1short,p2short,q0)
    factor = -1.0;
    val = ashort(L, b, p2short, p1short, factor, p1short-p2short, q0);
end

function val = w_WR(x)
    val = 0.5*(1 + tanh((x - 1.523)/0.1477));
end

function val = u1(q,L,b)
    val = Rgsquareshort(L,b)*q*q;
end

function val = u_WR(q,L,b)
    val = Rgsquare(L,b)*q*q;
end

function val = Sdebye_kernel(arg)
    val = 2.0*(exp(-arg) + arg -1.0)/(arg*arg);
end

function val  = Sdebye(q,L,b)
    arg = u_WR(q,L,b);
    val = Sdebye_kernel(arg);
end

function val = Sdebye1(q,L,b)
    arg = u1(q,L,b);
    val = Sdebye_kernel(arg);
end

function val = Sexv(q,L,b)
    C1=1.22;
    C2=0.4288;
    C3=-1.651;
    miu = 0.585;
    qRg = q*sqrt(Rgsquare(L,b));
    x = power(qRg, -1.0/miu);

    val = (1.0 - w_WR(qRg))*Sdebye(q,L,b) + w_WR(qRg)*x*(C1 + x*(C2 + x*C3));
end

function val = Sexvnew(q,L,b)
    C1 =1.22;
    C2 =0.4288;
    C3 =-1.651;
    miu = 0.585;
    del=1.05;
    qRg = q*sqrt(Rgsquare(L,b));
    x = power(qRg, -1.0/miu);

    %%calculating the derivative to decide on the corection (cutoff) term?
    %% I have modified this from WRs original code
    qdel = (Sexv(q*del,L,b)-Sexv(q,L,b))/(q*del - q);
    if qdel >= 0
       C_star2 = 0.0; 
    else
       C_star2 = 1.0;
    end
    val = (1.0 - w_WR(qRg))*Sdebye(q,L,b) +C_star2*w_WR(qRg)*...
        x*(C1 + x*(C2 + x*C3));
end

function val = Sk_WR(q,L,b)
    p1 = 4.12;
    p2 = 4.42;
    p1short = 5.36;
    p2short = 5.62;
    q0 = 3.1;
    
    if (1.9/sqrt(Rgsquareshort(L,b))) > 3.0
        q0short = (1.9/sqrt(Rgsquareshort(L,b)));
    else
        q0short = 3.0;
    end
    
    if L/b > 10.0
        C= 3.06/((L/b)^0.44);
    else
        C = 1.0;
    end
    
    val = zeros(size(q));

    if L > 4*b % Longer Chains
        for i=1:length(q)
            if (q(i)*b <= 3.1) 
                Sexvmodify = Sexvnew(q(i), L, b);
                val(i) = Sexvmodify + C * (4.0/15.0 + 7.0/(15.0*u_WR(q(i),L,b)) - ...
                    (11.0/15.0 + 7.0/(15.0*u_WR(q(i),L,b)))*...
                    exp(-u_WR(q(i),L,b)))*(b/L);
            else % q(i)*b > 3.1
                val(i) = a1long(L, b, p1, p2, q0)/(power((q(i)*b),p1)) + ...
                    a2long(L, b, p1, p2, q0)/(power((q(i)*b),p2)) + pi/(q(i)*L);
            end
        end
    else %L <= 4*b Shorter Chains
        for i=1:length(q)
            if (q(i)*b <= q0short ) 
                if (q(i)*b<=0.01) 
                    val(i) = 1.0 - Rgsquareshort(L,b)*(q(i)*q(i))/3.0;
                else
                    val(i) = Sdebye1(q(i),L,b);
                end
            else %  //q*b > max(1.9/sqrt(Rgsquareshort(q(i),L,b)),3)
                val(i) = a1short(L,b,p1short,p2short,q0short)/(power((q(i)*b),p1short))...
                    +a2short(L,b,p1short,p2short,q0short)/(power((q(i)*b),p2short))...
                    +pi/(q(i)*L);
            end
        end
    end
end
