function [rho,u,p,E,grid] = push_all(rho,u,p,E,grid)

%Push all quantities with the FCT described:
%https://pdf.sciencedirectassets.com/272570/1-s2.0-S0021999100X02859/1-s2.0-0021999179900512/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEO7%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIQDd8BdayoozwKUoM1CmwTdviKp6qo9Ia3h%2B1FiHFoEZjAIgMWTQLrWmaN8vWMbBPzvFbJrr29hxGP71T0I18rP9xYcqswUINhAFGgwwNTkwMDM1NDY4NjUiDGztTMnwtLBRAPIrdiqQBWNxqti6BKGTVdk3KpaQDWlrywAxn4Y1t%2FPjz9G%2BjjQ9qFpCyQkSryKD6coHMmTHjhAn795ND6AJdiO%2FeFt7RcWk77qP8LtqUtw8YNl03K78cXpaxa7n4emZhjNwBo7jgk7%2BUx6WcGgtz21rpDQ2uR%2FCSmunBzvxMfJwT2VKOVtcEbQ3RErWaEJEYxMrt%2BKsQU2fdDO2dPev%2BPt9khy%2FB0pTCUZ%2BpK1ZVdDUO6nc41M11LWKypeAKjIKmdyEeADL17xHeA8cWDG7f%2BSfIwDXv9491WKcEZs5NkdoLEsTv5dRVicgsrTNWSmESmsgxumg5JwlfWofYs5BtJaKDkGXitiOldeKfKFL1Tm%2Bmr1mTB3H%2BehapZu2LsUZm62zFnZUDifQDrkhWStBK9jCyPzCXFcei46gozWlf%2Bbf7Q%2FS7EZ7CkN3mLjhPGA0AyUey7ZWX6vT0tsq2W3Cozv1q1SRtb4MuGvh8ONo1C1cWeuKuG88uVSuxFMVXqKXsee0qKi5fqmNl2VukJnbeetv6NhUSdRfX5qG5AvD2%2BYAtu2i6l9xxE0hxsdQWbyyw25eqm9dDadx5CtFREoWhJO9EJ5tsUgL039s6mQ6ISzA5x6r%2BUVIRl38u8YG4Vb%2FxctTM6jTPXy2b3H3n4ytK5JT%2F4bKSrZEgQkx569kFgyaxYsEfEqewDIWyG7sgLgU%2BtBY1e3hU8dgNEwmYEAQb0P3wtQiowtrz1KEcTnWO6VPs0IzX5Nc8pRDTOXL2YMDh5fes4XEpilagy2XaG%2FSN%2BFFo2HCamTBYNDCgauySdTuokJgzZ0UzfGf8QpBlwBnNRjVWqh4ePsd3AhQ73pO%2FVFt5NZbyKpM%2BCG3UZXV9SL77XBdYvANMNSe%2BaMGOrEB%2B9j9Yp39Yy7E3R8K9cV2ECj04VPkTs9gD4a%2FeOD7HthXoR4JgdqjE5qrIj71M%2FZ37pKV5OaxGlX6j790LUU0CfaXRvq6akN8b77yWYMk6Bk%2FF6U9PH2tZfYx%2BbevHsGfrkuR%2BHgHdpMj0Bw1IqeodlfF4DuQDPxCJCgJ5KFFbQ%2BCXl8dkmNAvn7DO%2Fdai%2F%2B7vbaqDXm%2Bfgauhp4MaI39%2B1BTgKwmDsozXBPrKGkChubk&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230605T221714Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYX6UT6RNP%2F20230605%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=58e694cc967a550e40e4f6ce95030561576e66daaa1b80a263f620f6f0ed4e74&hash=0fa4626b4a6ef5040dee43c261ea3ac495898f88afc2b248f4f5b1725b9dc324&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=0021999179900512&tid=spdf-6248abf5-0a3f-43c2-950a-232e37cf6e4c&sid=da0f7889150d13473f2aef1-529c91f6051cgxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=131559025c0107570755&rr=7d2bc8d61e3630b8&cc=us

%Grab U, F vectors
U = [ rho ; rho.*u ; rho.*E ];
%F = [ rho.*u  ; p + rho.*u.*u ; u.*(rho.*E + p) ];
c = (grid.dt/grid.dx);

%%%%
% FCT:
% (1) Compute F with with a lower-order monotonic scheme
[FL_i_minus_half, FL_i_plus_half] = FL(U,grid);

% (2) Compute F with a higher-order scheme
[FH_i_minus_half, FH_i_plus_half] = FH(U,grid);

% (3) Compute "anti-diffusive flux":
A_i_plus_half = FH_i_plus_half - FL_i_plus_half;
A_i_minus_half = FH_i_minus_half - FL_i_minus_half;

% (4) Compute updated lower order solution (td= transported and diffused)
U_td = U - c*(FL_i_plus_half - FL_i_minus_half);

% (5) Limit A_i_plus_half such that the new step is free of extrema
[C_i_minus_half,C_i_plus_half] = C(A_i_plus_half,A_i_minus_half,U,U_td,grid);
A_C_i_plus_half = C_i_plus_half.*A_i_plus_half;
A_C_i_minus_half = C_i_minus_half.*A_i_minus_half;

% (6) Update the solution with the corrected fluxes
U = U_td - c*(A_C_i_plus_half - A_C_i_minus_half);
%%%%


%grab variables out of U
rho = U(1,:);
u = U(2,:)./rho;
E = U(3,:)./rho;

%Update p (Euler equations):
p = EOS(rho,u,E,grid,"calorically_ideal");

end


% FCT Functions
function [c_val_minus,c_val_plus] = C(A_i_plus_half,A_i_minus_half,U,U_td,grid)

%Compute Pi+, Qi+, and Ri+
Zero = zeros(size(U));
P_plus = max(Zero,A_i_minus_half) - min(Zero,A_i_plus_half);
Q_plus = (max_w(U,U_td,grid) - U_td)*(grid.dx/grid.dt);
R_plus = zeros(size(U));
for i = 1:max(size(U))
    for j = 1:3
        if P_plus(j,i) > 0
            R_plus(j,i) = min(1,Q_plus(j,i)/P_plus(j,i));
        elseif P_plus(j,i) == 0
            R_plus(j,i) = 0;
        end
    end
end

%Compute Pi-, Qi-, and Ri-
P_minus = max(Zero,A_i_plus_half) - min(Zero,A_i_minus_half);
Q_minus = (U_td - min_w(U,U_td,grid))*(grid.dx/grid.dt);
R_minus = zeros(size(U));
for i = 1:max(size(U))
    for j = 1:3
        if P_minus(j,i) > 0
            R_minus(j,i) = min(1,Q_minus(j,i)/P_minus(j,i));
        elseif P_minus(j,i) == 0
            R_minus(j,i) = 0;
        end
    end
end


%Compute C+ value
c_val_plus = zeros(1,grid.Nx);
for i = 1:max(size(U))
    for j = 1:3
        if A_i_plus_half(j,i) >= 0
        c_val_plus(i) = min( R_plus(j,grid.R(i)), R_minus(j,i) );
        elseif A_i_plus_half(j,i) < 0
        c_val_plus(i) = min( R_plus(j,i),R_minus(j,grid.R(i)) );
        end
    end
end

%Compute C- value
c_val_minus = zeros(1,grid.Nx);
for i = 1:max(size(U))
    for j = 1:3
        if A_i_minus_half(j,i) >= 0
        c_val_minus(i) = min( R_plus(j,i), R_minus(j,grid.L(i)) );
        elseif A_i_minus_half(j,i) < 0
        c_val_minus(i) = min( R_plus(j,grid.L(i)),R_minus(j,i) );
        end
    end
end

end


% max w
function [w_max] = max_w(w,w_td,grid)

%Compute improved fluxes
w_td = max(w,w_td);

%Find the max
w_max = max( w_td(:,grid.L) ,max(w_td,w_td(:,grid.R)) );

end

% min w
function [w_min] = min_w(w,w_td,grid)

%Compute improved fluxes
w_td = min(w,w_td);

%Find the max
w_min = min( w_td(:,grid.L) ,min(w_td,w_td(:,grid.R)) );

end

% Lower Order Flux
function [F_low_order_minus,F_low_order_plus] = FL(U,grid)

%Scheme: Lax-Friedrichs
UL = U(:,grid.L);
UR = U(:,grid.R);
F_low_order_minus = (1/2)*( Flux(UL,grid) + Flux(U,grid) ) ...
    - (grid.dx/(2*grid.dt))*(U - UL);
F_low_order_plus = (1/2)*( Flux(U,grid) + Flux(UR,grid) ) ...
    - (grid.dx/(2*grid.dt))*(UR - U);

end

% Higher Order Flux
function [F_high_order_minus,F_high_order_plus] = FH(U,grid)

%Scheme: Lax-Wendroff
% Ref: https://en.wikipedia.org/wiki/Lax–Wendroff_method

% Richtmyer 2-step Lax-Wendroff
c_half = (grid.dt/(2*grid.dx));

%Step 1:
UL = U(:,grid.L);
UR = U(:,grid.R);
U_n_half_ip = 0.5*(UR + U) - c_half*(Flux(UR,grid) - Flux(U,grid));
U_n_half_im = 0.5*(U + UL) - c_half*(Flux(U,grid) - Flux(UL,grid));

%Calculate the fluxes:
F_high_order_minus = Flux(U_n_half_im,grid);
F_high_order_plus = Flux(U_n_half_ip,grid);

end


%Flux
function [flux] = Flux(U,grid)

%grab variables out of U
rho = U(1,:);
u = U(2,:)./rho;
E = U(3,:)./rho;

%Update p (Euler equations):
p = EOS(rho,u,E,grid,"calorically_ideal");

%Calculate flux (Euler equations):
flux = [ rho.*u  ; p + rho.*u.*u ; u.*(rho.*E + p) ];

end