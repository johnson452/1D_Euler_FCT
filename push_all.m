function [rho,u,p,E,grid] = push_all(rho,u,p,E,grid)

%Push all quantities with the MUSCL scheme

%Update energy

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

%Update p:
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

%Update p:
p = EOS(rho,u,E,grid,"calorically_ideal");

%Calculate flux (Euler equations):
flux = [ rho.*u  ; p + rho.*u.*u ; u.*(rho.*E + p) ];

end