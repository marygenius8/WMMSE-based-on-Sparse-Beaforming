function [H] = channel_realization(params)
% Generate channel matrix
%%
%INPUT:
% params.K = 30;         % =# users in each cell
% params.L_MB = 1;       % =# Macro base station in each cell
% params.L_pB = 3;       % =# Pico cell in each cell
% params.M_MB = 4;       % =# antennas at each Macro base station
% params.M_pB = 2;       % =# antennas at each Pico-cell
% params.N_user = 2;     % =# antennas at each user
% params.R = 400/sin(pi/3);    % distance between cell 0.8km, radius of each cell 400m

%OUTPUT:
% H     = (K*N,LM = params.L_MB * M_MB + params.L_pB * M_pB) channel matrix
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = params.L_MB + params.L_pB; 
idx_LpB = params.L_MB * params.M_MB;
%%%%%%%%%%%%%%%Function: hexagon_user%%%%%%%%%%%%%%%%%%%%%%%
%% generate user position in a hexagon cell
% input: K = # of users, R = radius of cell, center position of the cell
function [U_position] = hexagon_user(K,R,center)
% Generate uniform points in the simplex
% convex combination of 3points in R^3: (1,0,0) (0,1,0) (0,0,1)
    m = 3;
    X = rand(m-1,K) .^ (1./(m-1:-1:1)'); % use bsxfun(@power,...) for old release
    X = cumprod([ones(1,K);X]).*[ones(m,K)-[X;zeros(1,K)]];
    % use X as a barycentric of the triangle (0,1,z6) in the complex plane
    % so point Z is uniform in this triangle
    z6 = exp(2i*pi/6);
    Z = [0, 1, z6]*X;
    % multiply by random 6th-roots of 1 to map into unit hexagonal
    Z = Z .* (z6.^floor(6*rand(1,K)));
    % U_position=params.MB_pos+[R*real(Z);R*imag(Z)];  %% user positions
    % BS_position=[params.MB_pos pB_position1 pB_position2 pB_position3];  %%pico-cell positions 
    U_position=center+[R*real(Z);R*imag(Z)];  %% user positions
end

%%%%%%%%%%%%%%%Function: hexagon_cell%%%%%%%%%%%%%%%%%%%%%%%
%% generate hexagon cell BS and pico cell or RRH positions, 
%% and corresponding user positions in each hexagon cell
% current version--create 7 hexagon cells (BS and user)next to each other
% in next version, draw the hexgon cell and user illustration can be inculded
% input: wireless network scenario setup
function [BS_pos, U_pos] = hexagon_cell(params)
% % n = params.CellN-1;                          % how many hexagon cells next to the center hexagon cell
% % theta represents the direction of surounding hexagon centers
% theta = linspace(pi/2,5*pi/2,n+1);    
% % alpha represent the direction of 3 uniformly distributed pico_cell
% alpha = linspace(0,2*pi,params.L_pB+1);   
% % compute the center position of other hexagons next to the center hexagon
% x0=sqrt(3)*R*cos(theta(1:n));
% y0=sqrt(3)*R*sin(theta(1:n));
% center = [x0;y0];
% center = [params.center center];    % permute all the Macro BS position
% BS_pos = [ ];
% % 3 pico cell position in each hexagon
% for i=1:params.CellN
%     pB_positionx=center(1,i)+R/2*cos(alpha(1:params.L_pB));
%     pB_positiony=center(2,i)+R/2*sin(alpha(1:params.L_pB));
%     pB_position = [pB_positionx; pB_positiony];
%     BS_pos = [BS_pos center(:,i) pB_position];
% end

% alpha represent the direction of 3 uniformly distributed pico_cell inside a hexagon
    alpha = linspace(0,2*pi,params.L_pB+1);
    % 3 pico cell position in each hexagon+1BS
%         pB_positionx=params.center(1,:)+params.R/2*cos(alpha(1:params.L_pB));
%         pB_positiony=params.center(2,:)+params.R/2*sin(alpha(1:params.L_pB));
%         pB_position = [pB_positionx; pB_positiony];
    pB_position=params.center+params.R/2*[cos(alpha(1:params.L_pB));sin(alpha(1:params.L_pB))];
    BS_pos = [params.center pB_position];
    U_pos = hexagon_user(params.K, params.R, params.center);
end

[BS_position, U_position] = hexagon_cell(params);
 %% %%%Generate Large-Scale Fading%%%%%%%%%%% %%
 for l=1:L
    for k=1:params.K
        if l <= params.L_MB
            d=(norm(BS_position(:,l)-U_position(:,k)));
            D(k,l)=10^(-0.015)*10^(normrnd(0,6.3096)/20)/(d^(1.88));
        else
            d=(norm(BS_position(:,l)-U_position(:,k)));
            D(k,l)=10^(-0.78)*10^(normrnd(0,6.3096)/20)/(d^(1.835));
        end
    end
end

%%%%%%Generate Small-Scale Fading%%%%%%%%%%%%%
U_BS_norm_r = normrnd(0,1/sqrt(2),params.N_user,params.M_MB);
U_BS_norm_c = normrnd(0,1/sqrt(2),params.N_user,params.M_MB);
U_pB_norm_r = normrnd(0,1/sqrt(2),params.N_user,params.M_pB);
U_pB_norm_c = normrnd(0,1/sqrt(2),params.N_user,params.M_pB);
for l=1:L    
    if l<=params.L_MB
        idx_cln = params.M_MB*(l-1)+(1:params.M_MB);
%         idx_LpB = l * params.M_MB;
    else
        idx_cln = idx_LpB+(1:params.M_pB);
        idx_LpB = idx_LpB+params.M_pB;
    end
    for k=1:params.K
%         idx_row_start = params.N_user*(k-1)+1;
%         idx_row_end = params.N_user*k;
        idx_row = params.N_user*(k-1)+(1:params.N_user);
        if l<=params.L_MB
%             H(idx_row_start:idx_row_end,idx_cln) = D(k,l)*(U_BS_norm_r+1i*U_BS_norm_c); 
            H(idx_row,idx_cln) = D(k,l)*(U_BS_norm_r+1i*U_BS_norm_c); 
        else
%             H(idx_row_start:idx_row_end,idx_cln) = D(k,l)*(U_pB_norm_r+1i*U_pB_norm_c);
            H(idx_row,idx_cln) = D(k,l)*(U_pB_norm_r+1i*U_pB_norm_c);
        end
    end
end
end
