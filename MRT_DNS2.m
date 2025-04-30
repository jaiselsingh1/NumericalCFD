function cylinderFlow_LBM_fullCircle_MRT_DNS_parallel()
  % Parallel MRT–LBM DNS around a full cylinder
  % D2Q9, 2× mesh refinement, full channel height,
  % inlet left, free-stream top/bottom, zero-grad outflow right,
  % skips collision inside the cylinder for speed,
  % draws velocity every 10 iterations.

  %% PARAMETERS
  factor  = 2;               % mesh refinement
  Nx_phys = 200; Ny_phys = 100;
  Nx = factor * Nx_phys;     Ny = factor * Ny_phys;
  cx0 = factor*(Nx_phys/3);  cy0 = factor*(Ny_phys/2);
  R   = factor * 20;         % cylinder radius

  U   = 0.05;                % free‐stream velocity
  Re  = 1000;                % Reynolds number
  nu  = U * 2*R / Re;        % kinematic viscosity
  tau = 3*nu + 0.5;          % relaxation time
  omega = 1/tau;             % collision frequency

  Nt  = 20000;               % timesteps

  %% D2Q9 LATTICE
  cx  = [0  1  0 -1  0  1 -1 -1  1];
  cy  = [0  0  1  0 -1  1  1 -1 -1];
  w   = [4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36];
  opp = [1  4  5  2  3  8  9  6  7];

  %% MRT MOMENT MATRICES
  M = [ ...
     1,  1,  1,  1,  1,  1,  1,  1,  1;
    -4, -1, -1, -1, -1,  2,  2,  2,  2;
     4, -2, -2, -2, -2,  1,  1,  1,  1;
     0,  1,  0, -1,  0,  1, -1, -1,  1;
     0, -2,  0,  2,  0,  1, -1, -1,  1;
     0,  0,  1,  0, -1,  1,  1, -1, -1;
     0,  0, -2,  0,  2,  1,  1, -1, -1;
     0,  1, -1,  1, -1,  0,  0,  0,  0;
     0,  0,  0,  0,  0,  1, -1,  1, -1 ];
  Minv = inv(M);

  % relaxation rates: conserve density & momentum, shear modes=omega
  S = diag([0,1.1,1.1,0,1.1,0,1.1,omega,omega]);

  %% INITIALIZE FIELDS
  rho = ones(Ny, Nx);
  u   = zeros(Ny, Nx);
  v   = zeros(Ny, Nx);
  f   = zeros(9, Ny, Nx);
  feq = zeros(9, Ny, Nx);

  for k = 1:9
    f(k,:,:) = w(k) * rho;  % start at equilibrium (rest)
  end

  % solid / fluid masks
  [X,Y] = meshgrid(1:Nx, 1:Ny);
  solid = (X - cx0).^2 + (Y - cy0).^2 <= R^2;
  fluid = ~solid;
  L     = Ny * Nx;
  idx   = find(fluid(:));      % linear indices of fluid nodes

  % start parallel pool if needed
  if isempty(gcp('nocreate'))
    parpool;
  end

  %% MAIN LOOP
  for t = 1:Nt
    %------------ Macroscopic ------------
    rho = squeeze(sum(f,1));
    um  = squeeze(f(2,:,:) + f(6,:,:) + f(9,:,:) ...
               - f(4,:,:) - f(7,:,:) - f(8,:,:));
    vm  = squeeze(f(3,:,:) + f(6,:,:) + f(7,:,:) ...
               - f(5,:,:) - f(8,:,:) - f(9,:,:));
    u(:) = 0;  v(:) = 0;
    u(fluid) = um(fluid) ./ rho(fluid);
    v(fluid) = vm(fluid) ./ rho(fluid);

    %------------ Equilibrium ------------
    for k = 1:9
      cu = 3*(cx(k)*u + cy(k)*v);
      feq(k,:,:) = w(k) * rho .* (1 + cu + 0.5*cu.^2 - 1.5*(u.^2+v.^2));
    end

    %------------ MRT Collision on Fluid ------------
    fvec = reshape(f, [9, L]);
    cold = fvec(:, idx);      % only fluid nodes
    m    = M * cold;

    rh = rho(fluid);
    uu = u(fluid);
    vv = v(fluid);
    tmp = rh .* (uu.^2 + vv.^2);
    meq = zeros(9, numel(rh));
    meq(1,:) = rh';
    meq(2,:) = (-2*rh + 3*tmp)';
    meq(3,:) = ( rh - 3*tmp)';
    meq(4,:) = ( rh.*uu )';
    meq(5,:) = (-rh.*uu )';
    meq(6,:) = ( rh.*vv )';
    meq(7,:) = (-rh.*vv )';
    meq(8,:) = ( rh.*(uu.^2 - vv.^2) )';
    meq(9,:) = ( rh.*(uu.*vv) )';

    m_post      = m - S * (m - meq);
    fvec(:,idx) = Minv * m_post;
    f = reshape(fvec, [9,Ny,Nx]);

    %------------ Parallel Streaming & Bounce-Back ------------
    % streaming
    ftemp = zeros(size(f));
    parfor k = 1:9
      ftemp(k,:,:) = circshift(f(k,:,:), [0, cy(k), cx(k)]);
    end
    % bounce-back
    fnew = zeros(size(f));
    parfor k = 1:9
      sl = squeeze(ftemp(k,:,:));
      sl(solid) = ftemp(opp(k), solid);
      fnew(k,:,:) = sl;
    end
    f = fnew;

    %------------ Boundary Conditions ------------
    % Inlet at x=1
    rho_in = ones(Ny,1);
    u_in   = U * rho_in;
    v_in   = zeros(Ny,1);
    for k = [2,6,9]
      cu = 3*(cx(k)*u_in + cy(k)*v_in);
      f(k,:,1) = w(k)*rho_in .* (1 + cu + 0.5*cu.^2 - 1.5*(u_in.^2+v_in.^2));
    end
    % Free-stream at y=1 and y=Ny
    for yb = [1, Ny]
      rho_b = rho(yb,:);
      u_b   = U * rho_b;
      v_b   = zeros(1, Nx);
      if yb == 1
        ks = [5,8,9];
      else
        ks = [3,6,7];
      end
      for k = ks
        cu = 3*(cx(k)*u_b + cy(k)*v_b);
        f(k,yb,:) = w(k)*rho_b .* (1 + cu + 0.5*cu.^2 - 1.5*(u_b.^2+v_b.^2));
      end
    end
    % Outflow at x=Nx
    f(:,:,Nx) = f(:,:,Nx-1);

    %------------ Visualization every 10 iterations ------------
    if mod(t,10) == 0
      vel = sqrt(u.^2 + v.^2);
      imagesc(vel);
      set(gca,'YDir','normal');
      axis equal tight; colorbar;
      title(sprintf('Iter %d / %d', t, Nt));
      drawnow;
    end
  end
end
