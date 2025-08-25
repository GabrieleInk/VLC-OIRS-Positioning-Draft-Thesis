function p = WLS(Coord_Ancore, distanze, zk, w) % Weighted Least Squares
    
    N = size(Coord_Ancore, 1); % Number of Coord_Ancore

    % Construct matrix A and vector b
    A = zeros(N, 2);
    b = zeros(N, 1);

    for i = 1:N
        xi = Coord_Ancore(i, 1);
        yi = Coord_Ancore(i, 2);
        zi = Coord_Ancore(i, 3);
        di = distanze(i);

        % Linearized system approximation
        A(i, :) = [-2 * xi, -2 * yi]; 
        b(i) = di^2 - xi^2 - yi^2 - zi^2; % Includes z-component of Coord_Ancore
    end
    % Pesi
    if nargin < 4 || isempty(w)
        w = ones(N,1);
    end
    W = diag(w);
    % Solve for position using Weighted Least Squares: (A'WA)X = A'Wb
    gu_ = (A' * W * A) \ (A' * W * b);
    %gu_ =lsqlin((A' * W * A),(A' * W * b),[],[],[],[],zeros(3,1),5*ones(3,1));
    % Estimated position with fixed z = 0
    p = [gu_' zk]; 

end
