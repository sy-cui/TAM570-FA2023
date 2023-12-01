function [dim,Lx,Lw,LD,Lwm,LJm,LJf] = set_mono_param(
    N,...               % Vector of length 2 or 3
    fourier=0,...       % If Fourier basis is used. Defaults to all 0.
    mfac=1.5,...        % Dealiasing factor. Defaults to 1.5
    ffac=2.0            % Plotting factor. Defaults to 2.0
    );    
    
    % Compute 2D or 3D reference domain parameters
    % Accept Legendre or Fourier basis
    % Return params as individual 1D versions
    %
    % Outputs (cell arrays):
    %   dim :   Dimension (integer 2 or 3) of the problem
    %   Lx  :   List of quadrature nodes
    %   Lw  :   List of quadrature weights
    %   LD  :   List of differentiation matrices
    %   Lwm :   List of dealiasing quadrature weights
    %   LJm :   List of dealiasing interpolation matrices 
    %   LJf :   List of plotting interpolation matrices

    %% Argument parsing
    validateattributes(N,{'numeric'},{'vector','integer','positive'});
    validateattributes(mfac,{'numeric'},{'scalar','>',1.0});
    validateattributes(ffac,{'numeric'},{'scalar','>',1.0});

    dim = length(N);
    if dim ~= 2 && dim ~= 3;
        error("N must be of length 2 or 3.");
    end;

    if length(fourier) == 1;
        fourier = repmat([fourier], 1, dim);
    elseif length(fourier) ~= dim;
        error("fourier must be a scalar or a vector of equal length to N");
    end; 

    %% Dealiasing and plotting number of points
    M = ceil(mfac * N); F = ceil(ffac * N);
    
    Lx  = cell(1,dim);  % Quadrature nodes
    Lw  = cell(1,dim);  % Quadrature weights
    LD  = cell(1,dim);  % Differentiation matrix
    Lwm = cell(1,dim);  % Dealiasing quadrature weights
    LJm = cell(1,dim);  % Dealiasing interpolation matrix 
    LJf = cell(1,dim);  % Plotting interpolation matrix 

    %% Obtain 1D operators
    for i = 1:dim;
        if fourier(i) == 0; % Legendre
            [z ,w ] = zwgll(N(i));
            [zm,wm] = zwgl (M(i));
            [zf,fp] = zwuni(F(i));
            D  = new_deriv_mat(z);
            Jm = interp_mat(zm,z);
            Jf = interp_mat(zf,z);

        else;               % Fourier
            [z ,w ] = zwuni(N(i));
            [zm,wm] = zwuni(M(i));
            [zf,wf] = zwuni(F(i));
            D  = f_deriv_mat(z);
            Jm = f_interp_mat(zm,z);
            Jf = f_interp_mat(zf,z);

        end;

        Lx{i}  = z;
        Lw{i}  = w;
        LD{i}  = D;
        Lwm{i} = wm;
        LJm{i} = Jm;
        LJf{i} = Jf;

    end;

end; % set_mono_param