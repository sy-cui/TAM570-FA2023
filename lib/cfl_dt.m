function [dt,varargout] = cfl_dt(
    CFL,        % Courant-Friedrichs-Lewy Number
    dx,         % Smallest grid size
    varargin    % Support 0, 1, or 2 variable inputs
    );
    % Compute time-step based on CFL constraint, and compute total
    % number of time steps on demand.
    %
    % Usage:
    %   dt = cfl_dt(CFL,dx)     Assume characteristic velocity U=1
    %                           and compute dt = CFL * dx
    %
    %   dt = cfl_dt(CFL,dx,u)   Compute dt = CFL * dx / u
    %
    %   [dt,nsteps] = cfl_dt(CFL,dx,Tend)
    %                           Compute dt assuming U=1, and total
    %                           number of time steps based on Tend
    %
    %   [dt,nsteps] = cfl_dt(CFL,dx,u,Tend)
    %                           Compute dt based on full CFL limit
    %                           and number of time steps 

    nin  = nargin  - 2;
    nout = nargout - 1;

    if nout ~= 0 && nout ~= 1;
        error("Unsupported number of outputs.");
    end;

    if nin ~= 0 && nin ~= 1 && nin ~= 2;
        error("Unsupported number of inputs.");
    end;

    
    if nin == 0;
        if nout == 0; 
            dt = CFL * dx; 
        else;
            error("Missing input: Tend"); 
        end;
        
    elseif nin == 1;
        if nout == 0; 
            u = varargin{1};
            dt = CFL * dx / u;
        else;
            Tend = varargin{1};
            dt = CFL * dx; 
            nsteps = ceil(Tend / dt);
            dt = Tend / nsteps;
            varargout{1} = nsteps;
        end;

    else;
        if nout == 0;
            error("Too many inputs");
        else;
            u = varargin{1};
            dt = CFL * dx / u;
            Tend = varargin{2};
            nsteps = ceil(Tend / dt);
            dt = Tend / nsteps;
            varargout{1} = nsteps;
        end;
    end;
end;