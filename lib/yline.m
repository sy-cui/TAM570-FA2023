function [] = yline(yval, varargin);

    plot(xlim, [yval yval], varargin{:});