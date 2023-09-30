function [] = xline(xval, varargin);

    plot([xval xval], ylim, varargin{:});
