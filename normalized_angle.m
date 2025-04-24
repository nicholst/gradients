function x = normalized_angle( x, doac )
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% Copyright (C) - 2024 - Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
if ~exist('doac', 'var')
    doac = 1;
end
%%  Add/check optional values
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
if doac == 1
    x = 1 - squareform(pdist(x, 'cosine'));
    x = 1 - acos(x) / pi;
else
    x = 1 - squareform(pdist(x, 'cosine'));
end
% x = 1 - acos(x) / (pi);
% 1 - acos((normalizeX2(x))'*(normalizeX2(x)))/pi;

end

