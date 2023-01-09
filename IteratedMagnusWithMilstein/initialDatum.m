function phi=initialDatum(x,v)
%%INITIALDATUM computes the initial datum of a call option (x-v)^+ and
% computes an invertible output (x-v)^+ + lambda.* I
%   Input:
%       x (nx x 1 array): position
%       v (1 x nv array): velocity
%   Output:
%       phi (nx x nv array): (x-v)^+
%     phi = (v-x);
    phi = (v+x);
%     phi=eye(length(x));
    phi(phi<=0)=0;
end