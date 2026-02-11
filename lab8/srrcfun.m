function SRRC = srrcfun(a, M, S)
%SRRCFUN Generate Square-Root Raised Cosine (SRRC) pulse
%   a : roll-off factor
%   M : oversampling factor
%   S : pulse span (in symbols)
%
%   Returns:
%       SRRC : SRRC pulse (normalized)
%       n    : sample indices

SRRC = zeros(1,2*M*S+1);
for n = -M*S:M*S
    idx = n+M*S+1;
    if n ~= 0
        SRRC(idx) = (4*a/pi) * (cos((1+a)*pi*n/M) + M*sin((1-a)*pi*n/M)/(4*a*n)) / (1-(4*a*n/M)^2);
    else
        SRRC(idx) = (1-a+(4*a/pi));
    end
end

DC_gain = sum(SRRC);
SRRC = SRRC/DC_gain;

end
