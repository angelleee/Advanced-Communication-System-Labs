function RC = rcfun(a, M, S)
%RCFUN Generate Raised Cosine (RC) pulse
%   a : roll-off factor
%   M : oversampling factor
%   S : pulse span (in symbols)
%
%   Returns:
%       RC : RC pulse (normalized)
%       n  : sample indices

RC = zeros(1,2*M*S+1);
for n = -M*S:M*S
    idx = n+M*S+1;
    RC(idx) = sinc(n/M)*(cos(pi*a*(n/M)) / (1 - 4*(a*n/M)^2));
end

DC_gain = sum(RC);
RC = RC/DC_gain;

end
