function gaufilter = gaufilter(M, BT, S)

gaufilter = zeros(1, 2*M*S+1);
for n = -M*S:M*S
    idx = n+M*S+1;
    gaufilter(idx) = exp(-(2*pi*pi/log(2)) * (BT/M)^2 * n^2);
end

gain = sum(gaufilter);
gaufilter = gaufilter/gain;

end