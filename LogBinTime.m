function out=LogBinTime(tmin,tmax,bs)
%returns linear binning of time in log scale

ltb=[0 1e-8 1e-7 1e-6 1e-5 1e-4 0.001 0.010 0.100 1 10 100 1000 10000];
if length(bs)==1
    bins=bs*[1 1 1 1 1 1 1 1 1 1 1 1 1];
else
    bins=bs;
end
tb=0;

for a=1:length(ltb)-1
    dummy=linspace(ltb(a),ltb(a+1),bins(a));
    tb=[tb dummy(2:end)];
end
p=(tb>=tmin).*(tb<=tmax);
out=tb(logical(p));

