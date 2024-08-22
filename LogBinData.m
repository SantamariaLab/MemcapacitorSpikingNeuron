function out=LogBinData(p,bs)
%bins data fromcorrelator2mat

if isempty(p.data)
    out.BinData=[0 0];
    return
end


%D=highV;%p.param.Duration;
rt=p.data(:,1);
rd=p.data(:,2);%refence to 0

tb=LogBinTime(p.lowV,p.highV,bs);   

for b=1:length(tb)
    if b==1
    dd=(rt<tb(b+1));
    elseif b<length(tb)
        dd=(rt>=tb(b)).*(rt<tb(b+1));
    elseif b==length(tb)
        dd=(rt>=tb(b));
    end
    if nnz(dd)
        BinData(b)=sum(rd(logical(dd)))./nnz(dd);
    else
        BinData(b)=0;
    end
end

out.BinData=[tb' BinData'];

