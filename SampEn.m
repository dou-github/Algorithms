function SAMPEN = SampEn(data,dim, r)

N = length(data);

for m = dim:dim+1
    X = [];
    B = [];
    
    for i = 1:m
        X(i,:) = data(i:N-m+i);
    end
    
    for i = 1:N-m+1
        dist = max(abs(X - repmat(X(:,i),1, N-m+1)));
        count = any((dist <= r),1);
        B(i) = (sum(count)-1)/(N-m);
    end

    res(m-dim+1) = sum(B)/(N-m+1);
end

SAMPEN = -log(res(2)/res(1));
