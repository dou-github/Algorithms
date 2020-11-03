function APPEN = ApEn(data, dim, r)

N = length(data);

for m = dim:dim+1
    X = [];
    C = [];
    
    for i = 1:m
        X(i,:) = data(i:N-m+i);
    end
    
    for i = 1:N-m+1
        dist = abs(X - repmat(X(:,i),1,N-m+1));
        count = any((dist <= r),1);
        C(i) = sum(count)/(N-m+1);
%         count = any((dist > r),1);
%         C(i) = sum(~count)/(N-m+1);
    end
    res(m-dim+1) = sum(log(C))/(N-m+1);
end

APPEN = res(1) - res(2);
   
