% fractal dimension value calculation (Higuchi)
function FD_H = FD_Higuchi(x,kmax)

% input x:  vector of length N

N = length(x);
Lmk = zeros(kmax,kmax);

for k = 1:kmax
    for m = 1:k
        num = floor((N-m)/k);
        z = 0;
        for i = 1:num
            z = z + abs(x(m+i*k)-x(m+(i-1)*k));
        end
        Re = (N-1)/(num*k);
        Lmk(m,k) = (z*Re)/k;
    end
Lk(k) = sum(Lmk(:,k))/k;
end

x_i = 1./(1:kmax);
p = polyfit(log(x_i),log(Lk),1);
FD_H = p(1);
