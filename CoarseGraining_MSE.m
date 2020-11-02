% Coarse Graining
function output = CoarseGraining_MSE(data, scale);

N = length(data);
r = fix(N/scale);

for i = 1:r
    Y(:,i) = data((i-1)*scale+1:i*scale);
end

output = mean(Y,1);