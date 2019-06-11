function [Sigscrol] = BSSCCA(Data)
A = Data(:,1:end-1);
B = Data(:,2:end);

z = [A;B];
C = cov(z.');

sx = size(A,1);
sy = size(B,1);
Cxx = C(1:sx, 1:sx) + 10^(-8)*eye(sx);
Cxy = C(1:sx, sx+1:sx+sy);
Cyx = Cxy';
Cyy = C(sx+1:sx+sy, sx+1:sx+sy) + 10^(-8)*eye(sy);
invCyy = inv(Cyy);

[Wx,r_1] = eig(inv(Cxx)*Cxy*invCyy*Cyx); % Basis in A
r_2 = sqrt(real(r_1));

V = fliplr(Wx);		% reverse order of eigenvectors
r_3 = flipud(diag(r_2));	% extract eigenvalues and reverse their order
[r,I]= sort((real(r_3)));	% sort reversed eigenvalues in ascending order
r = flipud(r);		% restore sorted eigenvalues into descending order
for j = 1:length(I)
  Wx(:,j) = V(:,I(j));  % sort reversed eigenvectors in ascending order
end
Wx = fliplr(Wx);

Wy_1 = invCyy*Cyx*Wx;     % Basis in B
Wy = Wy_1./repmat(sqrt(sum(abs(Wy_1).^2)),sy,1); % Normalize Wy

y = Wx'*A;

%[Sigscrol]=EMGscrol(Data,y,Wx);

D = pinv(Wx)';
%N = size(Data);
%Sigscrol{1}=Data;
% for i=1:N(1)
%     D(:,end-i+1:end)=0;
%     Sigscrol{i+1}=D*y+(ones(N(2)-1,1)*mean(Data'))';
% end

D(:,end-1:end) = 0;
Sigscrol = D*y;
