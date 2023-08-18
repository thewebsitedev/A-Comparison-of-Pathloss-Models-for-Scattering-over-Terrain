X = importdata("X.04");
R = zeros(200,2);
Epsilon_0 = 8.854e-12; 
Mu_0 = 12.56637061e-7;
GrossStep = 10.0;
f = 970e6;
c=(1.0/sqrt(Mu_0*Epsilon_0));
Lambda = c/f;
Omega = 2.0*pi*f;
Beta_0 = Omega*(sqrt(Mu_0*Epsilon_0));
Eta_0 = sqrt(Mu_0/Epsilon_0);
DeltaX = Lambda/4.0;
GrossNoSteps = 384; %70 - 250 - 384
source = [0.0, 442.0];
% E=zeros(0);
len=GrossStep*GrossNoSteps;
N = fix(len/DeltaX);
disp(N);
J1=zeros(N,2);
E=zeros(N,1);
% Split Z into two parts
N1 = floor(N/2);
N2 = ceil(N/2);
Z11 = zeros(N1, N1);
Z12 = zeros(N1, N2);
Z21 = zeros(N2, N1);
Z22 = zeros(N2, N2);
J=zeros(N,1);

% To calculate points in between
for i = 1:49697
    R(i,1)=sX(i,DeltaX);
    R(i,2)= sY(i,X,DeltaX,GrossStep);
end

R = sortrows(R, 1);

% To calculate value of E
for i=1:N
    E(i) = besselh(0, 2, Beta_0 * norm(R(i, :) -source));
end

% Compute Z11, Z12, Z21, and Z22
for i=1:N
    if i <= N1
        Z11(i, i) = DeltaX*Beta_0*Eta_0/4*(1-2*log(1.781*DeltaX*Beta_0/4/exp(1))/pi*1i);
        for j=1:i-1
            Z11(i, j) = DeltaX*Beta_0*Eta_0/4*besselh(0, 2, Beta_0 * norm(R(i,:)-R(j,:)));
            Z11(j, i) = Z11(i, j);
        end
        for j=N1+1:N
            Z12(i, j-N1) = DeltaX*Beta_0*Eta_0/4*besselh(0, 2, Beta_0 * norm(R(i,:)-R(j,:)));
        end
    else
        index = i - N1;
        Z22(index, index) = DeltaX*Beta_0*Eta_0/4*(1-2*log(1.781*DeltaX*Beta_0/4/exp(1))/pi*1i);
        for j=1:index-1
            Z22(index, j) = DeltaX*Beta_0*Eta_0/4*besselh(0, 2, Beta_0 * norm(R(i,:)-R(j+N1,:)));
            Z22(j, index) = Z22(index, j);
        end
        for j=1:N1
            Z21(index, j) = DeltaX*Beta_0*Eta_0/4*besselh(0, 2, Beta_0 * norm(R(i,:)-R(j,:)));
        end
    end
    disp(i);
end

% Concatenate Z11, Z12, Z21, and Z22
% Z = [Z11 Z12; Z21 Z22];

% Split E into E1 and E2
E1 = E(1:N1);
E2 = E(N1+1:end);

% Calculate J1 and J2
J1 = zeros(N1, 1);
J2 = zeros(N2, 1);

% To calculate value of J1
J1(1) = E1(1)/Z11(1,1);
for i= 2:N1
    sigma = 0;
    for j=1:i-1
        sigma= sigma+J1(j)*Z11(j,i);
    end
    J1(i) = (E1(i)-sigma)/Z11(i,i);
end

% To calculate value of J2
J2(1) = (E2(1)-Z21(1,:)*J1)/Z22(1,1);
for i= 2:N2
    sigma = 0;
    for j=1:i-1
        sigma= sigma+J2(j)*Z22(j,i);
    end
    J2(i) = (E2(i)-Z21(i,:)*J1-sigma)/Z22(i,i);
end

J = [J1; J2];

% To calculate value of J
% J(1) = E(1)/Z(1,1);
% for i= 2:N
%     sigma = 0;
%     for j=1:i-1
%         sigma= sigma+J(j)*Z(j,i);
%     end
%     J(i) = (E(i)-sigma)/Z(i,i);
% end

JJ(:,2)=abs(J);
for n=1:N
    JJ(n,1)=R(n,1);
end

% Start plotting J1 data
figure;
plot(JJ(:,1),JJ(:,2), 'DisplayName', 'current');

xlabel('Distance (m)');
ylabel('Electric Current (dB)');
legend('Location', 'southeast');

function s = sY(a, X,DeltaX,GrossStep)
    Temp = (double(a * DeltaX)) / GrossStep;
    Index = fix(Temp);
    Prop = Temp - Index;
    s = X(Index+1,2) + Prop * (X(Index+2,2) - X(Index+1,2));
end

function s = sX(a,DeltaX)
    s = a * DeltaX;
end
