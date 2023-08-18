X = importdata("/Users/gautamthapar/Documents/Trinity College/Dissertation/Meeting 1/X.04");
R1 = zeros(200,2);
R2 = zeros(200,2);
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
GrossNoSteps = 384;
source = [0.0, 442.0];
% E2=zeros(0); 
len=GrossStep*GrossNoSteps;
N = fix(len/DeltaX);
E2=zeros(N,1);
e_total=zeros(N,1);
e2=zeros(N,2);
height=2.4;

% To calculate points in between
for i = 1:49697
    R1(i,1)=sX(i,DeltaX);
    R1(i,2)= sY(i,X,DeltaX,GrossStep);
    R2(i,1)=sX(i,DeltaX);
    R2(i,2)= sY(i,X,DeltaX,GrossStep)+height;
    x = sX(i,DeltaX);
end

% To calculate value of E
for i=1:N
    E2(i) = besselh(0, 2, Beta_0 * norm(R2(i, :) -source));
end 

% Calculate Z2 row by row, update e_total, and discard the row
for i = 1:N
    Z2_row = zeros(1, N);
    Z2_row(i) = besselh(0, 2, height*Beta_0);
    for j = 1:i-1
        Z2_row(j) = DeltaX*Beta_0*Eta_0/4 * besselh(0, 2, Beta_0 * norm(R2(i,:) - R1(j,:)));
        Z2_row(i) = Z2_row(i) + DeltaX*Beta_0*Eta_0/4 * besselh(0, 2, Beta_0 * norm(R2(j,:) - R1(i,:)));
    end
    Sigma = Z2_row * J;
    e_total(i) = E2(i) - Sigma;
    disp(i);
end

% To calculate field above the current
for n=1:N
   e2(n,1)=R1(n,1);
   e2(n,2)=20*log10(sqrt(abs(e_total(n)/sqrt(norm(R2(n,:)-source))).^2));
end

% Calculate square root of each error point in e2(:,2)
sqrt_error = sqrt(abs(e2(:,2)));

% Sum all of them
sum_sqrt_error = sum(sqrt_error);

% Get average
avg_sqrt_error = sum_sqrt_error / N;

% Create a new vector that contains the average square root error for each x point
avg_sqrt_error_vector = avg_sqrt_error * ones(N, 1);

% Plot the field above the current
figure;
plot(e2(:,1),e2(:,2), 'DisplayName', 'field above current'); hold on;

% Plot the average square root error as a line on the same graph
plot(e2(:,1),avg_sqrt_error_vector, 'DisplayName', 'average sqrt error'); hold on;

legend

function s = sY(a, X,DeltaX,GrossStep)
    Temp = (double(a * DeltaX)) / GrossStep;
    Index = fix(Temp);
    Prop = Temp - Index;
    s = X(Index+1,2) + Prop * (X(Index+2,2) - X(Index+1,2));
end

function s = sX(a,DeltaX)
    s = a * DeltaX;
end