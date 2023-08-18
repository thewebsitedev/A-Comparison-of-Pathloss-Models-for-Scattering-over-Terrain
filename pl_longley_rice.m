X = importdata("X.04");
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

figure;
plot(e2(:,1),e2(:,2), 'DisplayName', 'current'); hold on; % Add a display name to the plot

% Define your frequency and antenna heights
f = 960; % (Frequency in MHz): This should be in the range of 20 to 20000 MHz. For example, you might choose a value in the cellular communication frequency band, such as 900 MHz or 1800 MHz.
h1 = 50; % transmitter height in meters. (Antenna heights in m): These should be in the range of 0.5 to 3000 m. You might choose values like 30 m for a cell tower (h1) and 1.5 m for a handheld device (h2).
h2 = 2;  % receiver height in meters. (Antenna heights in m): These should be in the range of 0.5 to 3000 m. You might choose values like 30 m for a cell tower (h1) and 1.5 m for a handheld device (h2).
hh = 30; % (Terrain irregularity parameter): This parameter describes the roughness of the terrain. Based on the given text, values could be in the range of 0 to 700. For example, you might use a value of 30 for plain terrain or 300 for rugged mountains.
% clim = ; % (Climate): This parameter is not clearly defined in the given text, so you might have to look up appropriate values in the relevant literature. Generally, climate parameters might describe different atmospheric conditions (humidity, temperature, pressure, etc.) and could be enumerated types (e.g., 1 for temperate, 2 for tropical, etc.).
% pol = "vertical"; % (Polarization): This parameter could take two values - 'vertical' or 'horizontal'. In many cases, vertical polarization is used for VHF and UHF bands.

% Define your areas
areas = "longley_rice";

PL_longley_rice = zeros(N,1);

% Loop through all points and calculate path loss
for area = areas
    disp(area)
    for i = 1:N
        % Calculate the distance from the source to the point
        d = norm(R2(i,:) - source) / 1000; % Convert from m to km. (Distance in km): This should be in the range of 1 to 2000 km. You might choose a value like 10 km to start.
        PL_longley_rice(i) = longley_rice(f, d, h1, h2, hh, Lambda);
    end
end

% Calculate the shift value for each area
shift_value = e2(1,2) - PL_longley_rice(1);

% Shift the path loss values for each area
PL_longley_rice_shifted = PL_longley_rice + shift_value;

% Invert the shifted path loss values
PL_longley_rice_shifted_inverted = max(PL_longley_rice_shifted) - PL_longley_rice_shifted;

% Adjust the starting point of the path loss curves
PL_longley_rice_shifted_inverted = PL_longley_rice_shifted_inverted - PL_longley_rice_shifted_inverted(1) + e2(1,2);

validIdx = ~isnan(PL_longley_rice_shifted_inverted);

plot(R2(validIdx,1), PL_longley_rice_shifted_inverted(validIdx), 'DisplayName', 'open'); hold on;

title('Path Loss Comparison');
xlabel('Distance (m)');
ylabel('Electric Field (dB) / Path Loss (dB)');
legend('Location', 'southeast');
hold off;

% Step 1: Get the error at each point
error = e2(:,2) - PL_longley_rice_shifted_inverted(validIdx);

% Step 2: Square each of these errors
error_squared = error.^2;

% Step 3: Sum them
sum_error_squared = sum(error_squared);

% Step 4: Divide your answer by the number of points taken
mean_error_squared = sum_error_squared / length(error_squared);

% Step 5: Get the square root of the result
std_dev_error = sqrt(mean_error_squared);

% Print the standard deviation of the error
fprintf('The standard deviation of the error is %f dB\n', std_dev_error);


function s = sY(a, X,DeltaX,GrossStep)
    Temp = (double(a * DeltaX)) / GrossStep;
    Index = fix(Temp);
    Prop = Temp - Index;
    s = X(Index+1,2) + Prop * (X(Index+2,2) - X(Index+1,2));
end

function s = sX(a,DeltaX)
        s = a * DeltaX;
end

% Longley-Rice model for propagation loss over irregular terrain
function Lp = longley_rice(f, d, h1, h2, hh, Lambda)
    % Constants
    c = 3e8;   % Speed of light (m/s)
    f = f * 1e6; % Frequency conversion from MHz to Hz
    
    % Calculate effective heights
    heb = h1 - hh;
    hem = h2 - hh;
    if heb < 5
        heb = 5;
    end
    if hem < 5
        hem = 5;
    end
    
    % Calculate horizon distances
    dLSb = sqrt(17 * heb);
    dLSm = sqrt(17 * hem);
    dL = dLSb + dLSm;
    
    % Calculate horizon elevation angles
    yeb = 0.0005 * sqrt(d/dLSb * 1.3 / (1 + hh / (4*heb*dL)));
    yem = 0.0005 * sqrt(d/dLSm * 1.3 / (1 + hh / (4*hem*dL)));
    
    % Calculate v-parameters for obstacles
    vb1 = 1.2915 * yeb * ((f*c*dLSb*(d-dL)/(d-dLSm))^(1/3));
    vm1 = 1.2915 * yem * ((f*c*dLSm*(d-dL)/(d-dLSb))^(1/3));
    vb2 = 1.2915 * yeb * ((f*c*dLSb*(d-dL)/(d-dLSm))^(1/3));
    vm2 = 1.2915 * yem * ((f*c*dLSm*(d-dL)/(d-dLSb))^(1/3));
    
    % Calculate diffraction losses A1 and A2
    A1 = A(vb1) + A(vm1);
    A2 = A(vb2) + A(vm2);
    
    % Compute diBraction loss
    LD = A1 + A2; % Simplified, actual computation will be more complex
    
    % Free space loss
    Lfs = 20 * log10((4*pi*d*1e3)/Lambda);
    
    % Median transmission loss
    Lp = LD + Lfs;

end

function loss = A(v)
    if 0 <= v && v < 0.8
        loss = 20 * log10(0.5 * 0.62 * v);
    elseif 0.8 <= v && v < 1
        loss = 20 * log10(0.5 * exp(-0.95 * v));
    elseif 1 <= v && v < 2.4
        loss = 20 * log10(0.4 * sqrt(0.1184 - (0.38 - 0.1 * v)^2));
    else
        loss = 20 * log10(0.225 * v);
    end
end
