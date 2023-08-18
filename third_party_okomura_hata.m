figure;
plot(e2(:,1),e2(:,2), 'DisplayName', 'current'); hold on; % Add a display name to the plot

PL_urban = zeros(N,1);
PL_suburban = zeros(N,1);
PL_open = zeros(N,1);

% Loop through all points and calculate path loss
for i = 1:N
    % Calculate the distance from the source to the point
    d = norm(R2(i,:) - source) / 1000; % Convert from m to km
    fc=960;% frequency of transmission in MHz
    fc2=150;
    hb=50;% effective height of transmitting base station antenna in meters
    hm=2;% effective receiving mobile device antenna height in meters
    
    hm1=(1.1*log10(fc)-0.7)*hm-(1.56*log10(fc)-0.8);%Open
    % hm2=8.29*(log10(1.54*hm))^2-1.1;%Metropolitan fc<200
    hm3=3.2*(log10(11.75*hm))^2-4.92;%Metropolitan fc>200
    
    C=-2*(log10(fc/28))^2-5.4;
    C2=-4.78*(log10(fc))^2+18.33*log10(fc)-40.98;
    C3=0;
    
%     d = (1:50:100);
    
    A=69.55+26.16*log10(fc)-13.82*log10(hb)-hm1;%Open 
    % A2=69.55+26.16*log10(fc2)-13.82*log10(hb)-hm2;%Metropolitan fc<=200
    A3=69.55+26.16*log10(fc)-13.82*log10(hb)-hm3;%Metropolitan fc>=200
    
    B=44.9-6.55*log10(hb);
    
    PL_suburban(i)=A+B*log10(d)+C; %Suburban
    PL_open(i)=A+B*log10(d)+C2; %Open
    %Plmp1=A2+B*log10(d)+C3 %Metropolitan fc<=200
    PL_urban(i)=A3+B*log10(d)+C3; %Metropolitan fc>=200

    disp(i)
end

% Calculate the shift value for each area
shift_value_urban = e2(1,2) - PL_urban(1);
shift_value_suburban = e2(1,2) - PL_suburban(1);
shift_value_open = e2(1,2) - PL_open(1);

% Shift the path loss values for each area
PL_urban_shifted = PL_urban + shift_value_urban;
PL_suburban_shifted = PL_suburban + shift_value_suburban;
PL_open_shifted = PL_open + shift_value_open;

% Invert the shifted path loss values
PL_urban_shifted_inverted = max(PL_urban_shifted) - PL_urban_shifted;
PL_suburban_shifted_inverted = max(PL_suburban_shifted) - PL_suburban_shifted;
PL_open_shifted_inverted = max(PL_open_shifted) - PL_open_shifted;

% Adjust the starting point of the path loss curves
PL_urban_shifted_inverted = PL_urban_shifted_inverted - PL_urban_shifted_inverted(1) + e2(1,2) + shift_value_urban;
PL_suburban_shifted_inverted = PL_suburban_shifted_inverted - PL_suburban_shifted_inverted(1) + e2(1,2) + shift_value_suburban;
PL_open_shifted_inverted = PL_open_shifted_inverted - PL_open_shifted_inverted(1) + e2(1,2);

% % Define a linear space between the first point of the open line and the last point of the respective lines
% % PL_urban_shifted_inverted_linspace = linspace(PL_open_shifted_inverted(1), PL_urban_shifted_inverted(end), N);
% % PL_suburban_shifted_inverted_linspace = linspace(PL_open_shifted_inverted(1), PL_suburban_shifted_inverted(end), N);
% 
% % Define a linear space from 1 to N
% x = 1:N;
% 
% % Find the maximum starting point value between current and open
% max_start_point = max(e2(1,2), PL_open_shifted_inverted(1));
% 
% % Adjust the starting point of the urban and suburban curves to the maximum starting point
% PL_urban_shifted_inverted(1) = max_start_point;
% PL_suburban_shifted_inverted(1) = max_start_point;
% 
% % Now generate a new linear space for urban and suburban that starts from max_start_point
% PL_urban_shifted_inverted_linspace = linspace(max_start_point, PL_urban_shifted_inverted(end), N);
% PL_suburban_shifted_inverted_linspace = linspace(max_start_point, PL_suburban_shifted_inverted(end), N);
% 
% % Then generate a new spline interpolation
% PL_urban_shifted_inverted_spline = spline(x, PL_urban_shifted_inverted_linspace, x);
% PL_suburban_shifted_inverted_spline = spline(x, PL_suburban_shifted_inverted_linspace, x);

% Indices where PL is not NaN
validIdx_urban = ~isnan(PL_urban_shifted_inverted);
validIdx_suburban = ~isnan(PL_suburban_shifted_inverted);
validIdx_open = ~isnan(PL_open_shifted_inverted);

% Plot the path loss
plot(R2(validIdx_urban,1), PL_urban_shifted_inverted(validIdx_urban), 'DisplayName', 'urban'); hold on;
plot(R2(validIdx_suburban,1), PL_suburban_shifted_inverted(validIdx_suburban), 'DisplayName', 'suburban'); hold on;
plot(R2(validIdx_open,1), PL_open_shifted_inverted(validIdx_open), 'DisplayName', 'open'); hold on;

title('Path Loss Comparison (MATLAB)');
xlabel('Distance (m)');
ylabel('Electric Field (dB) / Path Loss (dB)');
legend('Location', 'southeast');
hold off;

% Define chunk size
chunk_size = 1000; % This is just an example. Choose a chunk size suitable for your system memory.

% Preallocate the Error arrays
Error_urban = zeros(size(e2,1),1);
Error_suburban = zeros(size(e2,1),1);
Error_open = zeros(size(e2,1),1);

% Check if PL_ variables are column vectors and if not, transpose them
if size(PL_urban_shifted_inverted_spline, 2) > 1
    PL_urban_shifted_inverted_spline = PL_urban_shifted_inverted_spline';
end
if size(PL_suburban_shifted_inverted_spline, 2) > 1
    PL_suburban_shifted_inverted_spline = PL_suburban_shifted_inverted_spline';
end

% Calculate error for each area in chunks
for idx = 1:chunk_size:size(e2,1)
    endIdx = min(idx+chunk_size-1, size(e2,1));
    
    current_indices = idx:endIdx;
    
    % Ensure vectors being subtracted are the same size
    if length(PL_urban_shifted_inverted_spline(current_indices)') ~= length(e2(current_indices,2))
        error('Vectors are not the same size at indices %d to %d.', idx, endIdx);
    end

    Error_urban(current_indices) = e2(current_indices,2) - PL_urban_shifted_inverted_spline(current_indices);
    Error_suburban(current_indices) = e2(current_indices,2) - PL_suburban_shifted_inverted_spline(current_indices);
    Error_open(current_indices) = e2(current_indices,2) - PL_open_shifted_inverted(current_indices);
end


% Square the errors
SquaredError_urban = Error_urban.^2;
SquaredError_suburban = Error_suburban.^2;
SquaredError_open = Error_open.^2;

% Compute mean squared error for each area
MSE_urban = mean(SquaredError_urban(validIdx_urban));
MSE_suburban = mean(SquaredError_suburban(validIdx_suburban));
MSE_open = mean(SquaredError_open(validIdx_open));

% Compute root mean square error for each area
RMSE_urban = sqrt(MSE_urban);
RMSE_suburban = sqrt(MSE_suburban);
RMSE_open = sqrt(MSE_open);

% Display the RMSE for each area
disp(['RMSE for Urban: ', num2str(RMSE_urban)])
disp(['RMSE for Suburban: ', num2str(RMSE_suburban)])
disp(['RMSE for Open: ', num2str(RMSE_open)])
