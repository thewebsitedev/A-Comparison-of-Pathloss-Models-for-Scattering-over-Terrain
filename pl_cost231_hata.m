figure;
plot(e2(:,1),e2(:,2), 'DisplayName', 'current'); hold on; % Add a display name to the plot

% Define your frequency and antenna heights
f = 960;        % frequency in MHz
hb = 50;        % transmitter height in meters
hm = 2;       % receiver height in meters

% Define your areas
areas = ["urban", "suburban", "open"];

PL_urban = zeros(N,1);
PL_suburban = zeros(N,1);
PL_open = zeros(N,1);

% Loop through all points and calculate path loss
for area = areas
    disp(area)
    for i = 1:N
        % Calculate the distance from the source to the point
        d = norm(R2(i,:) - source) / 1000; % Convert from m to km
        if "urban" == area
            PL_urban(i) = Cost231Hata(f, hb, hm, d, area);
        elseif "suburban" == area
            PL_suburban(i) = Cost231Hata(f, hb, hm, d, area);
        else
            PL_open(i) = Cost231Hata(f, hb, hm, d, area);
        end
    end
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

% Define a linear space between the first point of the open line and the last point of the respective lines
% PL_urban_shifted_inverted_linspace = linspace(PL_open_shifted_inverted(1), PL_urban_shifted_inverted(end), N);
% PL_suburban_shifted_inverted_linspace = linspace(PL_open_shifted_inverted(1), PL_suburban_shifted_inverted(end), N);

% Define a linear space from 1 to N
x = 1:N;

% Find the maximum starting point value between current and open
max_start_point = max(e2(1,2), PL_open_shifted_inverted(1));

% Adjust the starting point of the urban and suburban curves to the maximum starting point
PL_urban_shifted_inverted(1) = max_start_point;
PL_suburban_shifted_inverted(1) = max_start_point;

% Now generate a new linear space for urban and suburban that starts from max_start_point
PL_urban_shifted_inverted_linspace = linspace(max_start_point, PL_urban_shifted_inverted(end), N);
PL_suburban_shifted_inverted_linspace = linspace(max_start_point, PL_suburban_shifted_inverted(end), N);

% Then generate a new spline interpolation
PL_urban_shifted_inverted_spline = spline(x, PL_urban_shifted_inverted_linspace, x);
PL_suburban_shifted_inverted_spline = spline(x, PL_suburban_shifted_inverted_linspace, x);

% Indices where PL is not NaN
validIdx_urban = ~isnan(PL_urban_shifted_inverted);
validIdx_suburban = ~isnan(PL_suburban_shifted_inverted);
validIdx_open = ~isnan(PL_open_shifted_inverted);

% Plot the path loss
plot(R2(validIdx_urban,1), PL_urban_shifted_inverted(validIdx_urban), 'DisplayName', 'urban'); hold on;
plot(R2(validIdx_suburban,1), PL_suburban_shifted_inverted(validIdx_suburban), 'DisplayName', 'suburban'); hold on;
plot(R2(validIdx_open,1), PL_open_shifted_inverted(validIdx_open), 'DisplayName', 'open'); hold on;

title('Path Loss Comparison');
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

% Cost231Hata: calculates the path loss using the Okumura-Hata model
%
% INPUTS:
% f  : frequency in MHz, must be between 150 and 1500
% hb : height of the transmitter in meters, must be between 30 and 200
% hm : height of the receiver in meters, must be between 1 and 10
% d  : distance between transmitter and receiver in km, must be between 1 and 20
%
% OUTPUTS:
% PL : path loss in dB
function PL = Cost231Hata(f, hb, hm, d, type)
    % Ensure inputs are within the correct range
    if f < 150 || f > 2000
        error('Frequency must be between 150 and 1500 MHz');
    end
    if hb < 30 || hb > 200
        error('Transmitter height must be between 30 and 200 meters');
    end
    if hm < 1 || hm > 10
        error('Receiver height must be between 1 and 10 meters');
    end
%     if d < 1
%         error('Distance must be minimum 1 km');
%     end

    G = 0;

    % calculate correction factor based on area type
    if "urban" == type
        % Urban Areas (medium to small cities)
        if f < 300
            cf = (1.1*log10(f) - 0.7)*hm - (1.56*log10(f) - 0.8);
        else
            cf = 3.2*((log10(11.75*hm))^2) - 4.97;
        end
        G = 3;
    elseif "suburban" == type
        % Suburban Areas
        cf = 2*((log10(f/28))^2) - 5.4;
    else
        % Open Areas
        cf = 4.78*(log10(f))^2 + 18.33*log10(f) + 40.94;
    end

    A = 46.3 + 33.9*log10(f) - 13.82*log10(hb);
    B = 44.9 - 6.55*log10(hb);
    
    % Calculate path loss
    PL = A + (B)*log10(d) - cf + G;
end