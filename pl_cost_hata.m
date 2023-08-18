f = 1800;        % frequency in MHz
hb = 50;        % transmitter height in meters
hm = 2;       % receiver height in meters
d = 1:1:20;          % distance in km
type = "urban"; % urban and suburban

PL_urban = Cost231Hata(f, hb, hm, d, "urban");
PL_suburban = Cost231Hata(f, hb, hm, d, "suburban");

% fprintf('The path loss is %.2f dB\n', PL);

figure;
plot(d, PL_urban, 'b', 'LineWidth', 2); hold on;
plot(d, PL_suburban, 'r', 'LineWidth', 2); hold on;

title('Path Loss for Different Area Types');
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
legend('Urban', 'Suburban', 'Location', 'northwest');
grid on;


% Cost231Hata: calculates the path loss using the Cost231-Hata model
% 
% L = F + BlogR - E + G
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
    if f < 1500 || f > 2000
        error('Frequency must be between 1500 and 2000 MHz');
    end
    if hb < 30 || hb > 200
        error('Transmitter height must be between 30 and 200 meters');
    end
    if hm < 1 || hm > 10
        error('Receiver height must be between 1 and 10 meters');
    end
    if d < 1
        error('Distance must be minimum 1 km');
    end
    if "urban" == type
        G = 3;
    else
        G = 0;
    end

    B = 44.9 - 6.55*log10(hb);
    E = (1.1*log10(f) - 0.7)*hm - (1.56*log10(f) - 0.8);
    F = 46.3 + 33.9*log10(f) - 13.82*log10(hb);
    
    % Calculate path loss
    PL = F + B*log10(d) - E + G;
end
