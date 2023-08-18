f = 960;        % frequency in MHz
hb = 50;        % transmitter height in meters
hm = 2;       % receiver height in meters
d = 1:1:20;          % distance in km

PL_urban = OkumuraHata(f, hb, hm, d, "urban");
PL_suburban = OkumuraHata(f, hb, hm, d, "suburban");
PL_open = OkumuraHata(f, hb, hm, d, "open");

% fprintf('The path loss for %s type is %.2f dB\n', type, PL);

figure;
plot(d, PL_urban, 'b', 'LineWidth', 2); hold on;
plot(d, PL_suburban, 'r', 'LineWidth', 2); hold on;
plot(d, PL_open, 'g', 'LineWidth', 2); hold off;

title('Path Loss for Different Area Types');
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
legend('Urban', 'Suburban', 'Open', 'Location', 'northwest');
grid on;

% OkumuraHata: calculates the path loss using the Okumura-Hata model
%
% INPUTS:
% f  : frequency in MHz, must be between 150 and 1500
% hb : height of the transmitter in meters, must be between 30 and 200
% hm : height of the receiver in meters, must be between 1 and 10
% d  : distance between transmitter and receiver in km, must be between 1 and 20
%
% OUTPUTS:
% PL : path loss in dB
function PL = OkumuraHata(f, hb, hm, d, type)
    % Ensure inputs are within the correct range
    if f < 150 || f > 1500
        error('Frequency must be between 150 and 1500 MHz');
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

    % calculate correction factor based on area type
    if "urban" == type
        % Urban Areas (medium to small cities)
        if f < 300
            cf = (1.1*log10(f) - 0.7)*hm - (1.56*log10(f) - 0.8);
        else
            cf = 3.2*((log10(11.75*hm))^2) - 4.97;
        end
    elseif "suburban" == type
        % Suburban Areas
        cf = 2*((log10(f/28))^2) - 5.4;
    else
        % Open Areas
        cf = 4.78*(log10(f))^2 + 18.33*log10(f) + 40.94;
    end

    A = 69.55 + 26.16*log10(f) - 13.82*log10(hb);
    B = 44.9 - 6.55*log10(hb);
    
    % Calculate path loss
    PL = A + (B)*log10(d) - cf;
end