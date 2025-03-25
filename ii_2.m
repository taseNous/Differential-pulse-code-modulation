function yq = my_quantizer(y, N, min_value, max_value)

    % Διασφάλιση ότι min_value < max_value
    assert(min_value < max_value, 'min_value must be less than max_value');
    
    % Διασφάλιση ότι N > 0
    assert(N > 0 && floor(N) == N, 'N must be a positive integer');
    
    % Βήμα 1: Περιορισμός δυναμικής περιοχής του σήματος στις τιμές [min_value,
    % max_value]
    y_confined = max(min(y, max_value), min_value);
    
    % Βήμα 2: Υπολογισμός του βήματος κβαντισμού Δ και των κεντρών κάθε
    % περιοχής
    L = 2^N;                 
    delta = (max_value - min_value) / L;
    centers = min_value + (0:L-1) * delta;
    
    % Βήμα 3: Εύρεση περιοχής που ανήκει το κάθε δείγμα
    epsilon = 1e-12;
    yq_index = round((y_confined - min_value + epsilon) / delta);
    yq_index = max(min(yq_index, L-1), 0); % Clamp to valid index range
    
    % Βήμα 4: Χαρτογράφηση του δείκτη στην αντίστοιχη κβαντισμένη τιμή
    yq = centers(yq_index + 1);

end

close all;

% Παράμετροι
p_values = [5, 10]; % Δύο τιμέσ για το p
N_values = [1, 2, 3];
min_value = -3.5;
max_value = 3.5;

% Φόρτωση δειγμάτων
data = load('source.mat');
variables = fieldnames(data);
x = data.(variables{1});

x_length = length(x);

% Λούπα για όλα τα p και Ν
figure;
for p_idx = 1:length(p_values)
    p = p_values(p_idx);
    
    R = zeros(p, p);
    r = zeros(p, 1);
    for i = 1:p
        sum_r = 0;
        for n = p+1:x_length
            sum_r = sum_r + x(n) * x(n-i);
        end
        r(i) = sum_r / (x_length - p);

        for j = 1:p
            sum_R = 0;
            for n = p+1:x_length
                sum_R = sum_R + x(n-j) * x(n-i);
            end
            R(i, j) = sum_R / (x_length - p);
        end
    end

    a = R \ r;
    for i = 1:p
        a(i) = my_quantizer(a(i), 8, -2, 2);
    end
    
    for N_idx = 1:length(N_values)
        N = N_values(N_idx);
        
        % Αρχικοποίηση παραμέτρων για DPCM
        memory = zeros(p, 1);
        y = zeros(x_length, 1);
        y_hat = zeros(x_length, 1);
        y_toned = 0;

        for i = 1:x_length
            y(i) = x(i) - y_toned;
            y_hat(i) = my_quantizer(y(i), N, min_value, max_value);
            y_hat_toned = y_hat(i) + y_toned; 
            memory = [y_hat_toned; memory(1:p-1)]; 
            y_toned = a' * memory;
        end

        % Plot αρχικού σήματος και πρόβλεψης σφάλματος
        subplot(length(p_values), length(N_values), (p_idx-1)*length(N_values) + N_idx);
        plot(1:x_length, x, 'b', 1:x_length, y, 'r');
        title(sprintf('p = %d, N = %d bits', p, N));
        legend('Original Signal', 'Prediction Error');
        xlabel('Sample Index');
        ylabel('Amplitude');
    end
end