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

% Define parameters
p_values = 5:10; % Predictor orders (from 5 to 10)
N_values = [1, 2, 3]; % Bits for quantization
min_value = -3.5;
max_value = 3.5;

% Load input data
data = load('source.mat');
variables = fieldnames(data);
x = data.(variables{1});

x_length = length(x);

% Initialize results
mse_results = zeros(length(p_values), length(N_values));
coefficients = cell(length(p_values), 1);

% Loop through p and N values
for p_idx = 1:length(p_values)
    p = p_values(p_idx);
    
    % Calculate R and r for this p
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

    % Solve for a and quantize coefficients
    a = R \ r; % Solve Yule-Walker equations
    
    % Apply threshold to remove insignificant coefficients
    threshold = 0.1 * max(abs(a)); % 10% of the maximum coefficient
    a(abs(a) < threshold) = 0;
    
    % Quantize coefficients
    for i = 1:p
        a(i) = my_quantizer(a(i), 10, -2, 2); % Quantize to 10 bits
    end
    coefficients{p_idx} = a; % Store quantized coefficients
    
    % Evaluate for each N
    for N_idx = 1:length(N_values)
        N = N_values(N_idx);

        % Initialize variables for DPCM
        memory = zeros(p, 1);
        y = zeros(x_length, 1);
        y_hat = zeros(x_length, 1);
        y_toned = 0;

        % DPCM Encoding
        for i = 1:x_length
            y(i) = x(i) - y_toned; % Prediction error
            y_hat(i) = my_quantizer(y(i), N, min_value, max_value); % Quantize error
            y_hat_toned = y_hat(i) + y_toned; % Add quantized error to prediction
            memory = [y_hat_toned; memory(1:p-1)]; % Update memory
            y_toned = a' * memory; % Update prediction
        end

        % Calculate MSE
        mse_results(p_idx, N_idx) = mean((x - y_hat).^2);
    end
end

% Plot MSE vs N for each p
figure;
for p_idx = 1:length(p_values)
    plot(N_values, mse_results(p_idx, :), '-o', 'DisplayName', sprintf('p = %d', p_values(p_idx)));
    hold on;
end
grid on;
title('MSE vs N for Different Predictor Orders');
xlabel('Number of Bits (N)');
ylabel('Mean Squared Error (MSE)');
legend;

% Display quantized predictor coefficients
disp('Quantized Predictor Coefficients:');
for p_idx = 1:length(p_values)
    fprintf('p = %d: ', p_values(p_idx));
    disp(coefficients{p_idx}');
end