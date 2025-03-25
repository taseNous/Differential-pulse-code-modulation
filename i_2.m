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

% Αρχικοποίηση Παραμέτρων
p = 8;
N = 3;
min_value = -3.5;
max_value = 3.5;

% Φόρτωση δειγμάτων
data = load('source.mat'); 
variables = fieldnames(data); 
x = data.(variables{1}); 

x_length = length(x);

% Βήμα 1: Υπολογισμός πίνακα αυτοσυσχέτισης R και διάνυσμα αυτοσυσχέτισης r
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

% Βήμα 2: Υπολογισμός συντελεστπών και κβαντισμός τους
a = R \ r; % Yule-Walker
for i = 1:p
    a(i) = my_quantizer(a(i), 8, -2, 2);
end

% Βήμα 3: DPCM Κωδικοποίηση
memory = zeros(p, 1);
y = zeros(x_length, 1);
y_hat = zeros(x_length, 1);
y_toned = 0;

for i = 1:x_length
    y(i) = x(i) - y_toned; % Πρόβλεψη σφάλματος
    y_hat(i) = my_quantizer(y(i), N, min_value, max_value); % Κβαντισμός error
    y_hat_toned = y_hat(i) + y_toned; % Προσθήκη σφάλμα κβαντισμού στη πρόβλεψη
    memory = [y_hat_toned; memory(1:p-1)];
    y_toned = a' * memory; 
end

% Βήμα 4: DPCM Αποκωδικοποίηση
memory = zeros(p, 1);
y_toned = 0;
y_final = zeros(x_length, 1);

for i = 1:x_length
    y_final(i) = y_hat(i) + y_toned; % Ανακατασκευή
    memory = [y_final(i); memory(1:p-1)];
    y_toned = a' * memory;
end

% Plot αρχικού και ανακατασκευασμένου σήματος
figure;
plot(1:x_length, x, 'b', 1:x_length, y_final, 'r');
legend('Original Signal', 'Reconstructed Signal');
title('Original vs Reconstructed Signal');
xlabel('Sample Index');
ylabel('Amplitude');

% Plot πρόβλεψη σφάλματος
figure;
plot(1:x_length, y, 'k');
title('Prediction Error');
xlabel('Sample Index');
ylabel('Amplitude');

plot(x) % Αρχικό σήμα
hold
plot(y) % Σήμα πρόβλεψης σφάλματος
xlabel('Input signal x') 
ylabel('Error Prediction Signal y') 