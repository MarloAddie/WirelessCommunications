clc, clear, close all


L = 5;
M = 16;
N = 4;
cp = L-1;

numBits = N*M*2; % # of QPSK info symbols

% Generate rand binary data
randSeq = rand(1, numBits);
binaryData = (randSeq > 0.5);



real = binaryData(1:2:end);
imag = binaryData(2:2:end);


% Create indexing for channel matrix based on knowlegde of signal
x = zeros(1, M);
for lk = 1:M
    x(lk) = lk;
end

new_x = zeros(1, length(x)+cp);
new_x(1:cp) = x((end-cp+1):end);
new_x(cp+1:end) = x(1:end);

% QPSK mapping
d_tilde = zeros(length(real),1);
for i = 1:length(real)
    if (real(i) == imag(i))
        re = 1;
    else
        re = -1;
    end
    if (real(i) < 1)
        im = 1i;
    else
        im = -1i;
    end
    d_tilde(i) = 1/sqrt(2)*(re+im);
end

h_k = [0.227,0.46,0.688,0.46,0.227];
% channel modeling designed for a signal with a cyclic prefix
h = zeros(M,M); % 16x16 matrix for the modeled channel
for i = cp+2:M+cp+1 % start after the cyclic prefix (adds and removes cp in one step)
    for k = 1:L % go through the taps
        jk = new_x(i - k);  
            h(i-5,jk) = h_k(k);
    end
end

 f_i_j = zeros(M,M);
% generate F_N matrix:
for i = 1:M
    for J = 1:M
        f_i_j(i,J) = exp((-1j*2*pi*(i-1)*(J-1))/M);
    end
end

F = (1/sqrt(M)) .* f_i_j; % F_N matrix

F_I = F'*F;





y_tilda = zeros(64,16);
for sub = 1:4
jklm = sub+(sub-1)*(16-1); % weird indexing shoudl fix
dd = d_tilde(jklm:jklm+15,1);
d = F'*dd; % transmitted signal
y = h*d; % received signal (time domain)
y_tilda(jklm:jklm+15,1) = F*y; % received signal (freq domain)

end

y_tilde = y_tilda(1:16,1);

% verify correctness using one-tap equalizer:
h_tilde = (F*h*F'); % freq domain channel
h_tilde_diag = diag(h_tilde,0); % one tap equalizer (diagonal of h_tilde)
d_tilde_rx = y_tilde./h_tilde_diag; % received symbols (equal to d_tilde)



