% function [x_Rx, Y_combined] = Rx(M, L, y_Tx1, y_Tx2, x_QAM_modulated, h1, h2, fft_size, No)
%     % Check which channel has bigger energy (smaller effect on signal)
%     % To equalize (remove channel effect) its received signal 
%     % And remove cyclic prefix
%     y1 = y_Tx1(L:fft_size+L-1);
%     y1 = y1 ./ sqrt(fft_size);
% 
%     y2 = y_Tx2(L:fft_size+L-1);
%     y2 = y2 ./ sqrt(fft_size);
% 
%     % Add noise to the channel impulse responses h1 and h2 in the time domain
%     noise_h1 = sqrt(No/2) .* (randn(size(h1)) + 1i * randn(size(h1)));
%     noise_h2 = sqrt(No/2) .* (randn(size(h2)) + 1i * randn(size(h2)));
% 
%     h1_noisy = h1 + noise_h1;
%     h2_noisy = h2 + noise_h2;
% 
%     % Compute the FFT of the channel impulse responses with noise
%     H1 = fft(h1_noisy, fft_size);
%     H2 = fft(h2_noisy, fft_size);
% 
%     % FFT of the received signals
%     Y1 = fft(y1, fft_size); % Freq. Domain
%     Y2 = fft(y2, fft_size); % Freq. Domain
% 
%     % MMSE Equalization
%     H_combined = H1 + H2; % Effective combined channel response
%     MMSE_filter = conj(H_combined) ./ (abs(H_combined).^2 + No); % MMSE filter
% 
%     % Apply MMSE filter to combine signals
%     Y_combined = (Y1 + Y2) .* MMSE_filter;
% 
%     % Demodulation
%     modulation_length = length(x_QAM_modulated);
%     x_demod = QAM_demodulation(Y_combined(1:modulation_length), M);
% 
%     % Channel Decoding
%     x_Rx = Channeldecoding(x_demod);
% end                                                         


function [x_Rx, Y_combined] = Rx(M, L, y_Tx1, y_Tx2, x_QAM_modulated, h1, h2, fft_size, No)
    % Check which channel has bigger energy (smaller effect on signal)
    % To equalize (remove channel effect) its received signal 
    % And remove cyclic prefix
    y1 = y_Tx1(L:fft_size+L-1);
    y1 = y1 ./ sqrt(fft_size);

    y2 = y_Tx2(L:fft_size+L-1);
    y2 = y2 ./ sqrt(fft_size);

    % Add noise to the channel impulse responses h1 and h2 in the time domain
    noise_h1 = 0; %sqrt(No/2) .* (randn(size(h1)) + 1i * randn(size(h1)));
    noise_h2 =0; % sqrt(No/2) .* (randn(size(h2)) + 1i * randn(size(h2)));

    h1_noisy = h1 + noise_h1;
    h2_noisy = h2 + noise_h2;

    % Compute the FFT of the channel impulse responses with noise
    H1 = fft(h1_noisy, fft_size);
    H2 = fft(h2_noisy, fft_size);

    % FFT of the received signals
    Y1 = fft(y1, fft_size); % Freq. Domain
    Y2 = fft(y2, fft_size); % Freq. Domain

    % Maximum Ratio Combining (MRC)
    Y_combined = (conj(H1) .* Y1 + conj(H2) .* Y2) ./ (abs(H1).^2 + abs(H2).^2);

    % Demodulation
    modulation_length = length(x_QAM_modulated);
    x_demod = QAM_demodulation(Y_combined(1:modulation_length), M);

    % Channel Decoding
    x_Rx = Channeldecoding(x_demod);
end