% -------------------------------------------------------------------------
%                          PRÁCTICA 1. HITO 1.
%
%           Receptores con diversidad: Combinador Selectivo y 
%                        Maximum Ratio Combining.
% -------------------------------------------------------------------------

% Parámetros simulación
    N = 100000;  % Numero de símbolos.
    Px = 1;   % Potencia de la señal transmitida: normalizada
    sigm_h = 1; % Varianza del canl: normalizado.
    vector_SNR = (0:2:20); % Relaciones señal a ruido
    s1 = 1; % Símbolos BPSK
    s2 = -1;

%Simulación

    % Inicialización vectores
    h = zeros(1,N);
    n = zeros(1,N);
    y = zeros(1,N);
    y_filt = zeros(1,N);
    ber = zeros(size(vector_SNR));
    
    % Bucle sobre SNR.
    k = 1;
    for SNR = vector_SNR,
    
        % Generación de señal (BPSK, Potencia Unidad)
        x = ((round(rand(1,N))* 2) - 1) * sqrt(Px);

        % Generación de los coeficientes del canal Rayleigh
        h1 = (randn(1,N) + 1i.*randn(1,N)) * sqrt(1/2);
        h2 = 0;     %% Ponemos 0 para que solo coja una antena.

        % Generación del ruido
        snr = 10^(SNR/10);
        n1 = (randn(size(x)) + 1i.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr);
        n2 = (randn(size(x)) + 1i.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr);

        % Introducimos la señal en el canal (es un canal plano)
        y1 = x.*h1 + n1;
        y2 = x.*h2 + n2;
        
         % COMBINADOR SC:
         idx1 = find(abs(h1)>abs(h2)); % Instantes en los que el canal h1 es mejor que el h2.
         h_COMBINADOR(idx1) = h1(idx1); % Cuando es mejor h1: h_COMBINADOR = h1
         y_SC(idx1) = y1(idx1);     % Cuando es mejor h1: el receptor se queda con y1.

        % DETECTOR ML:
        arg1 = abs(y_SC - h_COMBINADOR.*s1);
        arg2 = abs(y_SC - h_COMBINADOR.*s2);
        x_SC = ((arg1<arg2)*2)-1; % Un poco tricky, pero vale.

        % Cálculo de la BER
        ber_SC1(k) = mean((x_SC~=x));

    k=k+1;    
    end
    
    %% Para SC
    % Parámetros simulación
    N = 100000;  % Numero de símbolos.
    Px = 1;   % Potencia de la señal transmitida: normalizada
    sigm_h = 1; % Varianza del canl: normalizado.
    vector_SNR = (0:2:20); % Relaciones señal a ruido
    s1 = 1; % Símbolos BPSK
    s2 = -1;

%Simulación

    % Inicialización vectores
    h = zeros(1,N);
    n = zeros(1,N);
    y = zeros(1,N);
    y_filt = zeros(1,N);
    ber = zeros(size(vector_SNR));
    
    % Bucle sobre SNR.
    k = 1;
    for SNR = vector_SNR,
    
        % Generación de señal (BPSK, Potencia Unidad)
        x = ((round(rand(1,N))* 2) - 1) * sqrt(Px);

        % Generación de los coeficientes del canal Rayleigh
        h1 = (randn(1,N) + 1i.*randn(1,N)) * sqrt(1/2);
        h2 = (randn(1,N) + 1i.*randn(1,N)) * sqrt(1/2);

        % Generación del ruido
        snr = 10^(SNR/10);
        n1 = (randn(size(x)) + 1i.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr);
        n2 = (randn(size(x)) + 1i.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr);

        % Introducimos la señal en el canal (es un canal plano)
        y1 = x.*h1 + n1;
        y2 = x.*h2 + n2;
        
        %% DESCOMENTAR UNO DE LOS DOS COMBINADORES:

         % COMBINADOR SC:
             idx1 = find(abs(h1)>abs(h2)); % Instantes en los que el canal h1 es mejor que el h2.
             idx2 = find(abs(h2)>abs(h1)); % Instantes en los que el canal h2 es mejor que el h1.
             
             h_COMBINADOR(idx1) = h1(idx1); % Cuando es mejor h1: h_COMBINADOR = h1
             h_COMBINADOR(idx2) = h2(idx2); % Cuando es mejor h2: h_COMBINADOR = h2
             
             y_SC(idx1) = y1(idx1);     % Cuando es mejor h1: el receptor se queda con y1.
             y_SC(idx2) = y2(idx2);     % Cuando es mejor h2: el receptor se queda con y2.

            
            
        % DETECTOR ML:
        arg1 = abs(y_SC - h_COMBINADOR.*s1);
        arg2 = abs(y_SC - h_COMBINADOR.*s2);
        x_SC2 = ((arg1<arg2)*2)-1; % Un poco tricky, pero vale.

        % Cálculo de la BER
        ber_SC2(k) = mean((x_SC2~=x));

    k=k+1;    
    end

    % Parámetros simulación
    N = 100000;  % Numero de símbolos.
    Px = 1;   % Potencia de la señal transmitida: normalizada
    sigm_h = 1; % Varianza del canl: normalizado.
    vector_SNR = (0:2:20); % Relaciones señal a ruido
    s1 = 1; % Símbolos BPSK
    s2 = -1;

%Simulación

    % Inicialización vectores
    h = zeros(1,N);
    n = zeros(1,N);
    y = zeros(1,N);
    y_filt = zeros(1,N);
    ber = zeros(size(vector_SNR));
    
    % Bucle sobre SNR.
    k = 1;
    for SNR = vector_SNR,
    
        % Generación de señal (BPSK, Potencia Unidad)
        x = ((round(rand(1,N))* 2) - 1) * sqrt(Px);

        % Generación de los coeficientes del canal Rayleigh
        h1 = (randn(1,N) + 1i.*randn(1,N)) * sqrt(1/2);
        h2 = (randn(1,N) + 1i.*randn(1,N)) * sqrt(1/2);
        h3 = (randn(1,N) + 1i.*randn(1,N)) * sqrt(1/2);

        % Generación del ruido
        snr = 10^(SNR/10);
        n1 = (randn(size(x)) + 1i.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr);
        n2 = (randn(size(x)) + 1i.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr);
        n3 = (randn(size(x)) + 1i.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr);

        % Introducimos la señal en el canal (es un canal plano)
        y1 = x.*h1 + n1;
        y2 = x.*h2 + n2;
        y3 = x.*h3 + n3;
        
        %% DESCOMENTAR UNO DE LOS DOS COMBINADORES:
         
         % COMBINADOR SC:
             idx1 = find(abs(h1)>abs(h2) & abs(h1)>abs(h3)); % Instantes en los que el canal h1 es mejor que el h2 y h3.
             idx2 = find(abs(h2)>abs(h1) &  abs(h2)>abs(h3)); % Instantes en los que el canal h2 es mejor que el h1 y h3.
             idx3 = find(abs(h3)>abs(h1) & abs(h3)>abs(h2)); % Instantes en los que el canal h2 es mejor que el h1 y h2.
             
             h_COMBINADOR(idx1) = h1(idx1); % Cuando es mejor h1: h_COMBINADOR = h1
             h_COMBINADOR(idx2) = h2(idx2); % Cuando es mejor h2: h_COMBINADOR = h2
             h_COMBINADOR(idx3) = h3(idx3); % Cuando es mejor h2: h_COMBINADOR = h2
             
             y_SC(idx1) = y1(idx1);     % Cuando es mejor h1: el receptor se queda con y1.
             y_SC(idx2) = y2(idx2);     % Cuando es mejor h2: el receptor se queda con y2.
             y_SC(idx3) = y3(idx3);     % Cuando es mejor h3: el receptor se queda con y3.
            
            
        % DETECTOR ML:
        arg1 = abs(y_SC - h_COMBINADOR.*s1);
        arg2 = abs(y_SC - h_COMBINADOR.*s2);
        x_SC3 = ((arg1<arg2)*2)-1; % Un poco tricky, pero vale.

        % Cálculo de la BER
        ber_SC3(k) = mean((x_SC3~=x));

    k=k+1;    
    end
     
    % Parámetros simulación
    N = 100000;  % Numero de símbolos.
    Px = 1;   % Potencia de la señal transmitida: normalizada
    sigm_h = 1; % Varianza del canl: normalizado.
    vector_SNR = (0:2:20); % Relaciones señal a ruido
    s1 = 1; % Símbolos BPSK
    s2 = -1;

%Simulación

    % Inicialización vectores
    h = zeros(1,N);
    n = zeros(1,N);
    y = zeros(1,N);
    y_filt = zeros(1,N);
    ber = zeros(size(vector_SNR));
    
    % Bucle sobre SNR.
    k = 1;
    for SNR = vector_SNR,
    
        % Generación de señal (BPSK, Potencia Unidad)
        x = ((round(rand(1,N))* 2) - 1) * sqrt(Px);

        % Generación de los coeficientes del canal Rayleigh
        h1 = (randn(1,N) + 1i.*randn(1,N)) * sqrt(1/2);
        h2 = (randn(1,N) + 1i.*randn(1,N)) * sqrt(1/2);

        % Generación del ruido
        snr = 10^(SNR/10);
        n1 = (randn(size(x)) + 1i.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr);
        n2 = (randn(size(x)) + 1i.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr);

        % Introducimos la señal en el canal (es un canal plano)
        y1 = x.*h1 + n1;
        y2 = x.*h2 + n2;
        
        %% DESCOMENTAR UNO DE LOS DOS COMBINADORES:
        
         % COMBINADOR MRC:
             % Filtro (multiplicación por el conjugado del canal, en cada rama)
             y1_filt = y1 .* conj(h1);
             y2_filt = y2 .* conj(h2);
 
             % Suma de las señales recibidas filtradas.
             y_MRC = y1_filt + y2_filt; 
         
             % Canal equivalente MRC
             h_COMBINADOR = abs(h1).^2 + abs(h2).^2;
            
        % DETECTOR ML:
        arg1 = abs(y_MRC - h_COMBINADOR.*s1);
        arg2 = abs(y_MRC - h_COMBINADOR.*s2);
        x_MRC = ((arg1<arg2)*2)-1; % Un poco tricky, pero vale.

        % Cálculo de la BER
        ber_MRC(k) = mean((x_MRC~=x));

    k=k+1;    
    end
    
    % Parámetros simulación
    N = 100000;  % Numero de símbolos.
    Px = 1;   % Potencia de la señal transmitida: normalizada
    sigm_h = 1; % Varianza del canl: normalizado.
    vector_SNR = (0:2:20); % Relaciones señal a ruido
    s1 = 1; % Símbolos BPSK
    s2 = -1;

%Simulación

    % Inicialización vectores
    h = zeros(1,N);
    n = zeros(1,N);
    y = zeros(1,N);
    y_filt = zeros(1,N);
    ber = zeros(size(vector_SNR));
    
    % Bucle sobre SNR.
    k = 1;
    for SNR = vector_SNR,
    
        % Generación de señal (BPSK, Potencia Unidad)
        x = ((round(rand(1,N))* 2) - 1) * sqrt(Px);

        % Generación de los coeficientes del canal Rayleigh
        h1 = (randn(1,N) + 1i.*randn(1,N)) * sqrt(1/2);
        h2 = (randn(1,N) + 1i.*randn(1,N)) * sqrt(1/2);
        h3 = (randn(1,N) + 1i.*randn(1,N)) * sqrt(1/2);

        % Generación del ruido
        snr = 10^(SNR/10);
        n1 = (randn(size(x)) + 1i.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr);
        n2 = (randn(size(x)) + 1i.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr);
        n3 = (randn(size(x)) + 1i.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr);

        % Introducimos la señal en el canal (es un canal plano)
        y1 = x.*h1 + n1;
        y2 = x.*h2 + n2;
        y3 = x.*h3 + n3;
        
        %% DESCOMENTAR UNO DE LOS DOS COMBINADORES:
        
         % COMBINADOR MRC:
             % Filtro (multiplicación por el conjugado del canal, en cada rama)
             y1_filt = y1 .* conj(h1);
             y2_filt = y2 .* conj(h2);
             y3_filt = y3 .* conj(h3);
 
             % Suma de las señales recibidas filtradas.
             y_MRC = y1_filt + y2_filt + y3_filt; 
         
             % Canal equivalente MRC
             h_COMBINADOR = h1.*conj(h1) + h2.*conj(h2) + h3.*conj(h3);
            
        % DETECTOR ML:
        arg1 = abs(y_MRC - h_COMBINADOR.*s1);
        arg2 = abs(y_MRC - h_COMBINADOR.*s2);
        x_MRC3 = ((arg1<arg2)*2)-1; % Un poco tricky, pero vale.

        % Cálculo de la BER
        ber_MRC3(k) = mean((x_MRC3~=x));

    k=k+1;    
    end
    
    semilogy(vector_SNR, ber_SC1, '-b');
    hold on;
    semilogy(vector_SNR, ber_SC2, '-r');
    hold on;
    semilogy(vector_SNR, ber_SC3, '-g');
    hold on;
    semilogy(vector_SNR, ber_MRC, '-m');
    hold on;
    semilogy(vector_SNR, ber_MRC3, '-y');
    legend('Ber SC_1', 'Ber SC_2', 'Ber SC_3', 'Ber MRC', 'Ber MRC_3');
    xlabel('SNR(dB)');
    ylabel('BER');