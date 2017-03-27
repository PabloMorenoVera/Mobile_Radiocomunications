% -------------------------------------------------------------------------
%                          PRÁCTICA 1. HITO 2.
%
%                Sistema con codificación de canal BCH.
% -------------------------------------------------------------------------

clear all;

% Par?metros simulación
    n = 15; % Long. palabra del código BCH.
    k = 5; % Long. palabra de entrada para el código.
    L = 1e4;
    N = k*L; % Número de bits totales.

    Px = 1;   % Potencia de la señal transmitida: normalizada
    sigm_h = 1; % Varianza del canal: normalizado.
    vector_SNR = (0:3:9); % Relaciones señal a ruido
    s1 = 1; % Símbolos BPSK
    s2 = -1;


    ber = zeros(size(vector_SNR));
    m = 1;
    for SNR = vector_SNR,
    m
        % Generación de los bits
        b = randi([0 1],L,k); % Generamos una matriz de L filas y k columnas de bits.

        % Codificación BCH
        msg = gf(b);
        c = bchenc(msg,n,k);
        bits_cod = double(c.x);

        % Modulación
        x = ((bits_cod-0.5).*2) * sqrt(Px); % Modulador muy sencillo BPSK. Obsérvense las dimensiones de x.

        % Generación de los coeficientes del canal.
        h = (randn(size(x)) + 1i.*randn(size(x))) * sqrt(sigm_h/2); % h es Rayleigh y tiene las dimensiones de x.

        % CANAL PARA EL HITO 2.3
        % Ahora asumimos que h se mantiene constante durante un bloque de
        % codificación de 7 bits, es decir, en cada fila es siempre igual,
        % por lo que generamos 1 muestra de canal para cada bloque de
        % codificación (cada fila) y luego repetimos esa muestra n veces.
        %h = ( randn(  ) + j.*randn(  ) ) * sqrt(sigm_h/2); % Vector columna
        %h = repmat(h, 1, n); % Matriz L*n

        % Generación del ruido
        snr = 10^(SNR/10);
        w =  (randn(size(x)) + j.*randn(size(x))) * (1/sqrt(2)) * sqrt(1/snr); % el ruido tiene las mismas dimensiones que x;

        % Introducimos la señal en el canal.
        y = x.*h + w; % La señal recibida tendrá dimensiones L*n

        % Recepción ML (x_ML = argmin{|y-hx|^2}
        arg1 = abs(y - h.*s1);
        arg2 = abs(y - h.*s2);
        x_ML = ((arg1<arg2)*2)-1;
        bits_recibidos = (x_ML - s2) / ( s1 - s2 );

        % Decodificación BCH (MIRAR ENUNCIADO DE LA PRÁCTICA)
        c_dec = bchdec(gf(bits_recibidos),n,k);
        br_dec = double(c_dec.x);

        % Cálculo de la BER
        ber(m) =  mean((x_ML~=x));

    m=m+1 
    end

figure;    
semilogy(vector_SNR, ber, 'r');

% PUEDE AÑADIR AQUÍ EL CÓDIGO DEL OTRO HITO PARA LA SIMULACIÓN DEL SISTEMA SOBRE 
% CANAL RAYLEIGH SIN CODIFICACIÓN.