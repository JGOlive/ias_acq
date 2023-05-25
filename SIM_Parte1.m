
% Sinal
%sin
N0=120; % número de amostras
fs=12000; % frequência de amostragem
t0=(0:N0-1)/fs; % vetor do tempo

A=1.5; % amplitude
fase=45; % fase em graus
f=467; % frequência
%noise
n_rms=0.01; % valor eficaz ruído
ruido = n_rms*randn(1,N0);

u0=A*cos(2*pi*f*t0+fase*pi/180)+ ruido; %sinal

% Analise do sinal
% Necessário:   N0   (Número de amostras)
%               fs  (Frequência de amostragem)
%               u0   (Sinal)

% Para obter melhores resultados
delta_t = (t0(end)-t0(1))/(N0-1);
N_discarded = 0;
iterations = 2;
while iterations > 0
iterations = iterations - 1;
N = N0 - N_discarded;          % caso necessário serão descartadas amostras
t = t0(1:end -N_discarded);    %na segunda iteração
u = u0(1:end -N_discarded);

% Espectro
 % falar sobre as amostras terem que ser pares ou ímpares
 delta_f = fs/N;
 freqs = fs*(0:N/2)/N;
 
 bilateral = fft(u)/N; %FFT de u abs(fft(u)/N)
 unilateral = bilateral(1:floor(N/2)+1); % o floor é para quando N é impar
 unilateral(2:end-1) = 2*unilateral(2:end-1); % duplicar todos os termos 
 unilateral_ef = abs(unilateral)/sqrt(2);     %excepto o primeiro e último
 unilateral_db = 20*log10(unilateral_ef); % cada sinusoide tem de ser 
                                          %multiplicado por 1/sqrt(2)
% Frequência

 % encontrar o maximo
 peak = -inf;
 peak_id = 0;
 

 for i = 1:(N/2-1)
    if peak < unilateral_db(i)
        peak_id = i;
        peak = unilateral_db(i);
    end
 end
 
 peakA_id = peak_id;
 % Excluir primeira e última posição, indices vizinhos ultapassam o vetor 
 if peakA_id == 1
    anterior = -inf;
    posterior = unilateral_db(peakA_id+1);
 elseif peakA_id == (N/2)+1
    posterior = -inf;
    anterior = unilateral_db(peakA_id-1);
 else
    anterior = unilateral_db(peakA_id-1);
    posterior = unilateral_db(peakA_id+1);  
 end
 
 
 % peakB_id é da frequência mais alta dos dois mais altos picos
 if anterior < posterior
    peakB_id = peakA_id+1;
 else
    peakB_id = peakA_id;
    peakA_id = peakA_id-1;
 end
 
 % obter a frequência a partir do sinc
 omega_A = 2*pi*(peakA_id-1)*delta_f/fs;
 omega_B = 2*pi*(peakB_id-1)*delta_f/fs;
 XA = unilateral(peakA_id);
 UA = real(XA);             % parte real
 VA = imag(XA);             % parte imaginaria
 XB = unilateral(peakB_id);
 UB = real(XB);             % parte real
 VB = imag(XB);             % parte imaginaria
 
 Kopt = ((VB-VA)*sin(omega_A)+(UB-UA)*cos(omega_A))/(UB-UA);
 ZA = VA*(Kopt-cos(omega_A))/sin(omega_A)+UA;
 ZB = VB*(Kopt-cos(omega_B))/sin(omega_B)+UB;
 freq_estimada = (fs/(2*pi))*acos((ZB*cos(omega_B)-ZA*cos(omega_A))/(ZB-ZA));
 
 % info para fazer nova FFT
 
 pontos_p_periodo = fs/freq_estimada;
 N_discarded = round(rem(N,pontos_p_periodo)); %pontos por discartar
 
 end
 
 % Mean value
 valor_medio = sum(u)/N;   % u já com amostras descartadas
 
 % RMS
 valor_rms = sqrt(sum(u.^2)/N);
 
 % THD - tem que se fazer (rms das harmónicas)/(rms da fundamental)
 valor_thd = sqrt(sum(unilateral_ef(2*peak_id-1:peak_id-1:end).^2))/(unilateral_ef(peak_id));
 valor_thd_db = 20*log10(valor_thd);
 
 % Plots
 fig = figure(1);
 s1 = sprintf('Resolução temporal: %f [s]      ', delta_t);
 s2 = sprintf('Resolução espectral: %f [Hz]\n', delta_f);
 s3 = sprintf('Frequência estimada do sinal: %f [Hz]      ', freq_estimada);
 s4 = sprintf('Valor eficaz do sinal: %f [V] \n', valor_rms);
 s5 = sprintf('Valor médio do sinal: %f [V]      ', valor_medio);
 s6 = sprintf('THD: %f [dB]', valor_thd_db);
 string = sprintf('%s%s%s%s%s',s1,s2,s3,s4,s5,s6);
 subplot(2,1,1), plot(t(1:3*floor(pontos_p_periodo)),u(1:3*floor(pontos_p_periodo)));
 title(string);     % 3 períodos do sinal na figura
 xlim([0 t(3*floor(pontos_p_periodo))]);
 xlabel('Tempo [s]');
 ylabel('Amplitude [V]');
 subplot(2,1,2), plot(freqs,unilateral_db);
 xlabel('Frequência [Hz]');
 ylabel('Ampltitude [dB]');