
% Sinal
%sin
N0=120; % número de amostras
fs=12000; % frequência de amostragem
t0=(0:N0-1)/fs; % vetor do tempo
ns= 3; % número de aquisicoes

A=1.5; % amplitude
fase=45; % fase em graus
f=462; % frequência
%noise
n_rms=0.01; % valor eficaz ruído
ruido = n_rms*randn(1,N0);

ua=A*cos(2*pi*f*t0+fase*pi/180);

u0(1,:)=A*cos(2*pi*f*t0+fase*pi/180)+ ruido; %sinal
u0(2,:)=A*cos(2*pi*f*t0+(fase-15)*pi/180)+ ruido;
u0(3,:)=A*cos(2*pi*f*t0+(fase-20)*pi/180) + ruido;

% Analise do sinal
% Necessário:   N0   (Número de amostras)
%               fs  (Frequência de amostragem)
%               u0   (Sinal)

% Para obter melhores resultados
N_discarded = 0;
iterations = 3;
unilateral_total = 0;
for inout= 1:ns
u = u0(inout,:);
t = t0;
while iterations > 0

N = N0 - N_discarded;
t = t(1:end -N_discarded);
u = u(1:end -N_discarded);

% Espectro
 % falar sobre as amostras terem que ser pares ou ímpares
 delta_f = fs/N;
 freqs = fs*(0:(N/2))/N;
 
 bilateral = fft(u)/N; %FFT de u abs(fft(u)/N)
 unilateral = bilateral(1:floor(N/2)+1);
 unilateral(2:end-1) = 2*unilateral(2:end-1); % duplicar todos os termos excepto o primeiro e último
 unilateral_ef = abs(unilateral)/sqrt(2);
 unilateral_db = 20*log10(unilateral_ef); % cada sinusoide tem de ser multiplicado por 1/sqrt(2) para o valor eficaz

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
 % Excluir primeira e última posição
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
 
 
 % peakB_id é da frequência mais alta
 if anterior < posterior
    peakB_id = peakA_id+1;
 else
    peakB_id = peakA_id;
    peakA_id = peakA_id-1;
 end
 
 if iterations == 3 
 % estimativa do valor da frequencia, só é feito uma vez
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
 iterations = iterations - 1;
 end
 iterations = iterations - 1;
end
 iterations = 1;
 unilateral_total = unilateral_total+unilateral_ef; %soma dos espectros
end
 % realizar a média dos espectros
 unilateral_total = unilateral_total/ns;
 unilateral_total_db = 20*log10(unilateral_total);
 
 % RMS
 valor_rms = 0;
 for i = 1:ns
 valor_rms = sqrt(sum(u0(i,1:end-N_discarded).^2)/N)+valor_rms;
 end
 valor_rms = valor_rms/ns;
 
 % noise and distortion RMS
 % este valor é muito dependente do espalhamento, depende da freq 35-40 db
 nd_rms = sqrt(valor_rms^2-(unilateral_total(peak_id))^2);
 
 % SINAD
 sinad = valor_rms/nd_rms;
 sinad_db = 20*log10(sinad);
 
 % ENOB
 enob = (sinad-1.76)/6.02;
 
 % Plots
 fig = figure(1);
 s1 = sprintf('Frequência estimada do sinal: %f [Hz]      ', freq_estimada);
 s2 = sprintf('Valor eficaz do sinal: %f [V] \n', valor_rms);
 s3 = sprintf('Valor eficaz do ruído: %f [V]      ', nd_rms);
 s4 = sprintf('SINAD: %f [dB] \n', sinad_db);
 s5 = sprintf('ENOB: %f bits      ', enob);
 string = sprintf('%s%s%s%s%s',s1,s2,s3,s4,s5);
 subplot(2,1,1), plot(t(1:3*floor(pontos_p_periodo)),u(1:3*floor(pontos_p_periodo)));
 title(string);     % 3 periodos do sinal
 xlim([0 t(3*floor(pontos_p_periodo))]);
 xlabel('Tempo [s]');
 ylabel('Amplitude [V]');
 subplot(2,1,2), plot(freqs,unilateral_total_db);
 xlabel('Frequência [Hz]');
 ylabel('Ampltitude [dB]');