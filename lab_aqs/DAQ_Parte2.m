dispositivo = daqlist;
placa_id = dispositivo.DeviceID;
placa_vendor = dispositivo.VendorID;
sample_rate = 10000; % fs da daq
number_samples = 24000; % numero de amostras daq
dq = daq(dispositivo.VendorID);
dq.Rate = sample_rate;
addinput(dq, placa_id, "ai0", "Voltage");
ns = 5; %numero de aquisicoes
for a = 1:ns
[scanData, timeStamp] = read(dq,number_samples,"OutputFormat","Matrix");
ch0 = scanData(:,1);
u0(a,:) = ch0.';
t0(a,:) = timeStamp.';
end
N0 = number_samples;
fs = 1/(timeStamp(2)-timeStamp(1));
% Analise do sinal
% Necessario: N0 (Numero de amostras)
% fs (Frequencia de amostragem)
% u0 (Sinal)
% t0
% Para obter melhores resultados
N_discarded = 0;
iterations = 3;
unilateral_total = 0;
for inout= 1:ns
u = u0(inout,:);
t = t0(inout,:); %mudado
while iterations > 0
    
N = N0 - N_discarded;
t = t(1:end -N_discarded);
u = u(1:end -N_discarded);

% Espectro
% falar sobre as amostras terem que ser pares ou impares
 delta_f = fs/N;
 freqs = fs*(0:(N/2))/N;

 bilateral = fft(u)/N; %FFT de u abs(fft(u)/N)
 unilateral = bilateral(1:floor(N/2)+1);
 unilateral(2:end-1) = 2*unilateral(2:end-1); % duplicar todos os termos excepto o primeiro e ultimo
 unilateral_ef = abs(unilateral)/sqrt(2);
 unilateral_db = 20*log10(unilateral_ef); % cada sinusoide tem de ser multiplicado por 1/sqrt(2) para o valor eficaz

% Frequencia
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
% Excluir primeira e ultima posicao
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
 
% peakB_id e da frequencia mais alta
if anterior < posterior
 peakB_id = peakA_id+1;
else
 peakB_id = peakA_id;
 peakA_id = peakA_id-1;
end
if iterations == 3
    
% estimativa do valor da frequencia, so e feito uma vez
omega_A = 2*pi*(peakA_id-1)*delta_f/fs;
omega_B = 2*pi*(peakB_id-1)*delta_f/fs;
XA = unilateral(peakA_id);
UA = real(XA);              % parte real
VA = imag(XA);              % parte imaginaria
XB = unilateral(peakB_id);
UB = real(XB);              % parte real
VB = imag(XB);              % parte imaginaria

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
% realizar a media dos espectros
unilateral_total = unilateral_total/ns;
unilateral_total_db = 20*log10(unilateral_total);

% RMS
valor_rms = 0;
for i = 1:ns
valor_rms = sqrt(sum(u0(i,1:end-N_discarded).^2)/N)+valor_rms;
end
valor_rms = valor_rms/ns;

% noise and distortion RMS
% este valor e muito dependente do espalhamento, depende da freq 35-40 db
nd_rms = sqrt(valor_rms^2-(unilateral_total(peak_id))^2);

% SINAD
sinad = valor_rms/nd_rms;
sinad_db = 20*log10(sinad);

% ENOB
enob = (sinad-1.76)/6.02;

% Plots
fig = figure(1);
s1 = sprintf('Frequencia estimada do sinal: %f [Hz] ', freq_estimada);
s2 = sprintf('Valor eficaz do sinal: %f [V] \n', valor_rms);
s3 = sprintf('Valor eficaz do ruido: %f [V] ', nd_rms);
s4 = sprintf('SINAD: %f [dB] \n', sinad_db);
s5 = sprintf('ENOB: %f bits ', enob);
string = sprintf('%s%s%s%s%s',s1,s2,s3,s4,s5);
subplot(2,1,1), plot(t(1:6*floor(pontos_p_periodo)),u(1:6*floor(pontos_p_periodo)));
title(string); % 5 periodos do sinal mudou-se e 3
xlim([0 t(6*floor(pontos_p_periodo))]);
xlabel('Tempo [s]');
ylabel('Amplitude [V]');
subplot(2,1,2), plot(freqs,unilateral_total_db);
xlabel('Frequencia [Hz]');
ylabel('Ampltitude [dB]');