dispositivo = daqlist;
placa_id = dispositivo.DeviceID;
placa_vendor = dispositivo.VendorID;
M = 100; % numero de segmentos
sample_rate = 5000; % fs da daq
number_samples = 100; % numero de amostras daq
dq = daq(dispositivo.VendorID);
dq.Rate = sample_rate;
addinput(dq, placa_id, "ai0", "Voltage");
addinput(dq, placa_id, "ai1", "Voltage");
for a = 1:M
[scanData, timeStamp] = read(dq,number_samples,"OutputFormat","Matrix");
uin0(a,:) = scanData(:,1).';
uout0(a,:) = scanData(:,2).';
end
fs = 1/(timeStamp(2)-timeStamp(1));
t0 = timeStamp.';
%
N0 = 1000; % n de amostras
fase = 0;
N0 = number_samples;
delta_t = (t0(end)-t0(1))/(N0-1);
fs = 1/(timeStamp(2)-timeStamp(1));
% necessario fs, N0, t0, uin0, uout0, M
det_freq = 1; % a frequencia se deve ser determinada uma vez
for inout = 1:M
    N_discarded = 0;
    iteration = 1;
    while iteration > 0
    N = N0 - N_discarded;
    uin = uin0(inout,1:end-N_discarded);
    uout = uout0(inout,1:end-N_discarded);

    delta_f = fs/N;
    freqs = fs*(0:(N/2))/N;

    % Espectro da entrada
    bilateralin = fft(uin)/N; %FFT de u abs(fft(u)/N)
    unilateralin = bilateralin(1:floor(N/2)+1);
    unilateralin(2:end-1) = 2*unilateralin(2:end-1); % duplicar todos os termos excepto o primeiro e ultimo
    unilateral_efin = abs(unilateralin)/sqrt(2);
    unilateral_dbin = 20*log10(unilateral_efin); % cada sinusoide tem de ser multiplicado por 1/sqrt(2) para o valor eficaz
    % Espectro da saida

    % Frequencia
    bilateralout = fft(uout)/N; %FFT de u abs(fft(u)/N)
    unilateralout = bilateralout(1:floor(N/2)+1);
    unilateralout(2:end-1) = 2*unilateralout(2:end-1); % duplicar todos os termos excepto o primeiro e ultimo
    unilateral_efout = abs(unilateralout)/sqrt(2);
    unilateral_dbout = 20*log10(unilateral_efout); % cada sinusoide tem de ser multiplicado por 1/sqrt(2) para o valor eficaz

    %Frequencia
    if det_freq == 1
    peak = -inf;
    peak_id = 0;

    for i = 1:(N/2-1)
        if peak < unilateral_dbin(i)
        peak_id = i;
        peak = unilateral_dbin(i);
        end
    end

    peakA_id = peak_id;
    % Excluir primeira e ultima posicao
    if peakA_id == 1
    anterior = -inf;
    posterior = unilateral_dbin(peakA_id+1);
    elseif peakA_id == (N/2)+1
    posterior = -inf;
    anterior = unilateral_dbin(peakA_id-1);
    else
    anterior = unilateral_dbin(peakA_id-1);
    posterior = unilateral_dbin(peakA_id+1);
    end

    % peakB_id e da frequencia mais alta
    if anterior < posterior
    peakB_id = peakA_id+1;
    else
    peakB_id = peakA_id;
    peakA_id = peakA_id-1;
    end
    
    %then the sinc thingy
    omega_A = 2*pi*(peakA_id-1)*delta_f/fs;
    omega_B = 2*pi*(peakB_id-1)*delta_f/fs;
    XA = unilateralin(peakA_id);
    UA = real(XA); % parte real
    VA = imag(XA); % parte imaginaria
    XB = unilateralin(peakB_id);
    UB = real(XB); % parte real
    VB = imag(XB); % parte imaginaria
    Kopt = ((VB-VA)*sin(omega_A)+(UB-UA)*cos(omega_A))/(UB-UA);
    ZA = VA*(Kopt-cos(omega_A))/sin(omega_A)+UA;
    ZB = VB*(Kopt-cos(omega_B))/sin(omega_B)+UB;
    freq_estimada = (fs/(2*pi))*acos((ZB*cos(omega_B)-ZA*cos(omega_A))/(ZB-ZA));

     % info para fazer nova FFT
    pontos_p_periodo = fs/freq_estimada;
    N_discarded = round(rem(N,pontos_p_periodo)); %pontos por discartar
    det_freq = 0;
    iteration = iteration + 1;
    end
    iteration = iteration - 1;
    end
    %rms for gain
    rms_in = sqrt(sum(uin(1:end).^2)/N);
    rms_out = sqrt(sum(uout(1:end).^2)/N);
    ganho(inout) = rms_out/rms_in;
    % Phase
    a_in = rms_in*sqrt(2);
    a_out = rms_out*sqrt(2); %sinusoidal signal
    phase_in = (180/pi)*asin(uin(1)/a_in);
    phase_out = (180/pi)*asin(uout(1)/a_out);
    phase(inout) = real(phase_out-phase_in);

end

% Histograma
[vbinsg , nbinsg] = histcounts(ganho);
[vbinsf , nbinsf] = histcounts(phase);
fig = figure(1);
subplot(2,1,1),HG = histogram(ganho);
ylabel('Histograma do ganho');
subplot(2,1,2),HF = histogram(phase);