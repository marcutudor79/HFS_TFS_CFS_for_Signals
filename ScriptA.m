clc
clear
close all

%declaring the fundamental frequency and period of the signal
f0 = 200;
T = 1/f0;

%defining the step for the 
pas = T/100;

%defining the time vector for 4 periods
t  = 0:pas:4*T;

%generating the signal's values
Name_Signal = name_signal(t,f0);

%plotting the name signal on 4 periods
figure(1)
subplot(4,1,1)
plot(t,Name_Signal)
axis([min(t),max(t),3,22])
xlabel("Time[s]")
ylabel("Amplitude")
title("Name Signal on 4 periods")

%calculating the fourier series
%determining the radial frequency
omega0 = 2 * pi * f0;

%computing the DC component
a0pe2 = 1/T * integral(@(t)name_signal(t,f0), 0,T);

%computing the first N coefficients
N = 100;

a = zeros(1,N);
b = zeros(1,N);

for k = 1:N
    a(k) = 2/T * integral(@(t)name_signal(t,f0).*cos(k*omega0*t),0,T);
    b(k) = 2/T * integral(@(t)name_signal(t,f0).*sin(k*omega0*t),0,T);
end

%adding a threshold so that the 0 is 0
thr = 10^-5;

if abs(a0pe2) < thr
    a0pe2 = 0;
end

%if the values are very small, they are approximated to 0
for k = 1:N
    if abs(a(k)) < thr
        a(k) = 0;
    elseif abs(b(k)) < thr
        b(k) = 0;
    end
end

omega = (1:N) * omega0;

%plotting the Amplitude Spectrum of TFS
subplot(4,1,2)
title("Amplitude Spectrum TFS")
hold on
stem(0,a0pe2)
stem(omega,a)
stem(omega,b)
xlabel("Frequency [rad/s")
ylabel("Amplitude of a_k and b_k")
hold off

%computing the Harmonic Fourier Series
A = zeros(1,N+1);
phi = zeros(1,N+1);

A(1) = abs(a0pe2);

if a0pe2 >= 0
    phi(1) = 0;
else
    phi(1) = -pi;
end

for k = 1:N
    A(k+1) = sqrt(a(k)^2 + b(k)^2);

    if a(k) == 0 && b(k) == 0
        phi(k+1) = 0;
    else 
        if a(k) >= 0
            phi(k+1) = -atan(b(k)/a(k));
        else 
            phi(k+1) = -atan(b(k)/a(k))-pi;
        end
    end
end

omega = (0:N) .* omega0;

%plotting the amplitude and phase spectrum
subplot(4,1,3)
stem(omega,A)
xlabel('Frequency [rad/s]')
ylabel('Amplitudes A_k')
title('Amplitude spectrum HFS')

subplot(4,1,4)
stem(omega,phi)
xlabel('Frequency [rad/s]')
ylabel('Phase \phi_k')
title('Phase spectrum HFS')

%computing the Complex Fourier Series

C = zeros(1,2*N+1);

for n = -N:N
    C(n+N+1) = 1/T * integral(@(t)(name_signal(t,f0).*exp(-1j*n*omega0*t)),0,T);

    re = real(C(n+N+1));
    im = imag(C(n+N+1));

    if abs(re) < thr
        re = 0;
    end

    if abs(im) < thr
        im = 0;
    end

    C(n+N+1) = re + 1j*im;
end

%plotting the CFS
figure(2)
subplot(2,1,1)
stem((-N:N)*omega0,abs(C))
xlabel("Freqency [rad\s]")
ylabel("Amplitudes |C(n\omega_0)|")
title("Amplitudes spectrum for CFS")

subplot(2,1,2)
stem((-N:N)*omega0, angle(C));
xlabel("Frequency [rad\s]")
ylabel("Phases arg\{C (n\omega_0)\}")
title("Phase spectrum for CFS")


%Total power of the signal
Pt = 1/T*integral(@(t)(abs(name_signal(t,f0)).^2),0,T);
disp(['Total power of the name signal = ', int2str(Pt) , ' W'])

P99 = 0.99 * Pt;

fstart = 0;

Pn = A(1) ^ 2;
k = 1;

while Pn < P99
    if Pn <= 0
        fstart = k;
    end
    Pn = Pn + A(k+1)^2/2;
    fstop = k;
    k = k + 1;
end

disp(['Energy Band = [' int2str(fstart*f0) ',' int2str(fstop*f0) ']'])