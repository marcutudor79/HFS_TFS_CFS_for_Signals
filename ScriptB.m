clc
clear
close all

f0 = 200;
T = 1/f0;

pas = T/100;


t=0:pas:5*T;

fc = f(t,f0);

%plotting the signal for 5 periods
figure(1)
plot(t,fc)
xlabel("Time [S]")
ylabel("Amplitude")

%computing TFS
a0pe2 = 1/T * integral(@(t) f(t,f0), 0, T);
omega0 = 2 * pi * f0;

%computing the TFS a and b coefficients
N = 100;
a = zeros(1,N);
b = zeros(1,N);

for k = 1:N
    a(k) = 2/T * integral(@(t)f(t,f0).*cos(k*omega0*t),0,T);
    b(k) = 2/T * integral(@(t)f(t,f0).*sin(k*omega0*t),0,T);
end

%adjusting the values of the coefficients in order to make the really
%low ones 0

%defining the threshold value
thr = 10^-5;

for k = 1:N
    if a(k) < thr
        a(k) = 0;
    end
    
    if b(k) < thr
        b(k) = 0;
    end
end

%computing the HFS

A = zeros(1,N+1);
phi = zeros(1,N+1);

A(1) = a0pe2;

if a0pe2 >= 0
    phi(1) = 0;
else
    phi(1) = -pi;
end

for k = 1:N
    A(k+1) = sqrt(a(k)^2+b(k)^2);
    
    if(a(k) == 0 && b(k) == 0)
        phi(k+1) = 0;
    end

    if a(k) > 0
        phi(k+1) = -atan(b(k)/a(k));
    else
        phi(k+1) = -atan(b(k)/a(k)) - pi;
        if phi(k+1) > pi || phi(k+1) < pi
            phi(k+1) = -atan(b(k)/a(k)) +pi;
        end
    end
end

%plotting the HFS
figure(1)
subplot(2,1,1)
stem((0:N)*omega0,A)
title ("Spectrum of amplitudes for HFS")
xlabel("Frequency [rad/s]")
ylabel("Amplitudes of A_k")
axis([-1,N*omega0,0,0.4])

subplot(2,1,2)
stem((0:N)*omega0,phi)
title("Spectrum of phases for HFS")
xlabel("Freqeuncy [rad/s]")
ylabel("Phase \phi_k")
axis([-1,N*omega0,-pi,pi])

%computing the CFS

C = zeros(1,2*N+1);

for n = -N:N
    
    C(n+N+1) = 1/T * integral(@(t)(f(t,f0).*exp(-1j*n*omega0*t)),0,T);

    re = real(C(n+N+1));
    im = imag(C(n+N+1));
  
    if re <= thr
        re = 0;
    end
    if im <= thr
        im = 0;
    end

    C(n+N+1) = re + 1j*im;
end

%plotting the CFS
figure(2)
subplot(2,1,1)
stem((-N:N)*omega0,abs(C))
title("Spectrum of Amplitudes for CFS")
xlabel("Frequency [rad/s]")
ylabel("Amplitudes |C(n\omega_0)|")
axis([-N*omega0,N*omega0,0,0.3])

subplot(2,1,2)
stem((-N:N)*omega0, angle(C));
xlabel("Frequency [rad\s]")
ylabel("Phases arg\{C (n\omega_0)\}")
title("Phase spectrum for CFS")
axis([-N*omega0,N*omega0,-pi,pi])

