clc
close all
clear all
drawnow
%
tau_low=50/(3e8);
tau_high=150/(3e8);
fc=2.5e9;
%
MaxRealizationsNum=1e6;
% PHYplot.m
d_tau=tau_high-tau_low;
H=zeros(MaxRealizationsNum,1);
%Hi=Hr;
%
for realizationIndex=1:MaxRealizationsNum
taus=tau_low+d_tau*rand(10,1);
Hinst=sum(exp(1i*2*pi*taus*fc));
H(realizationIndex)= (Hinst);
%Hr(realizationIndex)=real(Hinst);
%Hi(realizationIndex)=imag(Hinst);
end
%
[hr,edges]=histcounts(real(H),round(sqrt(MaxRealizationsNum)));
[hi]=histcounts(imag(H),edges);
m=mean(H);
sr=std(real(H));
si=std(imag(H));
pr_core=1/(2)*erfc(-(edges-real(m))/(sr*sqrt(2)));
%pr_core = normcdf(edges,real(m),sr); % Requres statistics toolbox
pr=diff(pr_core);
pi_core=1/(2)*erfc(-(edges-imag(m))/(si*sqrt(2)));
% pi_core = normcdf(edges,imag(m),si); % Requres statistics toolbox
pi=diff(pi_core);
R=corrcoef(real(H),imag(H))
% Display
subplot(311)
plot(edges(1:end-1),100*hr/MaxRealizationsNum,'b')
hold on
plot(edges(1:end-1),100*pr,'b-')
hold off
grid on
ylabel('Probability [%]')
xlabel('Real(H)')
legend('Empirical real','Theoretical real')
title('Real component - empirical distribution & normal fit')
subplot(312)
plot(edges(1:end-1),hi/MaxRealizationsNum,'r')
hold on
plot(edges(1:end-1),pi,'r-')
hold off
legend('Empirical imaginary','Theoretical Imaginary')
grid on
ylabel('Probability [%]')
xlabel('Imag(H)')
title('Imaginary component - empirical distribution & normal fit')
subplot(313)
plotmatrix(real(H),imag(H))
xlabel('real(H)')
ylabel('imag(H)')
title({'Scatter plot of real and imaginary componnets of H', ['Corellation coeficient ' num2str(R(1,2)) ]})
