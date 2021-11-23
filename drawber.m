clc
clear;
close all;

data_rate1 = [9.62,13.91,16.46,17.84];
data_rate2 = [7.55,7.8,8.05,8.3,8.8];
BER1 = [9.3e-3, 1.03e-2, 1.68e-2, 3.17e-2];
BER2 = [1.018e-5, 3.55e-4, 8.84e-4, 3.3e-3,0.4047];
figure,semilogy(data_rate1,BER1,'O--'),hold on,semilogy(data_rate2,BER2,'O--'),hold on,
semilogy([7 18],[0.024 0.024],'k--');xlim([7 18]);
xlabel('Data rate (kb/s)');legend('PAM4','OOK');
ylabel('BER');
grid off;