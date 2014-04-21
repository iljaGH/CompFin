fid = fopen('task2.dat','r')
data =fscanf(fid,'%f');
fclose(fid);
[vals,bins]=hist(data,100);
plot(bins,vals/trapz(bins,vals),'*',bins,1/sqrt(2*pi)*exp(-bins.^2/2),'-');
legend('rejection sampling','p(x)');
xlabel('x');
ylabel('density');
pause;