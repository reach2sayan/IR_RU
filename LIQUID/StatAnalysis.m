function [val sterr]=StatAnalysis()

%loading the time series that we want to do statistical analysis.
s1=load('E0');
s2=load('EK');
s=s1+s2;
time=(1:length(s)).*3; %eahc timestep in MD is 3 fs
plot(time,s)
pause
a=s(200:end)/54; %a stores all the values in the time series
t=length(a) %this is the total number of data points (total lenght of chain or time series)

disp('loading completed');
tBs=[];
p=[];
stds=[];
ro=[];  % firstAutoCorrelation vector for different step time (gap time)

figure(1);
plot(a);

  
pause

for i=1:t/5

tB=i; %this is the block size. The average of each tB number of data points will be averaged. All these successive averages will be stored in a
        %vector called A_b
        
nB=t/tB; %this number indicates the number of times a block
A_b=zeros(nB,1);

for j=1:nB
    
    A_b(j,1)=mean(a((j-1)*tB+1:j*tB-1));
       
end

tBs=[tBs tB];
stds=[stds std(a(1:tB:t))];
b=a(1:tB:t);
covar=cov(b(1:length(b)-1),b(2:length(b)));
ro=[ro covar(1,2)/covar(1,1)];
p=[p tB*var(A_b)/var(a)];

end

figure(2);
plot(tBs(:),p(:))
xlabel('block time');
ylabel('decay time')

figure(3);
plot(tBs(:),stds(:))
xlabel('block time');
ylabel('standard deviation')


figure(4);
plot(tBs(:),ro(:))
xlabel('block time');
ylabel('first Autocorrelation')

pause

tB=25;

sterr=std(a(1:tB:t))/sqrt(length(a(1:tB:t)));
val=mean(a(1:tB:t));


end

