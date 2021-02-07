%% IBM data - Growth Analysis

clear;
close all;
clc;

%% Data extraction from .dat file
data = load('growth_ibm.dat');

global revenue;
global human_resourse;
global year;

revenue = data(:,1);        % Company Revenue
human_resourse = data(:,2); % Number of Employee (Human Resourse)
year = data(:,3);           % Year

%% PART1

% Lower Curve (Year-Revenue)
loglog(year,revenue,'o','color','b');   
hold on;
% Upper Curve (Year-Cummulative Revenue)
cumulutive = cumsum(revenue);
loglog(year,cumulutive,'o','color','k'); 
hold on;

% Variable Declaration
T = length(year);
dt = 0.01;
n = T/dt;
revenue_fit = zeros(1,n);   % Logistic best fit for lower
cumrevenue_fit = zeros(1,n);% Logistic best fit for upper
time = zeros(1,n);
time(1) = 1;

% Information from Paper
lambda_R = 0.145;
eta_R = 0.00001;
eta_C = 0.0000004;
disp('Initial Valur for Revenue fit curve ');
revenue_fit(1) = initial_value(lambda_R,eta_R,revenue,year,length(revenue));
disp('Initial Valur for Cumulative Revenue fit curve ');
cumrevenue_fit(1) = initial_value(lambda_R,eta_C,cumulutive,year,length(cumulutive));
for t=1:n-1
    revenue_fit(t+1) = revenue_fit(t) + (lambda_R*revenue_fit(t)*(1-eta_R*revenue_fit(t))*dt);
    cumrevenue_fit(t+1) = cumrevenue_fit(t) + (lambda_R*cumrevenue_fit(t)*(1-eta_C*cumrevenue_fit(t))*dt);
    time(t+1) = time(t)+dt;
end

% Plot of Revenue & cummulutive revenue versus time on log--log scale
loglog(time,revenue_fit,'color','b');
hold on;
loglog(time,cumrevenue_fit,'color','k');
title('Revenue & Cumulutive Revenue Vs Time in log-log scale');
xlabel('Time');
ylabel('Momney in Million $');
legend('Revenue','Cumulative Revenue');
figure;

loglog(year,human_resourse,'o','color','k');
hold on;

% Variable Declaration
hr_fit = zeros(1,n);

% Information from paper
lambda_H = 0.09;
alpha = 1;
eta_H = 0.000002;
disp('Initial Valur for Human Resourse fit curve ');
hr_fit(1) = initial_value(lambda_H,eta_H,human_resourse,year,10);

for t=1:n-1
    hr_fit(t+1) = hr_fit(t)+lambda_H*hr_fit(t)*(1-eta_H*hr_fit(t))*dt;
    time(t+1) = time(t)+dt;
end

loglog(time,hr_fit,'color','k');
title('Human Resourse Vs Time In log-log Scale');
xlabel('Time');
ylabel('Numbe rof Employees');
legend('Human Resource');
figure;

%% PART2 - Non-linearlity effact (exponential growth stops) - lambda slop

plot(time,log(revenue_fit),'color','b');
hold on;
plot(time,lambda_R*time + log(revenue_fit(1)),'--','color','r');
title('Revenue Vs Time');
xlabel('Time');
ylabel('Revenue in log scale');
hold on;
plot(time,log(cumrevenue_fit),'color','k');
hold on;
plot(time,lambda_R*time + log(cumrevenue_fit(1)),'--','color','r');
legend('Revenue','Cumulative Revenue');
figure;

plot(time,log(hr_fit),'color','k');
hold on;
plot(time,lambda_H*time + log(hr_fit(1)),'--','color','r');
title('Human Resourse Vs Time ');
xlabel('Time');
ylabel('Human Resourse in log scale');
legend('Human Resource');
figure;

%% PART 3 - dx/dt versus x (Phase Portrait)
x = [0:0.1:100000];
y = lambda_R.*x.*(1-eta_R.*x);
plot(x,y,'color','b');
title('Derivative of Revenue Vs Revenue (Downward Parabola)');
xlabel('Revenue');
ylabel('Derivatie of Revenue');
legend('Derivative of Revenue');
figure;

x = [1000:0.2:1000000];
y = lambda_H.*x.*(1-eta_H.*x);
plot(x,y,'color','k');
title('Derivative of Human Resource Vs Human Resource (Downward Parabola)');
xlabel('Resource');
ylabel('Derivatie of Human Resource');
legend('derivative of Human Resource');
figure;

%% PART 4 - Carrying Capacity & Net Annual Earning

pre = revenue_fit(length(revenue_fit)-1);
curr = revenue_fit(length(revenue_fit));
while curr-pre>0.0000001
    pre = curr;
    curr = pre + dt*(lambda_R*pre*(1-(eta_R*pre)));
end
disp(strcat('Carrying Capacity of Revenue is ',num2str(curr)));
disp('(with Accuracy of 0.0000001)');


pre = hr_fit(length(hr_fit)-1);
curr = hr_fit(length(hr_fit));
while curr-pre>0.00000001
    pre = curr;
    curr = pre + dt*(lambda_R*pre*(1-(eta_R*pre)));
end
disp(strcat('Carrying Capacity of Human Resourse is ',num2str(curr)));
disp('(with Accuracy of 0.00000001)');

data = load('prof_ibm.dat');

year = data(:,1);
profit = data(:,2);
plot(year,profit,'color','k');
title('Net Annual Earning Vs TIme');
xlabel('Time');
ylabel('Profit');
legend('Money');
