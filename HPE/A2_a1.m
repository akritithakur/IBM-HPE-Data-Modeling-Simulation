%% HP data - Growth Analysis

clear;
close all;
clc;

%% Data extraction from .dat file
data = load('hewpac.dat');

global revenue;
global human_resourse;
global year;

revenue = data(:,1);        % Company Revenue
human_resourse = data(:,2); % Number of Employee (Human Resourse)
year = data(:,3);           % Year

%% PART1 - Curve Fitting

% Variable Declaration
T = length(year);
dt = 0.01;
n = T/dt;
revenue_fit = zeros(1,n);   % Logistic best fit for lower
hr_fit = zeros(1,n);
k = zeros(1,n);             % Changing Carrying capacity
time = zeros(1,n);
time(1) = 1;

% Parameter for Carrying Capacity
lambda_K = 0.1775;
eta_K = 0.0000000000001;
k(1) = 1e6;

% Parameters for Annual Revenue
lambda_R = 1;
revenue_fit(1) = 1e4;

for t=1:n-1
    k(t+1) = k(t) + (dt*(lambda_K*k(t)*(1-(eta_K*k(t)))));
    revenue_fit(t+1) = revenue_fit(t) + (dt*(lambda_R*revenue_fit(t)*(1-(revenue_fit(t)/k(t+1)))));
    time(t+1) = time(t)+dt;
end

% (Year-Revenue)
loglog(year,revenue,'o','color','b');   
hold on;

% Plot of Revenue & cummulutive revenue versus time on log--log scale
loglog(time,revenue_fit,'color','b');
title('Revenue Vs Time in log-log scale');
xlabel('Time');
ylabel('Momney in Million $');
legend('Revenue');
figure;

% (Year-Human Resourse)
loglog(year,human_resourse,'o','color','k');
hold on;

% Parameters for human resourse
lambda_H = 0.3;
alpha_H = 0.8;
eta_H = 0.0001;
hr_fit(1) = 5;

for t=1:n-1
    hr_fit(t+1) = hr_fit(t) + (lambda_H*hr_fit(t)*(1-eta_H*(power(hr_fit(t),alpha_H)))*dt);
    time(t+1) = time(t)+dt;
end

loglog(time,hr_fit,'color','k');
title('Human Resourse Vs Time In log-log Scale');
xlabel('Time');
ylabel('Numbe rof Employees');
legend('Human Resource');

index = [6,7,8,12,14,16,15,17];

disp('Human Resourse Missing Values');
for i=1:length(index)
    disp(strcat(strcat(num2str(index(i)),' : '),num2str(hr_fit(1+ ((index(i)-1)/dt)))));
end

disp('Revenue Missing values:');
for i=1:length(index)
    disp(strcat(strcat(num2str(index(i)),' : '),num2str(revenue_fit(1+ ((index(i)-1)/dt)))));
end