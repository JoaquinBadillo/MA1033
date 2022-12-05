pkg load symbolic

clc; clear all; close all;

% Parameters (SI units)
m = 80; 
g = 9.81;
beta = 12;
h0 = 200;
l = 60;
k = 15;

% CODE
syms v1(t) s1(t) s2(t) s(t)

v1(t) = g*m/beta * (exp(-beta/m*t) - 1);
s1(t) = -g*(m^2)/(beta^2) * exp(-beta/m * t) - g*m/beta * t + g*m^2/beta^2 + h0;
t1 = solve(s1(t) == h0-l, t, "Real",true,"PrincipalValue",true);

lambda = - beta/(2*m);
mu = sqrt(abs(beta^2-4*m*k))/(2*m);
 
a11 = cos(mu*t1);
a12 = sin(mu*t1);
a21 = lambda*cos(mu*t1) - mu*sin(mu*t1);
a22 = lambda*sin(mu*t1) + mu*cos(mu*t1);
A = [a11, a12; a21, a22];
 
b1 = m*g/(k*exp(lambda*t1));
b2 = v1(t1)/exp(lambda*t1);
b = [b1; b2];
 
C = [a11, a12; a21, a22] \ [b1; b2];
C1 = simplify(C(1));
C2 = simplify(C(2));
s2(t) = (C1*cos(mu*t) + C2*sin(mu*t))*exp(lambda*t) + h0 - l - m*g/k;

% PLOT

s(t) = piecewise(t <= t1,s1(t),t > t1,s2(t));
figure
hold on
grid on
title("Position vs. time")
xlabel("Time (s)")
ylabel("Position (m)")
ezplot(s(t), [0, 80], "LineWidth",1.2)