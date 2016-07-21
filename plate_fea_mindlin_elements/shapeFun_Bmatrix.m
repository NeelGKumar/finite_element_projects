clc; clear all; close all; format compact; format long

% sk = (-1:.1:1);
% tk = (-1:.1:1);
% n = size(sk,2); 
% 
% N1 = zeros(n); N2 = zeros(n); N3 = zeros(n);
% N4 = zeros(n); N5 = zeros(n); N6 = zeros(n);
% N7 = zeros(n); N8 = zeros(n); N9 = zeros(n);
% 
% ns = 0; ts = 0;
% 
% for s = -1:.1:1
%     
%     ns = ns+1;
%     
%     for t = -1:.1:1
%         
%         ts = ts+1;
%         
%         N9(ns,ts) = (1-s^2)*(1-t^2);
%         N5(ns,ts) = 0.5*(1-s^2)*(1-t)-0.5*N9(ns,ts);
%         N6(ns,ts) = 0.5*(1+s)*(1-t^2)-0.5*N9(ns,ts);
%         N7(ns,ts) = 0.5*(1-s^2)*(1+t)-0.5*N9(ns,ts);
%         N8(ns,ts) = 0.5*(1-s)*(1-t^2)-0.5*N9(ns,ts);
%         N1(ns,ts) = 0.25*(1-s)*(1-t)-0.5*(N5(ns,ts)+N8(ns,ts))-0.25*N9(ns,ts);
%         N2(ns,ts) = 0.25*(1+s)*(1-t)-0.5*(N5(ns,ts)+N6(ns,ts))-0.25*N9(ns,ts);
%         N3(ns,ts) = 0.25*(1+s)*(1+t)-0.5*(N6(ns,ts)+N7(ns,ts))-0.25*N9(ns,ts);
%         N4(ns,ts) = 0.25*(1-s)*(1+t)-0.5*(N7(ns,ts)+N8(ns,ts))-0.25*N9(ns,ts);
%         
%     end
%     
%     ts = 0;
%     
% end

% subplot(3,3,1); surf(sk,tk,N1); 
% subplot(3,3,2); surf(sk,tk,N2);
% subplot(3,3,3); surf(sk,tk,N3);
% subplot(3,3,4); surf(sk,tk,N4); 
% subplot(3,3,5); surf(sk,tk,N5);
% subplot(3,3,6); surf(sk,tk,N6);
% subplot(3,3,7); surf(sk,tk,N7); 
% subplot(3,3,8); surf(sk,tk,N8);
% subplot(3,3,9); surf(sk,tk,N9);


syms s t;

N9 = (1-s^2)*(1-t^2);
N5 = 0.5*(1-s^2)*(1-t)-0.5*N9; N6 = 0.5*(1+s)*(1-t^2)-0.5*N9;
N7 = 0.5*(1-s^2)*(1+t)-0.5*N9; N8 = 0.5*(1-s)*(1-t^2)-0.5*N9;
N1 = 0.25*(1-s)*(1-t)-0.5*(N5+N8)-0.25*N9;
N2 = 0.25*(1+s)*(1-t)-0.5*(N5+N6)-0.25*N9;
N3 = 0.25*(1+s)*(1+t)-0.5*(N6+N7)-0.25*N9;
N4 = 0.25*(1-s)*(1+t)-0.5*(N7+N8)-0.25*N9;

%Derivatives

% N1s = diff(N1, s)
% N1t = diff(N1, t)
% N2s = diff(N2, s)
% N2t = diff(N2, t)
% N3s = diff(N3, s)
% N3t = diff(N3, t)
% N4s = diff(N4, s)
% N4t = diff(N4, t)
% N5s = diff(N5, s)
% N5t = diff(N5, t)
% N6s = diff(N6, s)
% N6t = diff(N6, t)
% N7s = diff(N7, s)
% N7t = diff(N7, t)
% N8s = diff(N8, s)
% N8t = diff(N8, t)
% N9s = diff(N9, s)
% N9t = diff(N9, t)

N1s = t/4 + (s*(t^2 - 1))/2 - (s*(t - 1))/2 - t^2/4;
N1t = s/4 - t*(s/2 - 1/2) + (t*(s^2 - 1))/2 - s^2/4;
N2s = (s*(t^2 - 1))/2 - t/4 - (s*(t - 1))/2 + t^2/4; 
N2t = t*(s/2 + 1/2) - s/4 + (t*(s^2 - 1))/2 - s^2/4;
N3s = t/4 + (s*(t^2 - 1))/2 + (s*(t + 1))/2 + t^2/4;
N3t = s/4 + t*(s/2 + 1/2) + (t*(s^2 - 1))/2 + s^2/4;
N4s = (s*(t^2 - 1))/2 - t/4 + (s*(t + 1))/2 - t^2/4;
N4t = (t*(s^2 - 1))/2 - t*(s/2 - 1/2) - s/4 + s^2/4;
N5s = s*(t - 1) - s*(t^2 - 1);
N5t = s^2/2 - t*(s^2 - 1) - 1/2;
N6s = 1/2 - t^2/2 - s*(t^2 - 1);
N6t = - 2*t*(s/2 + 1/2) - t*(s^2 - 1);
N7s = - s*(t^2 - 1) - s*(t + 1);
N7t = 1/2 - s^2/2 - t*(s^2 - 1);
N8s = t^2/2 - s*(t^2 - 1) - 1/2;
N8t = 2*t*(s/2 - 1/2) - t*(s^2 - 1);
N9s = 2*s*(t^2 - 1);
N9t = 2*t*(s^2 - 1);


