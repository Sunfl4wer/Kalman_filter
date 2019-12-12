clear,clc;
variance1 = 0.01;
variance2 = 120;
T = 0.01;
x = [30;10;0;1;0;0];
no_sample = 200;
Rmu = variance1;
Rv = variance2;
I = eye(6);
A = [[1,T,0,0,0,0];...
     [0,1,0,0,0,0];...
     [0,0,1,T,0,0];...
     [0,0,0,1,0,0];...
     [0,0,0,0,1,T];...
     [0,0,0,0,0,1]];
 
B = [[(T^2)/2,0,0];...
     [T,0,0];...
     [0,(T^2)/2,0];...
     [0,T,0];...
     [0,0,(T^2)/2];...
     [0,0,T]];
 
 H = [[1,0,0,0,0,0];...
      [0,0,1,0,0,0];...
      [0,0,0,0,1,0]];

 Ry = cov(x');
possition_no_noise = zeros(no_sample,3);
possition_with_noise = zeros(no_sample,3);
possition_kalman = zeros(no_sample,3);
x_no_noise = x;
x_kalman = x;
for i = 1:no_sample
     a = randn(3,1)*sqrt(variance1);
     noise_measure = randn(1)*sqrt(variance2)*ones(3,1);
     
     x_no_noise = A*x_no_noise+B*a;
     coordinate_no_noise = H*x_no_noise;

     coordinate_measure = coordinate_no_noise + noise_measure;
     
     x_kalman = A*x_kalman;
     %coordinate_kalman = H*x_kalman;  
     Ry = (A*Ry*(A')) + (B*Rmu*(B'));
     Rw = H*Ry*(H')+Rv;
     K = Ry*(H')/(Rw);
     x_kalman = x_kalman + K*(coordinate_measure-H*x_kalman);
     coordinate_kalman = H*x_kalman;
     Ry = (I-K*H)*Ry;
     for k = 1:3
         possition_no_noise(i,k) = coordinate_no_noise(k);
         possition_with_noise(i,k) = coordinate_measure(k);
         possition_kalman(i,k) = coordinate_kalman(k);
     end
end

X_no_noise = zeros(1,no_sample);
Y_no_noise = zeros(1,no_sample);
Z_no_noise = zeros(1,no_sample);
X_with_noise = zeros(1,no_sample);
Y_with_noise = zeros(1,no_sample);
Z_with_noise = zeros(1,no_sample);
X_kalman = zeros(1,no_sample);
Y_kalman = zeros(1,no_sample);
Z_kalman = zeros(1,no_sample);
for i = 1:no_sample
    for k = 1:3
        if(k==1)
            X_no_noise(i) = possition_no_noise(i,k);
            X_with_noise(i) = possition_with_noise(i,k);
            X_kalman(i) = possition_kalman(i,k);
        elseif(k==2)
            Y_no_noise(i) = possition_no_noise(i,k);
            Y_with_noise(i) = possition_with_noise(i,k);
            Y_kalman(i) = possition_kalman(i,k);
        else
            Z_no_noise(i) = possition_no_noise(i,k);
            Z_with_noise(i) = possition_with_noise(i,k);
            Z_kalman(i) = possition_kalman(i,k);
        end
    end
end
error = sqrt(X_no_noise.^2+Y_no_noise.^2+Z_no_noise.^2);
subplot(1,3,1);
plot(1:1:no_sample,X_no_noise);
hold on
plot(1:1:no_sample,X_with_noise);
hold on
plot(1:1:no_sample,X_kalman);
legend('x','x measure','x kalman');
hold off

subplot(1,3,2);
plot(1:1:no_sample,Y_no_noise);
hold on
plot(1:1:no_sample,Y_with_noise);
hold on
plot(1:1:no_sample,Y_kalman);
legend('y','y measure','y kalman');
hold off

subplot(1,3,3);
plot(1:1:no_sample,Z_no_noise);
hold on
plot(1:1:no_sample,Z_with_noise);
hold on
plot(1:1:no_sample,Z_kalman);
legend('z','z measure','z kalman');
hold off