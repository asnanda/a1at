clc
clear all
close all
files = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/5e3Data/Polymer/');
files2 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/5e3Data/Monomer/');
files3 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/5e3Data/Polymer/');


%cumsum(files(i).y-mean(files(i).y));



% can generate eigenvector defining 92.17% of the data for monomer
A = files(1).si_cor_norm;
B = files2(1).si_cor_norm;
C = files3(1).si_cor_norm;

for i = 2:length(files)
    m = files(i).si_cor_norm;
    
    %m = m/max(m);
    %plot(m)
    for j = 1:length(m)
        %if sign(m(j)) == -1
            %m(j) = 0;
        %end
    end
    A = [A;m];
   
end

for i = 2:length(files2)
    n = files2(i).si_cor_norm;
    o = files3(i).si_cor_norm;
    B = [B;n];
    C = [C;o];
end


[pc,score,latent,tsquare] = princomp(A);



M = transpose(A)*A;
e = eig(M);

[V,D] = eig(M);

[lengthA,~] = size(A);

r = A*[V(:,512) V(:,511)];
s = B*[V(:,512) V(:,511)];
t = C*[V(:,512) V(:,511)];

hold on
figure(6)
scatter(r(:,1),r(:,2),'k')
scatter(s(:,1),s(:,2),'r')
scatter(t(:,1),t(:,2),'b')



for k = 1:lengthA
    alpha(k) = A(k,:)*V(:,512);
    beta(k) = A(k,:)*V(:,511);
    gamma(k) = A(k,:)*V(:,510);
end
figure(1)
hold on;
scatter(alpha,beta)

P = polyfit(alpha,beta,1);
x = linspace(-1000,1000,2000);
y = P(1)*x+P(2);
plot(x,y,'r-.');

%figure(2)
%hold on
for l = 1:2000
    spectra = x(l)*V(:,512)+y(l)*V(:,511);
    peak = find(ismember(spectra,max(spectra)));
    max_dist = 512-peak;
    dist = min(max_dist,peak)-1;
    rmsd(l) = sqrt((sum(spectra(peak-dist:peak)-flipud(spectra(peak:peak+dist))).^2)/dist);
    %plot(spectra(peak-dist:peak));
    %plot(flipud(spectra(peak:peak+dist)));
    
    %calculate Hsy
end

figure(4)
plot(x,rmsd)
%rmsd = sort(rmsd);
%plot(x,rmsd)

DataInv = 1.01*max(rmsd) - rmsd;
[Minima,MinIdx] = findpeaks(DataInv);

purespec1 = x(MinIdx(1))*V(:,512)+y(MinIdx(1))*V(:,511);
purespec2 = x(MinIdx(2))*V(:,512)+y(MinIdx(2))*V(:,511);

figure(3)
hold on
plot(diff(purespec1)/max(diff(purespec1)))
plot(diff(purespec2)/max(diff(purespec2)))




