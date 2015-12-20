clc
clear all
close all
files = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/all/');
files2 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/monomer/');
files3 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/4b12/');
files4 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/Mutants/4b12/333/');

break
%{
for i = 1:2:length(files)
    figure(i)
    hold on
    plot(files(i).x,files(i).y_cor_norm/max(files(i).y_cor_norm),'k',files(i+1).x-0.0691,files(i+1).y_cor_norm/max(files(i+1).y_cor_norm),'r');
    plot(files4(1).x+0.078,files4(1).y_cor_norm/max(files4(1).y_cor_norm),'b',files4(2).x-0.019,files4(2).y_cor_norm/max(files4(2).y_cor_norm));
    %legend(files(i).mutant,files(i).mutant,'306apo','3064b12');
    
end


for i = 1:length(files2)
    var_mon(i) = mad(files2(i).y_cor_norm);
    var_pol(i) = mad(files3(i).y_cor_norm);
    var_ratio(i) = var_mon(i)/var_pol(i);
    var_theta(i) = acos(dot(files2(i).y_cor_norm,files3(i).y_cor_norm)/(norm(files2(i).y_cor_norm,2)*norm(files3(i).y_cor_norm,2)));
end

hold on
scatter(struct2mat_mutant(files2,'mutant'),log(var_mon))
scatter(struct2mat_mutant(files3,'mutant'),log(var_pol),'r')
plot(struct2mat_mutant(files3,'mutant'),(var_ratio),'k')
plot(struct2mat_mutant(files3,'mutant'),(ones(length(files3))),'-g');
%scatter(struct2mat_mutant(files3,'mutant'),(var_theta),'m')

figure(2)
hold on

for i = 1:length(files2)
   ecdf((((files2(i).y_cor_norm))));
   ecdf((((files3(i).y_cor_norm))));

end
break

%cumsum(files(i).y-mean(files(i).y));
%}

% can generate eigenvector defining 92.17% of the data for monomer
A = files(1).si_cor_norm_di;
B = files2(1).si_cor_norm_di;
C = files3(1).si_cor_norm_di;
hold on
for i = 2:length(files)
    m = files(i).si_cor_norm_di - mean(files(i).si_cor_norm_di);
    
    %m = m/max(m);
  
 
    for j = 1:length(m)
        if sign(m(j)) == -1
            m(j) = 0;
        end
    end

    A = [A;m];
     %plot(m)
end

for i = 2:length(files2)
    n = files2(i).si_cor_norm_di - mean(files2(i).si_cor_norm_di);
    o = files3(i).si_cor_norm_di - mean(files3(i).si_cor_norm_di);
    B = [B;n];
    C = [C;o];
end



[pc,score,latent,tsquare] = princomp(A);



M = transpose(A)*A;
e = eig(M);

[V,D] = eig(M);

[lengthA,~] = size(A);

r = A*[pc(:,1) pc(:,2)];
s = B*[pc(:,1) pc(:,2)];
t = C*[pc(:,1) pc(:,2)];

hold on
figure(1)
%scatter(r(:,1),r(:,2),'k')
scatter(s(:,1),s(:,2),'r')
scatter(t(:,1),t(:,2),'b')


for k = 1:lengthA
    alpha(k) = A(k,:)*V(:,512);
    beta(k) = A(k,:)*V(:,511);
    gamma(k) = A(k,:)*V(:,510);
end
figure(2)
hold on;
scatter(alpha,beta)
labelpoints(alpha,beta,struct2mat_mutant(files,'mutant'),'NW') 
%sufficent to only label one with arrows.


P = polyfit(alpha,beta,1);
x = linspace(-20000,200000,200000);
y = P(1)*x+P(2);
plot(x,y,'r-.');


%figure(2)
%hold on
for l = 1:200000
    spectra = x(l)*V(:,512)+y(l)*V(:,511);
    peak = find(ismember(spectra,max(spectra)));
    max_dist = 512-peak;
    dist = min(max_dist,peak)-1;
    rmsd(l) = sqrt((sum(spectra(peak-dist:peak)-flipud(spectra(peak:peak+dist))).^2)/dist);
    %plot(spectra);
    %plot(flipud(spectra(peak:peak+dist)));
    
    %calculate Hsy
end

figure(4)
plot(x,rmsd)
%rmsd = sort(rmsd);
%plot(x,rmsd)

DataInv = 1.01*max(rmsd) - rmsd;
[Minima,MinIdx] = findpeaks(DataInv);

break


purespec1 = x(MinIdx(1))*V(:,512)+y(MinIdx(1))*V(:,511);
purespec2 = x(MinIdx(2))*V(:,512)+y(MinIdx(2))*V(:,511);

figure(3)

plot(files(1).x(1:511),diff(baseline_correct_mutant(files(1).x,-1*purespec2','absorbance','poly')))
%xlabel('Magnetic Field (mT)','FontName','arial','fontsize',11)
set(gca,'YTickLabels','','YTick','','LineWidth',2,'TickDir','out','TickLength',[0.03 0.02],'XMinorTick','on','XTick',344:2:352);
set(gca,'YColor',[1 1 1]);
%set(gca,'XColor',[1 1 1]);

%set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 40, 50], 'PaperUnits', 'centimeters','PaperSize', [21, 29.7])
set (gcf,'InvertHardcopy','off','Color',[1 1 1])
grid off
box off



figure(4)

plot(files(1).x(1:511),diff(baseline_correct_mutant(files(1).x,purespec1','absorbance','poly')))


   %aesthetics:
%xlabel('Magnetic Field (mT)','FontName','arial','fontsize',11)
set(gca,'YTickLabels','','YTick','','LineWidth',2,'TickDir','out','TickLength',[0.03 0.02],'XMinorTick','on','XTick',344:2:352);
set(gca,'YColor',[1 1 1]);
%set(gca,'XColor',[1 1 1]);

%set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 40, 50], 'PaperUnits', 'centimeters','PaperSize', [21, 29.7])
set (gcf,'InvertHardcopy','off','Color',[1 1 1])
grid off
box off

