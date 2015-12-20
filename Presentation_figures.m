clear all
close all

%individual plots
files = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/monomer/','/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/4b12/');
files2 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/monomer/');
files3 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/4b12Data/4b12/');
files4 = generate_struct('/Users/Scott/Desktop/UCL Lab/EPR/StandardCurveV2/');

figure(1)

subplot(2,1,1)
plot_mutant(files,'306','figure','4b12figure')
set(gca,'XColor',[1 1 1]);

subplot(2,1,2)
plot_mutant(files,'301','figure','4b12figure')

xlabel('Magnetic Field [mT]','Fontsize',14,'FontName','arial');




figure(2)

subplot(2,1,1)
plot_mutant(files,'296','figure','4b12figure')
set(gca,'XColor',[1 1 1]);

subplot(2,1,2)
plot_mutant(files,'333','figure','4b12figure')
%set(gca,'XColor',[1 1 1]);
xlabel('Magnetic Field [mT]','Fontsize',14,'FontName','arial');


break

subplot(4,1,3)
plot_mutant(files,'333','figure','4b12figure')
set(gca,'XColor',[1 1 1]);

subplot(4,1,4)
plot_mutant(files,'301','figure','4b12figure')

xlabel('Magnetic Field [mT]','Fontsize',14,'FontName','arial');




break
mutantnames = unique(struct2mat_mutant(files,'mutant'));
resnames =  {'E','S','H','N','E','T','S','L','V','T','I','T'};

for i = 1:12
    subplot(3,4,i)
    plot_mutant(files,num2str(mutantnames(i)),'figure','4b12figure');
    %set(gcf,'PaperUnits','centimeters','Position', [0, 0, 1000, 1000]);
    text(0.8,0.6,strcat(resnames(i),num2str(mutantnames(i)),'C') ,'Units','normalized','Linestyle','None','Fontsize',14,'FontName','arial');
    %print(strcat('/Users/Scott/Desktop/UCL Lab/Graphics/Summary2/',num2str(mutantnames(i)),'c'),'-depsc')
    %saveas(gcf,strcat('/Users/Scott/Desktop/UCL Lab/Graphics/Summary2/',num2str(mutantnames(i)),'c'),'epsc');
    %clf
end

[ax2,h2]=suplabel('d \chi" / dB','y',0.001);
set(h2,'FontSize',18,'FontName','arial');
xlabel('Magnetic Field [mT]','Fontsize',14,'FontName','arial');


break
figure(1)
plot_mutant(files2,'360','figure','4b12figure')

hold on
plot(files(22).x-0.06,files(22).y_cor_norm,'m','linewidth',3)

break

figure(1)
set (gcf,'InvertHardcopy','off','Color',[1 1 1])
bar(files(1).x,files(1).y_cor_norm)
axis([-inf inf -0.04 0.04])
figure(2)
set (gcf,'InvertHardcopy','off','Color',[1 1 1])
bar(files(2).x,files(2).y_cor_norm)
axis([-inf inf -0.04 0.04])


break
figure(1)
subplot(2,1,1)
plot_mutant(files,'306','figure','4b12figure')
set(gca,'XColor',[1 1 1]);
subplot(2,1,2)
plot_mutant(files,'301','figure','4b12figure')
xlabel('Magnetic Field[mT]','Fontsize',14,'FontName','arial');

figure(2)
subplot(2,1,1)
plot_mutant(files,'296','figure','4b12figure')
set(gca,'XColor',[1 1 1]);
subplot(2,1,2)
plot_mutant(files,'333','figure','4b12figure')
xlabel('Magnetic Field[mT]','Fontsize',14,'FontName','arial');


set (gcf,'InvertHardcopy','off','Color',[1 1 1])
set(gca, 'box', 'off','ytick', []  )
break
%set(gca,'XColor',[0 0 0]);
%axes([-1 513 min(files(1).y_cor_norm) max(files(1).y_cor_norm)])

%plot_mutant(files2,'43','figure','5e3figure')

figure(1)
hist(files(1).y_cor_norm,512);
figure(2)
hist(files(2).y_cor_norm,512,'g')
%axis in here
%
xlabel('Magnetic Field[mT]','Fontsize',14,'FontName','arial');

break

mutantnames = unique(struct2mat_mutant(files,'mutant'));
resnames =  {'E','S','H','N','E','T','S','L','V','T','I','T'};

for i = 1:12
    subplot(3,4,i)
    plot_mutant(files,num2str(mutantnames(i)),'figure','4b12figure');
    %set(gcf,'PaperUnits','centimeters','Position', [0, 0, 1000, 1000]);
    text(0.8,0.6,strcat(resnames(i),num2str(mutantnames(i)),'C') ,'Units','normalized','Linestyle','None','Fontsize',14,'FontName','arial');
    %print(strcat('/Users/Scott/Desktop/UCL Lab/Graphics/Summary2/',num2str(mutantnames(i)),'c'),'-depsc')
    %saveas(gcf,strcat('/Users/Scott/Desktop/UCL Lab/Graphics/Summary2/',num2str(mutantnames(i)),'c'),'epsc');
    %clf
end

%set(gca,'XColor',[0 0 0]);
%text(0.3,-.35,'Magnetic Field[mT]','Units','normalized','Linestyle','None','Fontsize',14,'FontName','arial');




break
plot_mutant(files3,'333','figure','4b12figure')
xlabel('Magnetic Field[mT]','Fontsize',14,'FontName','arial');


%text(0.3,-.35,'Magnetic Field[mT]','Units','normalized','Linestyle','None','Fontsize',14,'FontName','arial');

%[ax2,h2]=suplabel('d \chi" / dB','y',0.001);
%set(h2,'FontSize',18,'FontName','arial');