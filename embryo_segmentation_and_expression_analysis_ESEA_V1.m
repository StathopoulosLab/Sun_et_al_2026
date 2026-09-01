%% 11/10/2025 Udit Chakraborty ESEA_v1 %for max intensity and t stacks not for z stacks
% instructions to run: 
% hit the green run button in MATLAB and input file name in quotes
% 1. horizontal oriented embryo ('has to be horizontal')
% 2. input file name in .tif format under quotes 
% 3. image file in max intensity projection
%% 11/11/2025 Fixed issues with threshold added %AP and relative scalling
%% 11/11/2025 multiple embryo analysis and plotting overlays
%% for embryo boundary
clear
new = input('are you adding to existing data? no=0/ yes=1: ');
if new == 0
    n = 1;
else
    dat = input('input data pool file name with .mat: ');
    load(dat)
    n = size(data,2)+1;
end
while n>0
Filenm = input('choose the image '': ');
nFilenm = string(Filenm) + ".tif";
for k =1:1:4
im(:,:,k) = imread(nFilenm, k);
end
comim = sum(im, 3);  
comim = mat2gray(comim); 
imshow(comim,[])
Thr = input('choose threshold '': ');
% imshow(comim, []);
comim (comim<Thr) = 0;

[BMsk,LMsk] = bwboundaries(comim,'noholes');
BM = BMsk;
for i =1:1:size(BMsk,1)
    BM(i,2) = {i};
end
BM = BM(cellfun(@(x) numel(x) >= 100, BMsk),1:2);

statsMsk = regionprops(LMsk,'PixelIdxList');
h1 = figure;
imshow(comim,[]); % plot out all cell candidates on DIC image
set(h1,'units','normalized','outerposition',[0.5 0 0.5 1])
hold on
for b1=1:size(BM,1)
    Mskpxlist = BM{b1,1};
    plot(Mskpxlist(:,2),Mskpxlist(:,1),'b');
    cell_y = Mskpxlist(1,1); cell_x = Mskpxlist(1,2);
    metric_string = int2str(BM{b1,2});
    %metric_string = int2str(b1);
    text(cell_x-5,cell_y-5,metric_string,'Color','r','FontSize',12,'FontWeight','bold');
end
hold off;

proceed = input('continue analysing? no=0/ yes=1: ');
if proceed ==1
% choose the mask
Mskchs = input('choose the Mask index ([idx]) : ');
cltormv = ones(size(BMsk,1),1); cltormv(Mskchs) = 0;
if isempty(statsMsk(cltormv>0))
Mskidx = cell2mat(struct2cell(statsMsk(cltormv>0)));
else
Mskidx = cell2mat(struct2cell(statsMsk(cltormv>0))');
end
close all

bw = imbinarize(comim,0);
imshow(bw)
cellMsk5 = bw;
if ~isempty(Mskidx)
    cellMsk5(Mskidx)= 0;
end
cellMsk_perim = bwperim(cellMsk5);

Mskpxlist = BMsk{Mskchs,1};
cell_x = Mskpxlist(:,2);
xx = min(cell_x):max(cell_x);
x = (xx - xx(1))/(xx(end) - xx(1));
% xs =0:0.0001:max(x);
%% for each channel expression data
%%channel 1
im1 = im(:,:,1);
im1(~cellMsk5) = 0;
subplot(2,4,1)
imshow((im1),[])
J1 = medfilt2(im1);
B1 = imgaussfilt(J1, 20);
J1 = J1-B1;
Exp1 = sum(J1);
Exp1 = Exp1(min(cell_x):max(cell_x));
Exp1 = Exp1 - (mean(Exp1));  %mean subtract
Exp1(Exp1<0) = 0;
% ys1 = spline(x,Exp1,xs);


%%channel 2
im2 = im(:,:,2);
im2(~cellMsk5) = 0;
subplot(2,4,2)
imshow((im2),[])
J2 = medfilt2(im2);
B2 = imgaussfilt(J2, 50);
J2 = J2-B2;
Exp2 = sum(J2);
Exp2 = Exp2(min(cell_x):max(cell_x));
Exp2 = Exp2 - (mean(Exp2));  %mean subtract
% Exp2 = Exp2 - (std(Exp2));
Exp2(Exp2<0) = 0;
% ys2 = spline(x,Exp2,xs);


%%channel 3
im3 = im(:,:,3);
im3(~cellMsk5) = 0;
subplot(2,4,3)
imshow((im3),[])
J3 = medfilt2(im3);
B3 = imgaussfilt(J3, 50);
J3 = J3-B3;
Exp3 = sum(J3);
Exp3 = Exp3(min(cell_x):max(cell_x));
Exp3 = Exp3 - (mean(Exp3));  %mean subtract
% Exp3 = Exp3 - (std(Exp3));
Exp3(Exp3<0) = 0;
% ys3 = spline(x,Exp3,xs);


%%channel 4
im4 = im(:,:,4);
im4(~cellMsk5) = 0;
subplot(2,4,4)
imshow((im4),[])
J4 = medfilt2(im4);
B4 = imgaussfilt(J4, 50);
J4 = J4-B4;
Exp4 = sum(J4);
Exp4 = Exp4(min(cell_x):max(cell_x));
Exp4 = Exp4 - (mean(Exp4)); %mean subtract
% Exp4 = Exp4 - (std(Exp4));
Exp4(Exp4<0) = 0;
% ys4 = spline(x,Exp4,xs);

%plotting
subplot(2,4,5)
plot(x, Exp1(end,:),'cyan', 'lineWidth', 2);
hold on
plot(x, smooth(Exp1(end,:),20),'k', 'lineWidth', 1);
ylabel('expression levels')
xlabel('% of AP')
xlim([0 max(x)])
ylim([0 max([Exp1(:); Exp2(:); Exp3(:); Exp4(:)])])

subplot(2,4,6)
plot(x, Exp2(end,:),'k', 'lineWidth', 2);
hold on
plot(x, smooth(Exp2(end,:),20),'k', 'lineWidth', 1);
ylabel('expression levels')
xlabel('% of AP')
xlim([0 max(x)])
ylim([0 max([Exp1(:); Exp2(:); Exp3(:); Exp4(:)])])

subplot(2,4,7)
plot(x, Exp3(end,:),'g', 'lineWidth', 2);
hold on
plot(x, smooth(Exp3(end,:),20),'k', 'lineWidth', 1);
ylabel('expression levels')
xlabel('% of AP')
xlim([0 max(x)])
ylim([0 max([Exp1(:); Exp2(:); Exp3(:); Exp4(:)])])

subplot(2,4,8)
plot(x, Exp4(end,:),'magenta', 'lineWidth', 2);
hold on
plot(x, smooth(Exp4(end,:),20),'k', 'lineWidth', 1); 
ylabel('expression levels')
xlabel('% of AP')
xlim([0 max(x)])
ylim([0 max([Exp1(:); Exp2(:); Exp3(:); Exp4(:)])])

figure(2)
img = comim;
img(~cellMsk5) = 0;
% subplot(2,1,1)
% imshow((img),[])
% subplot(2,1,2)
hold on
plot(x, Exp4(end,:),'magenta', 'lineWidth', 4);
plot(x, smooth(Exp4(end,:),20),'k', 'lineWidth', 2);
idx4 = find(smooth(Exp4(end,:),20)==max(smooth(Exp4(end,:),20)));
maxX4 = x(idx4);
text(0.25,1000,['max is at =', (num2str(x(idx4)))])
xline (maxX4, 'lineWidth', 2)

plot(x, Exp3(end,:),'g', 'lineWidth', 4);
idx3 = find(smooth(Exp3(end,:),20)==max(smooth(Exp3(end,:),20)));
maxX3 = x(idx3);
text(0.5,2000,['max is at =', num2str(x(idx3))])
xline (maxX3, 'lineWidth', 2)
plot(x, smooth(Exp3(end,:),20),'k', 'lineWidth', 2);

plot(x, Exp2(end,:),'k', 'lineWidth', 4);
plot(x, smooth(Exp2(end,:),20),'k', 'lineWidth', 2)
ylabel('Expression levels')
xlabel('% of AP')
% ylabel('Normalized expression levels')
%plot(x, Exp4(end,:)./max([Exp1(:); Exp2(:); Exp3(:); Exp4(:)]) ,'magenta', 'lineWidth', 2);
%plot(x, Exp3(end,:)./max([Exp1(:); Exp2(:); Exp3(:); Exp4(:)]),'g', 'lineWidth', 2);
%plot(x, Exp2(end,:)./max([Exp1(:); Exp2(:); Exp3(:); Exp4(:)]),'k', 'lineWidth', 2);
%plot(x, Exp1(end,:)./max([Exp1(:); Exp2(:); Exp3(:); Exp4(:)]),'cyan', 'lineWidth', 2);
xlim([0 max(x)])
hold off

figure(3)
img = comim;
img(~cellMsk5) = 0;
hold on
plot(x, Exp4(end,:),'magenta', 'lineWidth', 4);
plot(x, smooth(Exp4(end,:),20),'k', 'lineWidth', 2);
idx4 = find(smooth(Exp4(end,:),20)==max(smooth(Exp4(end,:),20)));
maxX4 = x(idx4);
text(0.25,1000,['max is at =', (num2str(x(idx4)))])
xline (maxX4, 'lineWidth', 2)

plot(x, Exp3(end,:),'g', 'lineWidth', 4);
plot(x, smooth(Exp3(end,:),20),'k', 'lineWidth', 2);
idx3 = find(smooth(Exp3(end,:),20)==max(smooth(Exp3(end,:),20)));
maxX3 = x(idx3);
text(0.5,2000,['max is at =', num2str(x(idx3))])
xline (maxX3, 'lineWidth', 2)

plot(x, Exp1(end,:),'cyan', 'lineWidth', 4);
plot(x, smooth(Exp1(end,:),20),'k', 'lineWidth', 2)
ylabel('Expression levels')
xlabel('% of AP')
% ylabel('Normalized expression levels')
%plot(x, Exp4(end,:)./max([Exp1(:); Exp2(:); Exp3(:); Exp4(:)]) ,'magenta', 'lineWidth', 2);
%plot(x, Exp3(end,:)./max([Exp1(:); Exp2(:); Exp3(:); Exp4(:)]),'g', 'lineWidth', 2);
%plot(x, Exp1(end,:)./max([Exp1(:); Exp2(:); Exp3(:); Exp4(:)]),'cyan', 'lineWidth', 2);
xlim([0 max(x)])
hold off


data{n} = [x',Exp1',Exp2',Exp3',Exp4'];
sum_d{n} = [sum(Exp1), sum(Exp2), sum(Exp3), sum(Exp4)]
imj{n} = im;
keep = input('add last embryo to list? no=0/ yes=1: ');
if keep == 0
    data(n) = [];
    imj(n) = [];
    n = n-1;
end
n = n+1;
end
h1 = figure(1);
nch = string(Filenm)+"_1";
saveas(h1,nch, 'fig');

h2 = figure(2);
nch = string(Filenm)+"_2";
saveas(h2,nch, 'fig');

h3 = figure(3);
nch = string(Filenm)+"_3";
saveas(h3,nch, 'fig');

exit = input('want to exit experiment? no=0/ yes=1:: ');
clearvars -except data imj n exit
if exit ==1
    break;
end

end

%plotting overlay multiple replicates/ stage data
figure(4)
for i =1:1:size(data,2)
plot(data{1,i}(:,1), data{1,i}(:,2), 'cyan', 'lineWidth', 2);
hold on
title('channel 1')
ylabel('expression levels')
xlabel('% of AP')
end

figure(5)
for i =1:1:size(data,2)
plot(data{1,i}(:,1), data{1,i}(:,3), 'k', 'lineWidth', 2);
hold on
title('channel 2')
ylabel('expression levels')
xlabel('% of AP')
end

figure(6)
for i =1:1:size(data,2)
plot(data{1,i}(:,1), data{1,i}(:,4), 'g', 'lineWidth', 2);
hold on
title('channel 3')
ylabel('expression levels')
xlabel('% of AP')
end

figure(7)
for i =1:1:size(data,2)
plot(data{1,i}(:,1), data{1,i}(:,5), 'magenta', 'lineWidth', 2);
hold on
title('channel 4')
ylabel('expression levels')
xlabel('% of AP')
end