%% fisher's test - model 1

n_genes = 100;  % Number of genes in the GRN
n_cells = 50;
p = 0.5;
alpha = 10.^-4;
k = 3;  % Number of expected connections per gene

ode1 = @(t, x, w, w2) -x + w2'.*(x ./ (1 + x)) + sum(w.*(x ./ (1 + x)),1)';


t_span = [0 20];                % Time interval for ODE simulation
wself = rand(1,n_genes) * 2;

% Creating a base model
w_ = rand(n_genes) * 2;         % The weights for each gene (range: 0–2)
L1 = randi(n_genes,n_genes)<=k; % The structure connections between the genes
L = triu(L1,1)+triu(L1,1)';     % Make the matrix symmetric

data = zeros(n_genes,n_cells);  % Matrix to store final simulated data

    % creating the data
for i = 1:n_cells
    L2 = rand(n_genes)<=p;          % Random mask
    L2 = triu(L2,1)+triu(L2,1)';    % Symmetrize
    L2 = logical(L2.*L);            % Apply the mask on the base structure
    
    w = w_.*L;                      % Apply base weights to the base structure
    w(L2) = rand(sum(L2,"all"),1) * 2; % Randomize weights according to mask
    
    x0 = rand(1, n_genes);
    offgenes = randi(n_genes,1,5);
    x0(offgenes) = 0;
    w(offgenes,:) = 0;
    w(:,offgenes) = 0;
    
    [t, x] = ode45(@ (t, x) ode1(t,x,w,wself), t_span, x0); % Solve ODE
    data(:,i) = x(end,:); % Save final state
end
data(data<10^-7)=0;

[~,A] = GCN(data,alpha);
compatibility = (A-1)*2+L;
TP = sum(compatibility==1,'all');
FN  = sum(compatibility==-1,'all');
FP = sum(compatibility==0,'all');
TN = sum(compatibility==-2,'all');

[~,h] = fishertest([TP FN ; FP TN]);


%plotting
figure("Name",n_cells+"X"+n_genes); % open new figure named the type of data

imagesc(compatibility)
title("fishertest = "+h)
colorbar('Ticks',[-2,-1,0,1],'TickLabels',{'TN','FN','FP','TP'})
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);



%% idetify clusters

n_genes = 100;  % Number of genes in the GRN
n_cells = 50;
p = 0.5;
alpha = 10.^-4;
k = 3;  % Number of expected connections per gene

ode1 = @(t, x, w, w2) -x + w2'.*(x ./ (1 + x)) + sum(w.*(x ./ (1 + x)),1)';


t_span = [0 20];                % Time interval for ODE simulation
wself = rand(1,n_genes) * 2;

% Creating a base model
w_ = rand(n_genes) * 2;         % The weights for each gene (range: 0–2)
L1up = randi(n_genes,n_genes/2)<=k; % The structure connections between the genes
L1down = randi(n_genes,n_genes/2)<=k; % The structure connections between the genes
L1 = [L1up zeros(n_genes/2); zeros(n_genes/2) L1down];
L = triu(L1,1)+triu(L1,1)';     % Make the matrix symmetric

data = zeros(n_genes,n_cells);  % Matrix to store final simulated data

    % creating the data
for i = 1:n_cells
    L2 = rand(n_genes)<=p;          % Random mask
    L2 = triu(L2,1)+triu(L2,1)';    % Symmetrize
    L2 = logical(L2.*L);            % Apply the mask on the base structure
    
    w = w_.*L;                      % Apply base weights to the base structure
    w(L2) = rand(sum(L2,"all"),1) * 2; % Randomize weights according to mask
    
    x0 = rand(1, n_genes);
    offgenes = randi(n_genes,1,5);
    x0(offgenes) = 0;
    w(offgenes,:) = 0;
    w(:,offgenes) = 0;
    
    [t, x] = ode45(@ (t, x) ode1(t,x,w,wself), t_span, x0); % Solve ODE
    data(:,i) = x(end,:); % Save final state
end
data(data<10^-7)=0;

[~,A] = GCN(data,alpha);
compatibility = (A-1)*2+L;

%plotting
figure("Name",n_cells+"X"+n_genes); % open new figure named the type of data

imagesc(compatibility)
title("Identify clusters")
colorbar('Ticks',[-2,-1,0,1],'TickLabels',{'TN','FN','FP','TP'})
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);


%% Robustness fig 1

data = load('/Users/Daniel1/תואר שני/תזה/אריזה/data/2000_Young,LTHSC,C57BL6,Non Cycling.mat');
data = data.data2;
[data,~] = reduce(data,data,3000);
numcells = size(data,2);

%bootstrap characteristics
numsub = 30; % number of times
k = 10:10:numcells;
resultsgcn = zeros(length(k),2);
resultsgcl = zeros(length(k),2);

alpha = 10^-4;
gcn = @(data) GCN(data,alpha);

for i=1:length(k)
   
    persub = k(i)/numcells;  
    %gcn
    BSgcn = Bootstrap(gcn,data,numsub,persub); %bootstrap gcn real data
    resultsgcn(i,:) = [mean(BSgcn) std(BSgcn)];
    %gcl
    BSgcl = Bootstrap(@GCL,data,numsub,persub); %bootstrap gcn real data
    resultsgcl(i,:) = [mean(BSgcl) std(BSgcl)];
end

resultsgcncut = zeros(length(k),6);
threshold = [0.3 0.5 0.7];

for j=1:length(threshold)
    for i=1:length(k)
        %gcncutoff
        persub = k(i)/numcells;
        gcncutoff = @(data) GCNcutoff(data,threshold(j));
        BSgcncut = Bootstrap(gcncutoff,data,numsub,persub); %bootstrap gcn real data
        resultsgcncut(i,(2*j-1):2*j) = [mean(BSgcncut) std(BSgcncut)];
    end
end
%%
%plotting
figure('Name',"Robutness"); % open new figure named the type of data
subplot(1,3,1)
errorbar(k,resultsgcl(:,1),resultsgcl(:,2),'LineWidth',2)
title("GCL")
xlabel('# cells')
ylabel('GCL')
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

subplot(1,3,2)
errorbar(k,resultsgcn(:,1),resultsgcn(:,2),'LineWidth',2)
title("MND")
xlabel('# cells')
ylabel('MND')
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

subplot(1,3,3)
%plooting cutoff
hold on
errorbar(k,resultsgcncut(:,1),resultsgcncut(:,2),'LineWidth',2)
hold on
errorbar(k,resultsgcncut(:,3),resultsgcncut(:,4),'LineWidth',2)
hold on
errorbar(k,resultsgcncut(:,5),resultsgcncut(:,6),'LineWidth',2)
title("MND with threshold cutoff")
xlabel('# cells')
ylabel('MND')
legend('0.3 cutoff','0.5 cutoff','0.7 cutoff')
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

%% Consistency
%the goal for this step is to compare the GCN method to the GCL method.
% we will also compare a random data to normal the comparison.

files = dir('data/*Young*.mat');

alpha = 10.^-4;
gcn = @(data) GCN(data,alpha);

for i=1:length(files)
% first we will start with preprocessing of the data.
% the data is already normalized, we will limit it to 3000, and the
% threshold DEG is irrelevant for now.
path = [files(i).folder '/'];
[type,young,old] = preprocessing(files(i).name,path,3000);

%bootstrap characteristics
numsub = 30; % number of times
persub = 0.75; % precent from the data
BSgcl = zeros(2,numsub);
BSgcn = zeros(2,numsub);

%GCL
BSgcl(1,:) = Bootstrap(@GCL,young,numsub,persub); %bootstrap gcl young data
BSgcl(2,:) = Bootstrap(@GCL,old,numsub,persub); %bootstrap gcl old data
% overlap = [max(min(BSgcl,[],2)), min(max(BSgcl,[],2))];
% sum(BSgcl >= overlap(1) & BSgcl <= overlap(2),"all")/(2*numsub);
Cdgcl = -diff(mean(BSgcl,2))/sqrt(sum(std(BSgcl,[],2).^2)/2);

%GCN
BSgcn(1,:) = Bootstrap(gcn,young,numsub,persub); %bootstrap gcn young data
BSgcn(2,:) = Bootstrap(gcn,old,numsub,persub); %bootstrap gcn old data
% overlap = [max(min(BSgcn,[],2)), min(max(BSgcn,[],2))];
% sum(BSgcn >= overlap(1) & BSgcn <= overlap(2),"all")/(2*numsub);
Cdgcn = -diff(mean(BSgcn,2))/sqrt(sum(std(BSgcn,[],2).^2)/2);

%plotting
fig = figure('Name',type,'Color','w','Position', [796,360,240,420]); % open new figure named the type of data
t = tiledlayout(1, 2);  % (rows, columns)

% Optional: control spacing
t.TileSpacing = 'compact';
t.Padding = 'compact';

% hold all;
% Define the colors
custom_colors = {[0 0 1],[1 1 0]};  % white, yellow, blue

%plooting the gcl
nexttile(1);
boxplot(BSgcl',["young","old"],'Widths',1)
title("GCL")

h = findobj(gca, 'Tag', 'Box');
for j = 1:2
    patch(get(h(j), 'XData'), get(h(j), 'YData'), custom_colors{j}, 'FaceAlpha', 0.5);
end

text(0.98, 0.95, num2str(Cdgcl,3),'Units', 'normalized', 'HorizontalAlignment','right');

%plotting the gcn
nexttile(2);
boxplot(BSgcn',["young","old"],'Widths',1)
title("MND")

h = findobj(gca, 'Tag', 'Box');
for j = 1:2
    patch(get(h(j), 'XData'), get(h(j), 'YData'), custom_colors{j}, 'FaceAlpha', 0.5);
end

text(0.98, 0.95, num2str(Cdgcn,3),'Units', 'normalized','HorizontalAlignment','right');

% saveas(fig,['compare2/' type '.fig']);
% close(fig);

end


%% sensitivity to noise - comparison with replacing
% the goal for this step is to see how the significant genes affect the 
% results. we will construct a subset that excludes these genes and compare
% the two approaches on that subset

path = '/Users/Daniel1/תואר שני/תזה/אריזה/data/';
name = '2000_Young,LTHSC,C57BL6,Non Cycling.mat';
% name = '2000_Young,STHSC,DBA,Non Cycling.mat';
% name = '2000_Young,MPP,CellRep2.mat';
[type,young,old] = preprocessing(name,path,5000);

N = size(young,1);
alpha = 10^-4;
gcn = @(data) GCN(data,alpha);

for k = [1000,500,100]

    [young,old] = reduce(young,old,3000+k);
    N = size(young,1);

    % finding the number of links for each gene.
    [~,A] = GCN(young,alpha); %for young
    slinksy = sum(A,1);
    [~,A] = GCN(old,alpha); %for old
    sllinkso = sum(A,1);
    % finding all the genes that has links and therefor significant
    % sorting the genes for the most present ones.
    [~,Iy] = sort(slinksy,'descend');
    [~,Io] = sort(sllinkso,'descend');
    
    WOTy = Iy(k+1:N);
    WOTo = Io(k+1:N);
    WOTindex = intersect(WOTy,WOTo);

    young_WOT = young(WOTindex,:);
    old_WOT = old(WOTindex,:);

    WOBy = Iy(1:N-k);
    WOBo = Io(1:N-k);
    WOBindex = intersect(WOBy,WOBo);
   
    young_WOB = young(WOBindex,:);
    old_WOB = old(WOBindex,:);
    
    %bootstrap characteristics
    numsub = 30; % number of times
    persub = 0.75; % precent from the data
    BSgclWOB = zeros(2,numsub);
    BSgcnWOB = zeros(2,numsub);
    BSgclWOT = zeros(2,numsub);
    BSgcnWOT = zeros(2,numsub);
    
    %GCL
    BSgclWOB(1,:) = Bootstrap(@GCL,young_WOB,numsub,persub); %bootstrap gcl young data
    BSgclWOB(2,:) = Bootstrap(@GCL,old_WOB,numsub,persub); %bootstrap gcl old data
    BSgclWOT(1,:) = Bootstrap(@GCL,young_WOT,numsub,persub); %bootstrap gcl young data
    BSgclWOT(2,:) = Bootstrap(@GCL,old_WOT,numsub,persub); %bootstrap gcl old data

    CdgclWOB = -diff(mean(BSgclWOB,2))/sqrt(sum(std(BSgclWOB,[],2).^2)/2);
    CdgclWOT = -diff(mean(BSgclWOT,2))/sqrt(sum(std(BSgclWOT,[],2).^2)/2);

    %GCN
    BSgcnWOB(1,:) = Bootstrap(gcn,young_WOB,numsub,persub); %bootstrap gcn young data
    BSgcnWOB(2,:) = Bootstrap(gcn,old_WOB,numsub,persub); %bootstrap gcn old data
    BSgcnWOT(1,:) = Bootstrap(gcn,young_WOT,numsub,persub); %bootstrap gcn young data
    BSgcnWOT(2,:) = Bootstrap(gcn,old_WOT,numsub,persub); %bootstrap gcn old data
    
    CdgcnWOB = -diff(mean(BSgcnWOB,2))/sqrt(sum(std(BSgcnWOB,[],2).^2)/2);
    CdgcnWOT = -diff(mean(BSgcnWOT,2))/sqrt(sum(std(BSgcnWOT,[],2).^2)/2);

    %plotting
    figure('Name',[type,' #cells: ',num2str(k)],'Color','w'); % open new figure named the type of data
    tiledlayout(2, 1);  % (rows, columns) 
    
    % Define the colors
    custom_colors = {[0 0 1],[1 1 0],[0 0 1],[1 1 0]};  % yellow, blue
    
    %plooting the gcl
    nexttile;
    boxplot([BSgclWOB;BSgclWOT]',["young WOB","old WOB","young WOT","old WOT"],'Widths',1)
    title("GCL")
    text(0.48, 0.95, num2str(CdgclWOB),'Units', 'normalized', ...
    'HorizontalAlignment','right');
    text(0.98, 0.95, num2str(CdgclWOT),'Units', 'normalized', ...
    'HorizontalAlignment','right');
    
    h = findobj(gca, 'Tag', 'Box');
    for j = 1:4
        patch(get(h(j), 'XData'), get(h(j), 'YData'), custom_colors{j}, 'FaceAlpha', 0.5);
    end
    
    nexttile;
    boxplot([BSgcnWOB;BSgcnWOT]',["young WOB","old WOB","young WOT","old WOT"],'Widths',1)
    title("MND")
    text(0.48, 0.95, num2str(CdgcnWOB),'Units', 'normalized', ...
    'HorizontalAlignment','right');

    text(0.98, 0.95, num2str(CdgcnWOT),'Units', 'normalized', ...
    'HorizontalAlignment','right');

    h = findobj(gca, 'Tag', 'Box');
    for j = 1:4
        patch(get(h(j), 'XData'), get(h(j), 'YData'), custom_colors{j}, 'FaceAlpha', 0.5);
    end

    % saveas(fig,['compare/Sensitivity/' type '.fig']);
    % close(fig);
end

%% p change cooperation

n_genes = 1000;  % Number of genes in the GRN
n_cells = 100;
p = 0.1:0.1:1;
alpha = 10^-4;
k = 3;  % Number of expected connections per gene

ode1 = @(t, x, w, w2) -x + w2'.*(x ./ (1 + x)) + sum(w.*(x ./ (1 + x)),1)';
ode3 = @(t, x, w3, w2) -x + w2'.*(x ./ (1 + x)) + coop(x,w3);


t_span = [0 20];                % Time interval for ODE simulation
wself = rand(1,n_genes) * 2;

% Creating a base model
w_ = rand(n_genes) * 2;         % The weights for each gene (range: 0–2)
L1 = randi(n_genes,n_genes)<=k; % The structure connections between the genes
L = triu(L1,1)+triu(L1,1)';     % Make the matrix symmetric

% Creating a base model
w3_ = rand(n_genes,n_genes,n_genes) * 2;            % The weights for each gene (range: 0–2)
L3 = randi(n_genes^2,n_genes,n_genes,n_genes)<=k;    % The structure connections between the genes


% bootstrap characteristics
numsub = 30; % number of times
persub = 0.75; % precent from the data
model1 = zeros(4,length(p));
model2 = zeros(4,length(p));
gcn = @(data) GCN(data,alpha);

for j = 1:length(p)
    % creating the data
    data = zeros(n_genes,n_cells);  % Matrix to store final simulated data
    coopdata = zeros(n_genes,n_cells);  % Matrix to store final simulated data

    for i = 1:n_cells
        % model 1
        L2 = rand(n_genes)<=p(j);          % Random mask
        L2 = triu(L2,1)+triu(L2,1)';    % Symmetrize
        L2 = logical(L2.*L);            % Apply the mask on the base structure
        
        w = w_.*L;                      % Apply base weights to the base structure
        w(L2) = rand(sum(L2,"all"),1) * 2; % Randomize weights according to mask
        
        x0 = rand(1, n_genes);
        offgenes = randi(n_genes,1,5);
        x0(offgenes) = 0;
        w(offgenes,:) = 0;
        w(:,offgenes) = 0;
    
        [t, x] = ode45(@ (t, x) ode1(t,x,w,wself), t_span, x0); % Solve ODE
        data(:,i) = x(end,:); % Save final state
    
        % model 2
        LL = rand(n_genes,n_genes,n_genes)<=p(j);          % Random mask
        LL = logical(LL.*L3);           % Apply the mask on the base structure
        w3 = w3_.*L3;                   % Apply base weights to the base structure
        w3(LL) = rand(sum(LL,"all"),1) * 2; % Randomize weights according to mask
        
                
        w3(offgenes,:,:) = 0;
        w3(:,offgenes,:) = 0;
        w3(:,:,offgenes) = 0;
    
        [t, x] = ode45(@ (t, x) ode3(t,x,w3,wself), t_span, x0); % Solve ODE
        coopdata(:,i) = x(end,:); % Save final state
    end
    data(data<10^-7)=0;
    coopdata(coopdata<10^-7)=0;
    
    BSgcl = Bootstrap(@GCL,data,numsub,persub); %bootstrap gcl real data
    BSgcn = Bootstrap(gcn,data,numsub,persub); %bootstrap gcn real data
    model1(:,j) = [mean(BSgcl),std(BSgcl),mean(BSgcn),std(BSgcn)];
    
    BSgcl = Bootstrap(@GCL,coopdata,numsub,persub); %bootstrap gcl real data
    BSgcn = Bootstrap(gcn,coopdata,numsub,persub); %bootstrap gcn real data
    model2(:,j) = [mean(BSgcl),std(BSgcl),mean(BSgcn),std(BSgcn)];
end

%plotting
figure("Name",n_cells+"X"+n_genes); % open new figure named the type of data
hold all;
%plooting the gcl
subplot(2,1,1);
errorbar(p,model1(1,:),model1(2,:),'LineWidth',2)
hold on
errorbar(p,model2(1,:),model2(2,:),'LineWidth',2)
title("GCL")
legend("activation","cooperation")
xlabel('p')
ylabel('GCL ')
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

%plotting the gcn
subplot(2,1,2);
errorbar(p,model1(3,:),model1(4,:),'LineWidth',2)
hold on
errorbar(p,model2(3,:),model2(4,:),'LineWidth',2)
title("MND")
legend("activation","cooperation")
xlabel('p')
ylabel('MND')
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

% %plooting the gcl
% subplot(2,2,3);
% errorbar(p,model2(1,:),model2(2,:))
% title("GCL cooperation")
% xlabel('p')
% ylabel('GCL')
% set(gcf, 'Color', 'w');
% box on;
% set(gca, 'FontSize', 16);
% 
% %plotting the gcn
% subplot(2,2,4);
% errorbar(p,model2(3,:),model2(4,:))
% title("GCN cooperation")
% xlabel('p')
% ylabel('GCN')
% set(gcf, 'Color', 'w');
% box on;
% set(gca, 'FontSize', 16);

%% data size effect cooperation

n_genes = 1000;  % Number of genes in the GRN
n_cells = 100;
p = 0.5;
alpha = 10^-4;
k = 3;  % Number of expected connections per gene

ode1 = @(t, x, w, w2) -x + w2'.*(x ./ (1 + x)) + sum(w.*(x ./ (1 + x)),1)';
ode3 = @(t, x, w3, w2) -x + w2'.*(x ./ (1 + x)) + coop(x,w3);


t_span = [0 20];                % Time interval for ODE simulation
wself = rand(1,n_genes) * 2;

% Creating a base model
w_ = rand(n_genes) * 2;         % The weights for each gene (range: 0–2)
L1 = randi(n_genes,n_genes)<=k; % The structure connections between the genes
L = triu(L1,1)+triu(L1,1)';     % Make the matrix symmetric

% Creating a base model
w3_ = rand(n_genes,n_genes,n_genes) * 2;            % The weights for each gene (range: 0–2)
L3 = randi(n_genes^2,n_genes,n_genes,n_genes)<=k;    % The structure connections between the genes


data = zeros(n_genes,n_cells);  % Matrix to store final simulated data
coopdata = zeros(n_genes,n_cells);  % Matrix to store final simulated data

    % creating the data
for i = 1:n_cells
    % model 1
    L2 = rand(n_genes)<=p;          % Random mask
    L2 = triu(L2,1)+triu(L2,1)';    % Symmetrize
    L2 = logical(L2.*L);            % Apply the mask on the base structure
    
    w = w_.*L;                      % Apply base weights to the base structure
    w(L2) = rand(sum(L2,"all"),1) * 2; % Randomize weights according to mask
    
    x0 = rand(1, n_genes);
    offgenes = randi(n_genes,1,5);
    x0(offgenes) = 0;
    w(offgenes,:) = 0;
    w(:,offgenes) = 0;

    [t, x] = ode45(@ (t, x) ode1(t,x,w,wself), t_span, x0); % Solve ODE
    data(:,i) = x(end,:); % Save final state

    % model 2
    LL = rand(n_genes,n_genes,n_genes)<=p;          % Random mask
    LL = logical(LL.*L3);           % Apply the mask on the base structure
    w3 = w3_.*L3;                   % Apply base weights to the base structure
    w3(LL) = rand(sum(LL,"all"),1) * 2; % Randomize weights according to mask
    
            
    w3(offgenes,:,:) = 0;
    w3(:,offgenes,:) = 0;
    w3(:,:,offgenes) = 0;

    [t, x] = ode45(@ (t, x) ode3(t,x,w3,wself), t_span, x0); % Solve ODE
    coopdata(:,i) = x(end,:); % Save final state
end
data(data<10^-7)=0;
coopdata(coopdata<10^-7)=0;

%bootstrap characteristics
numsub = 30; % number of times
persub = 0.1:0.1:1; % precent from the data
model1 = zeros(4,length(persub));
model2 = zeros(4,length(persub));
gcn = @(data) GCN(data,alpha);

for i=1:length(persub)
    BSgcl = Bootstrap(@GCL,data,numsub,persub(i)); %bootstrap gcl real data
    BSgcn = Bootstrap(gcn,data,numsub,persub(i)); %bootstrap gcn real data
    model1(:,i) = [mean(BSgcl),std(BSgcl),mean(BSgcn),std(BSgcn)];
    
    BSgcl = Bootstrap(@GCL,coopdata,numsub,persub(i)); %bootstrap gcl real data
    BSgcn = Bootstrap(gcn,coopdata,numsub,persub(i)); %bootstrap gcn real data
    model2(:,i) = [mean(BSgcl),std(BSgcl),mean(BSgcn),std(BSgcn)];

end

%plotting
figure("Name",n_cells+"X"+n_genes); % open new figure named the type of data
hold all;
%plooting the gcl
subplot(2,1,1);
errorbar(floor(persub.*n_cells),model1(1,:),model1(2,:),'LineWidth',2)
hold on
errorbar(floor(persub.*n_cells),model2(1,:),model2(2,:),'LineWidth',2)
title("GCL")
legend("activation","cooperation")
xlabel('# of cells')
ylabel('GCL ')
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

%plotting the gcn
subplot(2,1,2);
errorbar(floor(persub.*n_cells),model1(3,:),model1(4,:),'LineWidth',2)
hold on
errorbar(floor(persub.*n_cells),model2(3,:),model2(4,:),'LineWidth',2)
title("MND")
legend("activation","cooperation")
xlabel('# of cells')
ylabel('MND')
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

% %plooting the gcl
% subplot(2,2,3);
% errorbar(floor(persub.*n_cells),model2(1,:),model2(2,:))
% title("GCL cooperation")
% xlabel('# of cells')
% ylabel('GCL')
% set(gcf, 'Color', 'w');
% box on;
% set(gca, 'FontSize', 16);
% 
% %plotting the gcn
% subplot(2,2,4);
% errorbar(floor(persub.*n_cells),model2(3,:),model2(4,:))
% title("GCN cooperation")
% xlabel('# of cells')
% ylabel('GCN')
% set(gcf, 'Color', 'w');
% box on;
% set(gca, 'FontSize', 16);



%% plotting
    fig = figure('Name',type,'Color','w','Position', [796,360,240,420]); % open new figure named the type of data
    t = tiledlayout(1, 2);  % (rows, columns)
    
    % Optional: control spacing
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    
    % hold all;
    % Define the colors
    custom_colors = {[0 0 1],[1 1 0]};  % white, yellow, blue
    
    %plooting the gcl
    nexttile(1);
    boxplot(BSgcl',["young","old"],'Widths',1)
    title("GCL")
    
    h = findobj(gca, 'Tag', 'Box');
    for j = 1:2
        patch(get(h(j), 'XData'), get(h(j), 'YData'), custom_colors{j}, 'FaceAlpha', 0.5);
    end
    
    %plotting the gcn
    nexttile(2);
    boxplot(BSgcn',["young","old"],'Widths',1)
    title("MND")
    
    h = findobj(gca, 'Tag', 'Box');
    for j = 1:2
        patch(get(h(j), 'XData'), get(h(j), 'YData'), custom_colors{j}, 'FaceAlpha', 0.5);
    end
    
   

