%% Introduction
% the goal for this step is to use new GCN method to compare the young data
% to the old data. we will also compare a random data to ensure our results

% first we will start with preprocessing of the data.
% the data is already normalized, we will limit it to 3000, and the
% threshold DEG is irrelevant for now.
path = '/Users/Daniel1/תואר שני/תזה/אריזה/data/';
name = '2000_Young,LTHSC,C57BL6,Non Cycling.mat';
[type,young,old] = preprocessing(name,path,3000);


% I chose the threshold to be 0.7 arbitrary.
threshold = 0.7;
[~,Ay] = GCNcutoff(young,threshold); %for young
linksy = sum(Ay,1);
indexy = find(linksy>1); % the indexes of the genes with more than 1 connection
Sy = Ay(indexy,indexy); % taking the the genes

[~,Ao] = GCNcutoff(old,threshold); %for old
linkso = sum(Ao,1);
indexo = find(linkso>1); % the indexes of the genes with more than 1 connection
So = Ao(indexo,indexo); % taking the the genes

%plotting
figure('Name',type); % open new figure named the type of data
N = size(Sy,1); 
%creating a circle with N nodes
cor = [sin((1:N)*2*pi/N);cos((1:N)*2*pi/N)]';
subplot(1,2,1); % present the graph in left side of the figure
gplot(Sy,cor); %plotting the relations on the circle
% display the number of links next to each node
text(cor(:,1)*1.05,cor(:,2)*1.05,string(sum(Sy,1)))
title('young','Units', 'normalized','Position', [0.5, 1.05,0])
axis off
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

N = size(So,1); 
%creating a circle with N nodes
cor = [sin((1:N)*2*pi/N);cos((1:N)*2*pi/N)]';
subplot(1,2,2); % present the graph in right side of the figure
gplot(So,cor); %plotting the relations on the circle
% display the number of links next to each node
text(cor(:,1)*1.05,cor(:,2)*1.05,string(sum(So,1)))
title('old','Units', 'normalized','Position', [0.5,1.05,0])
axis off
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

%% identifing strong genes
% the goal for this step is to identify the genes with strong pairwise
% correlation.

% first we will start with preprocessing of the data.
% the data is already normalized, we will limit it to 3000, and the
% threshold DEG is irrelevant for now.
path = '/Users/Daniel1/תואר שני/תזה/אריזה/data/';
name = '2000_Young,LTHSC,C57BL6,Non Cycling.mat';
[type,young,old] = preprocessing(name,path,3000);

alpha = 10^-4;

% finding the number of links for each gene.
[~,Ay] = GCN(young,alpha); %for young
linksy = sum(Ay,1);

[~,Ao] = GCN(old,alpha); %for old
linkso = sum(Ao,1);

% finding all the genes that has links and therefor significant
index = find(linkso+linksy>0); %both young and old 

figure('Name',type); % open new figure named the type of data
hold all;

scatter(linksy(index),linkso(index),'LineWidth',2) %ploting all the significant genes 
L = min([linksy,linkso]):max([linksy,linkso],'LineWidth',2); 
plot(L,L) %plotting the barier between young and old
title("# of significant links of genes")
xlabel("young")
ylabel("old")
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

%% ROC curve

n_genes = 1000;  % Number of genes in the GRN
n_cells = 100;
alpha = 10.^(-3:-0.5:-7);
% alpha = (0.1:0.1:1).*10^-4;
p = 0.5;

PD = zeros(length(alpha),1);
FAR = zeros(length(alpha),1);

t_span = [0 20];                % Time interval for ODE simulation
data = zeros(n_genes,n_cells);  % Matrix to store final simulated data
wself = rand(1,n_genes) * 2;
k = 3;  % Number of expected connections per gene

% Creating a base model
w_ = rand(n_genes) * 2;         % The weights for each gene (range: 0–2)
L1 = randi(n_genes,n_genes)<=k; % The structure connections between the genes
L = triu(L1,1)+triu(L1,1)';     % Make the matrix symmetric

ode = @(t, x, w, w2) -x + w2'.*(x ./ (1 + x)) + sum(w.*(x ./ (1 + x)),1)';

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

    [t, x] = ode45(@ (t, x) ode(t,x,w,wself), t_span, x0); % Solve ODE
    data(:,i) = x(end,:); % Save final state
end
data(data<10^-7)=0;
    
% % bootstrap characteristics
% numsub = 20; % number of times
% persub = 0.75; % precent from the data

for j = 1:length(alpha)
    
    % confo = zeros(2,numsub);
    % 
    % parfor j=1:numsub
    %     index = randperm(n_cells); % randomly picking the cells.
    %     data2 = data(:,index(1:floor(n_cells*persub))); % lessen the data.
    %     [~,A] = GCN(data2,alpha(i)); % perform the function
    %     compatibility = (A-1)*2+L;
    %     TP = sum(compatibility==1,'all');
    %     FN  = sum(compatibility==-1,'all');
    %     FP = sum(compatibility==0,'all');
    %     TN = sum(compatibility==-2,'all');
    % 
    %     confo(:,j) = [TP./(TP + FN),FP./(FP + TN)];
    % end
    % 
    % PD(i,:) = [mean(confo(1,:)) std(confo(1,:))];
    % FAR(i,:) = [mean(confo(2,:)) std(confo(2,:))];


    [~,A] = GCN(data,alpha(j));
    compatibility = (A-1)*2+L;
    TP = sum(compatibility==1,'all');
    FN  = sum(compatibility==-1,'all');
    FP = sum(compatibility==0,'all');
    TN = sum(compatibility==-2,'all');

    PD(j) = TP./(TP + FN);
    FAR(j) = FP./(FP + TN);

end

figure
plot(FAR,PD,'LineWidth',2)
% errorbar(FAR(:,1),PD(:,1),PD(:,2),PD(:,2),FAR(:,2),FAR(:,2))
title("ROC curve")
xlabel("1 - Specificity")
ylabel("Sensitivity")
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

%% data size effect

data = load('/Users/Daniel1/תואר שני/תזה/אריזה/data/2000_Young,LTHSC,C57BL6,Non Cycling.mat');
data = data.data2;
[data,~] = reduce(data,data,3000);
numcells = size(data,2);

%bootstrap characteristics
numsub = 10; % number of times
k = 10:10:numcells;
resultsgcn = zeros(length(k),2);
resultsgcncut = zeros(length(k),6);
threshold = [0.3 0.5 0.7];

alpha = 10^-4;

for i=1:length(k)
    %gcn
    persub = k(i)/numcells;
    gcn = @(data) GCN(data,alpha);
    BSgcn = Bootstrap(gcn,data,numsub,persub); %bootstrap gcn real data
    MND = BSgcn;
    resultsgcn(i,:) = [mean(MND) std(MND)];
end

% for j=1:length(threshold)
%     for i=1:length(k)
%         %gcncutoff
%         persub = k(i)/numcells;
%         gcncutoff = @(data) GCNcutoff(data,threshold(j));
%         BSgcncut = Bootstrap(gcncutoff,data,numsub,persub); %bootstrap gcn real data
%         resultsgcncut(i,(2*j-1):2*j) = [mean(BSgcncut) std(BSgcncut)];
%     end
% end

%plotting
figure('Name',"MND"); % open new figure named the type of data
errorbar(k,resultsgcn(:,1),resultsgcn(:,2),'LineWidth',2)
title("MND")
xlabel('# cells')
ylabel('MND')
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

% %plooting cutoff
% figure('Name',"MND cutoff"); % open new figure named the type of data
% hold on
% errorbar(k,resultsgcncut(:,1),resultsgcncut(:,2),'LineWidth',2)
% hold on
% errorbar(k,resultsgcncut(:,3),resultsgcncut(:,4),'LineWidth',2)
% hold on
% errorbar(k,resultsgcncut(:,5),resultsgcncut(:,6),'LineWidth',2)
% title("MND with threshold cutoff")
% xlabel('# cells')
% ylabel('MND')
% legend('0.3 cutoff','0.5 cutoff','0.7 cutoff')
% set(gcf, 'Color', 'w');
% box on;
% set(gca, 'FontSize', 16);

%% data size effect cases


young = load('/Users/Daniel1/תואר שני/תזה/אריזה/data/2000_Young,LTHSC,C57BL6,Non Cycling.mat');
young = young.data2;
[young,~] = reduce(young,young,3000);
N = size(young,1);
k = [100,60,20];
alpha = 10^-4;

fig = figure('Name',"significant correlation");
hold all;
for i=1:length(k)

    data = young(:,randperm(113,k(i)));
    [~,A,C] = GCN(data,alpha);

    [N,edges] = histcounts(C(A));
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    plot(edges, N,"LineWidth",6);
end

title('pairwise genes significant correlations')
legend("100 cells","60 cells","20 cells");
xlabel("correlation values")
ylabel('number of pairs of genes')
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

%% data size effect model 1

n_genes = 1000;  % Number of genes in the GRN
n_cells = 100;
p = 0.5;
alpha = 10^-4;
k = 3;  % Number of expected connections per gene

ode1 = @(t, x, w, w2) -x + w2'.*(x ./ (1 + x)) + sum(w.*(x ./ (1 + x)),1)';


% bootstrap characteristics
numsub = 30; % number of times
persub = 0.1:0.1:1; % precent from the data

PD = zeros(length(persub),2);
FAR = zeros(length(persub),2);
results = zeros(length(persub),2);

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


for i=1:length(persub)
    
    %GCN
    % gcn = @(data) GCN(data,alpha);
    % BSgcn = Bootstrap(gcn,data,numsub,persub(i));
    
    BSgcn = zeros(1,numsub); % the vector of results
    confo = zeros(2,numsub);

    parfor j=1:numsub
        index = randperm(n_cells); % randomly picking the cells.
        data2 = data(:,index(1:floor(n_cells*persub(i)))); % lessen the data.
        [BSgcn(j),A] = GCN(data2,alpha); % perform the function
        compatibility = (A-1)*2+L;
        TP = sum(compatibility==1,'all');
        FN  = sum(compatibility==-1,'all');
        FP = sum(compatibility==0,'all');
        
        confo(:,j) = [TP./(TP + FN),TP./(TP + FP)];
    end

    results(i,:) = [mean(BSgcn) std(BSgcn)] ;
    
    % index = randperm(n_cells); % randomly picking the cells.
    % data2 = data(:,index(1:floor(persub(i)*n_cells))); % lessen the data.
    % [~,A] = GCN(data2,alpha);
    % compatibility = (A-1)*2+L;
    % TP = sum(compatibility==1,'all');
    % FN  = sum(compatibility==-1,'all');
    % FP = sum(compatibility==0,'all');
    % TN = sum(compatibility==-2,'all');
        
    PD(i,:) = [mean(confo(1,:)) std(confo(1,:))];
    FAR(i,:) = [mean(confo(2,:)) std(confo(2,:))];


end

%plotting
figure("Name",n_cells+"X"+n_genes); % open new figure named the type of data
hold all;

%plotting the gcn
subplot(2,1,1);
errorbar(floor(persub.*n_cells),results(:,1),results(:,2),'LineWidth',2)
title("MND")
xlabel("# cells")
ylabel('MND')
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

%plooting the recall
subplot(2,2,3);
errorbar(floor(persub.*n_cells),PD(:,1),PD(:,2),'LineWidth',2)
title("Recall")
xlabel("# cells")
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

%plotting the precision
subplot(2,2,4);
errorbar(floor(persub.*n_cells),FAR(:,1),FAR(:,2),'LineWidth',2)
title("Precision")
xlabel("# cells")
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);


