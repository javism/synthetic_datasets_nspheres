%SYNTHETIC_DATA_NSPHERES Generate a synthetic dataset based on n-spheres
%
% Examples:
%   synthetic_datasets_nspheres('synthetic_data_nspheres_conf')
% 
% This function generates a set of synthetic dataset based on a set of 
% isotropic Gaussian distributions. The function needs a configuration file
% containing a set of variables describing synthetic data features and
% data plot options. Please, check parameters comments at example file 
% 'synthetic_data_nspheres_conf.m'
%
% The function is a prototype and does not perform error checking regarding
% configuration parameters. Please, send errors, suggestions and improvements 
% to jsanchezm at uco dot es
%
% The data generation process is described in the following conference
% article (if you use this software, please cite the article as):
%
% J. Sánchez-Monedero, P.A. Gutiérrez, M. Pérez-Ortiz, and C. Hervás-Martínez
% An n-spheres based synthetic data generator for supervised classification
% International Work-Conference on Artificial Neural Networks (IWANN 2013)
% 
% Abstract: 
% 
% Synthetic datasets can be useful in a variety of situations, specifically 
% when new machine learning models and training algorithms are developed or 
% when trying to seek the weaknesses of an specific method. In contrast to 
% real-world data, synthetic datasets provide a controlled environment for 
% analysing concrete critic points such as outlier tolerance, data dimensionality 
% influence and class imbalance, among others. In this paper, a framework for 
% synthetic data generation is developed with special attention to pattern 
% order in the space, data dimensionality, class overlapping and data multimodality. 
% Variables such as position, width and overlapping of data distributions in 
% the n-dimensional space are controlled by considering them as $n$-spheres. 
% The method is tested in the context of ordinal regression, a paradigm of 
% classification where there is an order arrangement between categories. 
% The contribution of the paper is the full control over data topology and 
% over a set of relevant statistical properties of the data. 
% 
% Bibtex: 
% 
% @INPROCEEDINGS{Sanchez-Monedero2013iwann,
%  author = {J. S\'anchez-Monedero, P.A. Guti\'errez, M. P\'erez-Ortiz, and C.	Herv\'as-Mart\'inez},
%  title = {An n-spheres based synthetic data generator for supervised classification},
%  booktitle = {International Work-Conference on Artificial Neural Networks (IWANN)},
%  year = {2013},
%  volume = {Accepted},}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) Javier Sánchez Monedero (jsanchezm at uco dot es)
% 
% This code implements the synthetic data generator described in the following publication:
%
% J. Sánchez-Monedero, P.A. Gutiérrez, M. Pérez-Ortiz, and C. Hervás-Martínez
% An n-spheres based synthetic data generator for supervised classification
% International Work-Conference on Artificial Neural Networks (IWANN 2013)
% 
% The code has been tested with Ubuntu 12.04 x86_64 / Debian Stable 6.0 and Matlab R2009a
% 
% If you use this code, please cite the associated paper
% Code updates and citing information:
% http://www.uco.es/grupos/ayrna/iwann2013-syntheticdatagenerator
% 
% AYRNA Research group's website:
% http://www.uco.es/ayrna 
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
% Licence available at: http://www.gnu.org/licenses/gpl-3.0.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function synthetic_datasets_nspheres(conf_file)

    if nargin < 1
        disp('A configuration file should be provided without ''.m'' extension')
        disp(', using default file name ''synthetic_data_nspheres_conf''')
        conf_file = 'synthetic_data_nspheres_conf';
    end
    
    % Add external dependencies
    addpath(['external_tools' filesep 'export_fig'])
    addpath(['external_tools' filesep 'plot2svg'])
    addpath('external_tools')

    % Load setup (just loading variables)
    eval(conf_file)

    for Kk = 1:size(Kset,2)
        for Ss = 1:size(Sset,2)
            K = Kset(1,Kk);
            sigma = Sset(1,Ss);
            synthetic_datasets_nspheres_priv(conf_file, muq, sigma, K,alphaq,lambdaq)
        end
    end

    % If you do not want to close all figures at once comment the following
    % line
    close all 
end 

function synthetic_datasets_nspheres_priv(config,muq,sigma,K,alphaq,lambdaq)

eval(config)

% Create necessary dirs
if ~exist(outputArffDir,'dir')
    mkdir(outputArffDir)
end

if ~exist(outputDatDir,'dir')
    mkdir(outputDatDir)
end

if ~exist(outputPicsDir,'dir')
    mkdir(outputPicsDir)
end

% Get number of classes 
J = size(N,1);


%%%%%%%%%% Generate data %%%%%%%%%% 
% Adjust the number of patterns per class with the number of modes
NP=N.*nOfModes';

% Number of patterns for each Gaussian
NN = repmat(N(1,1),nOfModes(1,1),1);
for j=2:J
    NN = [NN;repmat(N(j,1),nOfModes(1,j),1)];
end

nNM = [0; nOfModes'];
nNP2 = [0; NP];

% Adjust centers separation with the number of dimensions
alphaq = alphaq/sqrt(K);
lambdaq = lambdaq/sqrt(K);

% Setup mus
mus = zeros(sum(nOfModes),K);
for j = 1:J
    m = sum(nNM(1:j))+1;
    muq = muq + alphaq;
    mus(m,:) = muq;
    
    for m = sum(nNM(1:j))+2:(sum(nNM(1:j+1)))
        hypersphere_points=randn(K,hyper_sphere_surface_points);
        hypersphere_points2=lambdaq*normr(hypersphere_points');
        random_mu = hypersphere_points2(randi(hyper_sphere_surface_points,1),:);
        mus(m,:) = muq + random_mu;
    end
    
end

% Setup sigma/width
sigmas = repmat(sigma, J, K);
RN = zeros(1,K);

% Random sampling from multivariate Normal distribution
for j=1:(sum(nOfModes))
    RNtmp = mvnrnd(repmat(mus(j,: ),NN(j,1),1), eye(K)*sigma^2);
    RN = [RN;RNtmp];
end
clear RNtmp
RN(1,:) = [];

% Add targets

TN = ones(NP(1,1),1)*1;

for j=2:J
    TN = [TN; ones(NP(j,1),1)*j];
end

DataSet.patterns = RN; 
DataSet.targets = TN;

clear T TN

if plotPoints || plotSD
     figure; hold on;
     axis equal;
end


%%%%%%%%%% Plot data %%%%%%%%%% 

% NOTE: Some loops code are intentionally separated to rightly paint the plot

if K == 1
    for j = 1:J
        X = RN(sum(nNP2(1:j))+1:sum(nNP2(1:j+1)),1);

        if plotPoints 
            legs(j) = plot(X, zeros(size(X)),...
                 markersymbols{j}, 'MarkerSize',markersizes{j},'Color',markercolours{j});
        end
        
        if plotSD
            modas = mus(sum(nNM(1:j))+1:sum(nNM(1:j+1)),:);
            
            for m = 1:size(modas,1)
                % Mus
                legs(J+1) = plot(modas(m,1),0,'.','Color',mucolor,'MarkerSize',15);
                
                % Gaussian p.d.f.
                XG = min(RN):0.001:max(RN);
                YG = gaussmf(XG,[sigmas(m,1) modas(m,1)]);
                plot(XG,YG,'Color',markercolours{j})
            end     
        end     
    end
    
    % Build legend string
    legendStr = cell(J+1,1);
    for j=1:J
        legendStr{j,1} = sprintf('Class %d points',j); 
    end
    
    legendStr{J+1,1} = sprintf('\\mu_q');
    legend(legs,legendStr,'Location','NorthEastOutside')
    
elseif (K == 2)
    
    % Save handles for legends    
    for j = 1:J
        X = RN(sum(nNP2(1:j))+1:sum(nNP2(1:j+1)),1);
        Y = RN(sum(nNP2(1:j))+1:sum(nNP2(1:j+1)),2);       

        if plotPoints 
            legs(j*2-1) = plot(X, Y,...
                markersymbols{j}, 'MarkerSize',markersizes{j},'Color',markercolours{j});
        end

        if plotSD
            modas = mus(sum(nNM(1:j))+1:sum(nNM(1:j+1)),:);

            for m = 1:size(modas,1)
                legs(j*2) = plotEllipse(modas(m,:), sigmas(j,:)*numberSD,markercolours{j});
            end
        end
    end
    
    for j = 1:J
        modas = mus(sum(nNM(1:j))+1:sum(nNM(1:j+1)),:);

        for m = 1:size(modas,1)
            legs(J*2+1) = plot(modas(m,1),modas(m,2),'.','Color',mucolor,'MarkerSize',15);
        end
    end
    
    % Build legend string
    legendStr = cell(2*J+1,1);
    for j=1:J
        legendStr{2*j-1,1} = sprintf('Class %d points',j); 
        legendStr{2*j,1} = sprintf('Class %d %d\\sigma',j,numberSD);
    end
    
    legendStr{2*J+1,1} = sprintf('\\mu_q');
    legend(legs,legendStr,'Location','NorthEastOutside')
    
elseif K >= 3
    
    for j = 1:J
        if plotPoints 
            legs(j) = plot3(RN(sum(nNP2(1:j))+1:sum(nNP2(1:j+1)),1), RN(sum(nNP2(1:j))+1:sum(nNP2(1:j+1)),2),RN(sum(nNP2(1:j))+1:sum(nNP2(1:j+1)),3),...
                markersymbols{j}, 'MarkerSize',markersizes{j},'Color',markercolours{j});
            h = legs(j);
        end

        if plotSD
            MU = mus(sum(nNM(1:j))+1:sum(nNM(1:j+1)),:);
            SIG = sigmas(j,:)*numberSD;

            for m = 1:size(MU,1)
                [x, y, z] = ellipsoid(MU(m,1),MU(m,2),MU(m,3),SIG(1),SIG(2),SIG(3),20);
                h = surfl(x, y, z);
                set(h,'edgecolor','black','facealpha',0,'EdgeAlpha',ealpha,'linewidth',1)
            end
        end
            
    end
    
    % Build legend string
    legendStr = cell(J,1);
    for j=1:J
        legendStr{j,1} = sprintf('Class %d points',j); 
    end
    
    legend(legs,legendStr,'Location','NorthEastOutside')
end

if plotPoints || plotSD
    if K >= 3
        view(defaultview)
        zoom(defaultzoom)
    end
    hold off;
end

%%%%%%%%%%%%%%% output files %%%%%%%%%%%%%%%%%%%%%
outputFile = sprintf('Synthetic-Gaussian-K-%d-sig-%.3f-modes%d', K, sigma,max(nOfModes));

if outputArff || outputPDF || outputSVG
    disp(sprintf('Output files base name: %s', outputFile));
end


if outputArff
    outputFileArff = [outputArffDir filesep outputFile '.arff'];
    
    fid = fopen(outputFileArff,'w');
    
    fprintf(fid,'%% Synthetic data generated with the following parameters: \n');
    fprintf(fid,'%%  - Number of dimensions (k): %d\n', K);
    fprintf(fid,'%%  - Variance (sigma): %.4f\n', sigma);
    fprintf(fid,'%%  - Number of modes: %s\n', mat2str(nOfModes));
    fprintf(fid,'%% \n\n');
    
    fprintf(fid,'@relation %s', sprintf('Synthetic-Gaussian-K-%d-sig-%f\n', K, sigma));
    
    for k = 1:K
        fprintf(fid,'@attribute X%d numeric\n', k);
    end
    
    aux = '@attribute Y {1,';
        
    for j = 2:J-1
        aux = [aux sprintf('%d,',j)];
    end
    
    aux = [aux sprintf('%d}\n\n',J)];

    fprintf(fid,sprintf('%s', aux));
    
    fprintf(fid,'@data\n');
    
    for i=1:size(DataSet.patterns,1)
        s = sprintf('%f', DataSet.patterns(i,1));
        for k=2:K
            s = sprintf('%s,%f',s, DataSet.patterns(i,k));
        end
        s = sprintf('%s,%d', s,DataSet.targets(i,1));
        fprintf(fid,s);
        fprintf(fid,'\n');
    end
    
    fclose(fid);
end

if outputDat
    dlmwrite([outputDatDir filesep outputFile '.dat'], [DataSet.patterns DataSet.targets], 'delimiter', ' ');
end
    
if outputPDF || outputSVG% || outputPNG
    if K >= 3
        if K>3 
            disp('Number of dimensions > 3, only the 3 first dimensions are represented')
        end
        if outputPDF
            saveas(h, [outputPicsDir filesep outputFile '.pdf']);
        end
        if outputSVG
            plot2svg([outputPicsDir filesep outputFile '.svg']);
        end
        % If you want to export to PNG uncomment the following lines, but
        % consider that this code can cause issues with some system configuration. 
        % Issues can be restart of X in Linux systems, which can make you to loose data
        %
        %if outputPNG
        %    export_fig([outputPicsDir filesep outputFile '.png'],'-png','-m3');
        %end
    elseif K <= 2
        if outputPDF
            export_fig([outputPicsDir filesep outputFile '.pdf'],'-pdf','-transparent');
        end
        if outputSVG
            plot2svg([outputPicsDir filesep outputFile '.svg']);
        end
        % If you want to export to PNG uncomment the following lines, but
        % consider that this code can cause issues with some system configuration. 
        % Issues can be restart of X in Linux systems, which can make you to loose data
        %
        %if outputPNG
        %    export_fig([outputPicsDir filesep outputFile '.png'],'-png','-m3');
        %end
    end
end


end
