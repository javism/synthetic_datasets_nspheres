%SYNTHETIC_DATA_NSPHERES_CONF Configuration file example for SYNTHETIC_DATA_NSPHERES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) Javier Sánchez Monedero (jsanchezm at uco dot es)
% 
% This code implements the synthetic data generatod described in the following publicacion:
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

%%%%%%%%%%%%%%%%%%%%%%%% Gaussian distribution setup %%%%%%%%%%%%%%%%%%%%%%

% NOTE: N and nOfModes must have to have the same number of elements!!!!!

% Patterns distribution for each class
% and Number of modes for each class
% Example with 4 classes and 1 mode per class:
% N = [5;50;80;10];
% nOfModes = [1,1,1,1];
% Example with 7 classes and several number of modes:
N = [30;30];
nOfModes = [3,3];


% Number of dimensions. You can specify one or a set of dimensions so that
% several datasets will be generated, one per dimension in combination
% with the standard deviation values to use (Sset parameter).
% Example: Kset = [2, 10, 50, 100]; will generate one dataset for each 
% dimension in combination with Sset
Kset = [2];

% Variance value. You can specify a set of values in order to generate
% several datasets. 
% Example: Sset = [1/6,1/3,1/2];
Sset = [1/8,1/4,1/3];


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extra parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial value for the mean of the first Gaussian
muq = 0;

% Separation between the 'main' Gaussian of each class with the following
% class
alphaq = 1;

% Radius of the n-sphere in which hyper-surface there will be placed
% additional Gaussian functions for each class
lambdaq = 0.75;

% For multi-modal case, there should be generated a number of points
% in the n-sphere surface, and nOfModes(j) Gaussian functions will be centered in 
% those points, being j the class number
hyper_sphere_surface_points = 100;


%%%%%%%%%%%%%%%%%%%%% Select which information to plot %%%%%%%%%%%%%%%%%%%%

% Plot patterns
plotPoints = true;

% Plot 3 standard deviations, this is, the n-sphere of the isotropic
% Gaussian
plotSD = false;
numberSD = 3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% arff format for Weka
outputArff = true;

% raw format that can be read with Matlab's 'load' function
outputDat = true;

% PDF plot of the dataset
outputPDF = true;
% SVG plot of the dataset
outputSVG = true;

% UNSTABLE: PNG plot of the dataset
% If you want to export to PNG uncomment the lines at the end of synthetic_datasets_nspheres
% but consider that this code can cause issues with some system configuration. 
% Issues can be restart of X in Linux systems, which can make you to loose data
%outputPNG = false;

% Datasets output dir
outputArffDir = 'arff';

% Datasets output dir
outputDatDir = 'dat';

% Pics output dir
outputPicsDir = 'pics';


%%%%%%%%%%%%%%%%%%%%% Patterns colours and symbols setup %%%%%%%%%%%%%%%%%%
mucolor = [0.7,0.7,0.7];
softgray = [0.7,0.7,0.7];
darkgreen = [0.2,0.7,0.2];

% Sice of the patterns points
markersize = 5;

% colours+symbols combinations limited to 6 number of classes
markersymbols = {'+'; '.'; '^';'*';'d';'o';'v';'s';'>';'<'};
markercolours = {'r'; 'b'; darkgreen;'k';softgray;'m';'y';'r';'b';'g'};

% correct marker sizes
markersizes = {markersize;markersize*2.5;markersize;markersize;markersize...
                ;markersize;markersize;markersize;markersize;markersize};
            
% Default view angles and zoom for 3D
defaultview= [20,20];
defaultzoom = 1.2;

% Transparency for figures, actually 'EdgeAlpha' matlab's property value 
ealpha = 0.1;

% Plot resolution. The larger, the better-looking, but also the slower
RESOL = 10;
% Overrange: fraction of plot bordering the data
OR = 0.1;

% Histogram precision (for 1D plots)
histprecision = 40;
