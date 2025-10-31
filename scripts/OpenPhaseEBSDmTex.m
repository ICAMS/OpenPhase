%
% This file is a part of the OpenPhase software project.
% For more details visit www.openphase.de
%
% Created:    2014
%
% Authors:    Philipp Engels
%
% Copyright (c) 2009-2014 Interdisciplinary Centre for Advanced Materials
%              Simulation (ICAMS). Ruhr-Universitaet Bochum. Germany
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% This script was tested for mtex-4.0.beta1. Run startup_mtex.m in advance.
% Note: Script will search for ebsd data in the directory /RawData.

clc; clear;

%% Specify input type

plotinterface = false;
inputdata = 'euler';
%inputdata = 'bunge';

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {crystalSymmetry('m-3m','mineral','Interface','color', 'black'),...
      crystalSymmetry('m-3m','mineral','Cubic1','color', 'blue'),...
      crystalSymmetry('m-3m','mineral','Cubic2','color', 'red')};

% plotting convention
plotphase = 'Cubic1';
plotmarkersmax = 500;

setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');
plotx2east;

%% Create figure loop

for tStep = 0:0
    index = sprintf('%08d', tStep);
    %index = '001';
    %% Specify File Names
    subfilename = [pwd, '/RawData/ebsd_',index];
    fname = [subfilename,'.dat'];
    disp(['ebsd_', num2str(index)]);

    %% Import the Data

    % create an EBSD variable containing the data
    if strcmp(inputdata,'quaternion')  
       ebsd = loadEBSD_generic(fname,'CS', CS,...
      'ColumnNames', {'Index' 'Phase' 'X' 'Y' 'Z' 'Quat real' 'Quat i' 'Quat j' 'Quat k'},...
      'Quaternion');
    else
       ebsd = loadEBSD_generic(fname,'CS', CS,...
       'ColumnNames', {'Index' 'Phase' 'X' 'Y' 'Z' 'phi1' 'Phi' 'phi2'},...
       'Columns', [1 2 3 4 5 6 7 8], 'Bunge', 'radians');
    end
      
    % Calculate grains (3D support only implemented in old mtex
    % versions...)
    %grains = calcGrains(ebsd,'angle',10.0*degree);

    %% Plot

    %plot(ebsd('plotphase'))
    %hold on
    %plot(grains,'linewidth',1.5);
    % 
    % Calculate odf, variant 1
    %psi = calcKernel(grains('Al').meanOrientation);
    %odf = calcODF(ebsd('Al').orientations,'kernel',psi);
    % Calculate odf, variant 2
    if (plotinterface == true)
      odf = 0.5*calcODF(ebsd(plotphase).orientations) + 0.5*calcODF(ebsd('Interface').orientations);
    else
      odf = calcODF(ebsd(plotphase).orientations); 
    end
    
    % Plot of microstructure (only 2D in mtex4.0)
    ori =  [Miller(0,0,1,ebsd(plotphase).CS)...
            Miller(0,1,1,ebsd(plotphase).CS)...
            Miller(1,1,1,ebsd(plotphase).CS)];
    plotPDF(odf, ori, 'antipodal','grid','grid_res',15*degree);
    %mtexColorMap white2black
    %setColorRange('equal','all');
    setColorRange([0 15], 'current');
    %CLim(gcm,'equal');
    colorbar(gcm,'location','South')

    hold all
    if length(ebsd(plotphase).orientations) < plotmarkersmax
        nplotmarkers = 'all';
    else
        nplotmarkers = plotmarkersmax;
    end
    plotPDF(ebsd(plotphase).orientations, ori,'antipodal', 'MarkerSize',2,...
        'points', nplotmarkers, 'MarkerColor','black');
    if (plotinterface == true)
    plotPDF(ebsd('Interface').orientations, ori,'antipodal', 'MarkerSize',2,...
    'points', nplotmarkers, 'MarkerColor','red');
    end
%     plotPDF(ebsd('Interface').orientations, ori,'antipodal', 'MarkerSize',2,...
%         'points', nplotmarkers, 'MarkerColor','red');
    annotate([xvector, yvector, zvector], 'label', {'x', 'y', 'z'},...
        'BackgroundColor', 'w');
    hold off;
    export_fig(gcf,subfilename,'-png','-r300', '-nocrop'); 

    % Plot of inverse pole figure
%     plotIPDF(odf, xvector, 'antipodal', 'grid',...
%         'points',1000, 'MarkerSize',3);
%     plotIPDF(ebsd('Al').orientations, xvector, 'antipodal', 'grid',...
%        'points',1000, 'MarkerSize',3);
    %mark crystal axes
    %hold all
    %plot(Miller(0,0,1,ebsd(plotphase).orientations),'symmetrised','labeled');


end