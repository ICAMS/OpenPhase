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

clear;
clc;

fprintf('HexToReg EBSD parser\n')
fprintf('--------------------\n')

% -------------------------------------------------------------------------
% Input
filename = 'PR21-6008h-0_ strain.txt';
filexcoordcolumn = 4;
fileycoordcolumn = 5;
filephasecoordcolumn = 8;
dx = 0.1;
% -------------------------------------------------------------------------

fprintf(['Read from: ',filename, '\n'])
fprintf(['Use discretization: ', num2str(dx), '\n']);
fprintf(['Read x coordinate from col: ', num2str(filexcoordcolumn), '\n']);
fprintf(['Read y coordinate from col: ', num2str(fileycoordcolumn), '\n']);
fprintf(['Read phase index from col: ', num2str(filephasecoordcolumn), '\n']);

% -------------------------------------------------------------------------

rawdata = importdata(filename);

x = rawdata(:,filexcoordcolumn);
y = rawdata(:,fileycoordcolumn);
z = rawdata(:,filephasecoordcolumn);

maxX = max(x);
maxY = max(y);

Nx = round(maxX / dx);
Ny = round(maxY / dx);

fprintf(['Used Nx: ', num2str(Nx), '\n']);
fprintf(['Used Ny: ', num2str(Ny), '\n']);

[Xq,Yq] = meshgrid(0:maxX/Nx:max(x),0:maxY/Ny:max(y));

intval = griddata(x,y,z,Xq,Yq);
intvalint = round(intval);

Res = ones(Nx*Ny,3);

for j = 1:Ny,
    for i = 1:Nx,
        loopvar = i + (j-1)*Nx;
        
        Res(loopvar,1) = Xq(1,i);
        Res(loopvar,2) = Yq(j,1);
        Res(loopvar,3) = intvalint(j,i);
    end
end

outputfilename = [filename, '_squared.dat'];
dlmwrite(outputfilename, Res, ' ');

fprintf(['Written to: ',outputfilename, '\n'])
