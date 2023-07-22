% Copyright Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay, 2020

if verLessThan('matlab','8.1')
    cxxFlags = ['CXXFLAGS="-std=c++17" ' ...
        'CXXLIBS="\$CXXLIBS -lc++" ' ]; % This flag is required on some platforms, but must be commented on others...
    outputFlag = '-o ';
else
    cxxFlags = 'CXXFLAGS="-std=c++17" ';
    outputFlag = '-output ';
end

compileHFM = @(binary_Dir,name,flag) eval(['mex ' ...
    outputFlag 'HamiltonGeodesicDistanceComputation' ' HamiltonGeodesicDistanceComputation.cxx' ...
    ' COMPFLAGS="/std:c++17"' ... % needed on windows 
    ' MACOSX_DEPLOYMENT_TARGET=12.3'...% needed on macOS to remove warningw. 
    ' -outdir ' binary_Dir ...
    ' ' cxxFlags  '-D' flag ...
    ' -I' '../GeoMetrics' ...
    ' -I' '../JMM_CPPLibs' ...
    ' -I' 'FastMarchingBase' ...
    ]);


 compileModelsHFM = @(binary_Dir,modelNames) ...
    cellfun(@(name) compileHFM(binary_Dir,name,['ModelName=' name]), modelNames);
standardModelsHFM = {'generalModels'};


fprintf(['\nPlease execute the function compileModelsHFM(binary_Dir,standardModelsHFM) to build \n'...
'the Hamilton Fast Marching executables in directory binary_Dir, \n' ...
'In case of need, replace standardModelsHFM with experimentalModelsHFM, or any list of desired models.\n']);

%For me : binary_Dir = '/Users/mirebeau/Dropbox/Programmes/MATLAB/MexBin';
% or : binary_Dir = '/Users/dachen/Desktop/Binary';

