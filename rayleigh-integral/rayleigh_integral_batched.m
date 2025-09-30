function [source,sensor] = rayleigh_integral_batched(...
    source,sensor,medium,batchSize)
% RAYLEIGH_INTEGRAL_BATCHED computes rayleigh_integral_frequency in batches
% with each batchSize sensor points or source points.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

Nsource = size(source.points,1);       % Total number of source points
Nsensor = size(sensor.points,1);       % Total number of sensor points
Nfreq   = size(source.frequencies,2);  % Number of frequencies

if strcmp(medium.direction,'forward')
    i1 = 1:batchSize:Nsensor; % Start index for each batch
    i2 = i1 + batchSize - 1;  % End index for each batch
    i2(end) = Nsensor;
else
    i1 = 1:batchSize:Nsource; % Start index for each batch
    i2 = i1 + batchSize - 1;  % End index for each batch
    i2(end) = Nsource;
end

Nbatch = length(i1);      % Number of batches

if strcmp(medium.direction,'forward')
    sensor.pressures = zeros(Nsensor,Nfreq);
    sensorBatch = sensor;
elseif strcmp(source.type,'pressure')
    source.pressures = zeros(Nsource,Nfreq);
    sourceBatch = source;
else
    source.velocities = zeros(Nsource,Nfreq);
    sourceBatch = source;
end

% Progress message
fprintf('Computing Rayleigh integral\n');
fieldWidth = length(num2str(Nbatch));
fprintf('Evaluating batch ')
fprintf(repmat(' ',1,fieldWidth*2+1)) % Spaces to be overwritten

for m = 1:Nbatch
    
    fprintf(repmat('\b',1,fieldWidth*2+1)) % Overwrite previous characters
    fprintf(['%' num2str(fieldWidth) '.f/' ...
        '%' num2str(fieldWidth) '.f'],m,Nbatch) % Print progress
    
    I1 = i1(m); % Start index current batch
    I2 = i2(m); % End index current batch
    
    if strcmp(medium.direction,'forward')
        % Sensor struct with subset of the sensor points:
        sensorBatch.points = sensor.points(I1:I2,:);
        
        if isfield(sensor,'pressures')
            sensorBatch.pressures = sensor.pressures(I1:I2,:);
        end
        
        if isfield(sensor,'normal')
            sensorBatch.normal = sensor.normal(I1:I2,:);
        end
        
        if isfield(sensor,'weights')
            sensorBatch.weights = sensor.weights(I1:I2,:);
        end

        % Evaluate the rayleigh integral for the batch:
        [~,sensorBatch] = rayleigh_integral_frequency(...
            source, sensorBatch, medium);

        
    else
        % Source struct with subset of the source points:
        sourceBatch.points = source.points(I1:I2,:);
        
        if isfield(source,'velocities')
            sourceBatch.velocities = source.velocities(I1:I2,:);
        end
        
        if isfield(source,'pressures')
            sourceBatch.pressures = source.pressures(I1:I2,:);
        end
        
        if isfield(source,'normal')
            sourceBatch.normal = source.normal(I1:I2,:);
        end
        
        if isfield(source,'weights')
            sourceBatch.weights = source.weights(I1:I2,:);
        end

        % Evaluate the rayleigh integral for the batch:
        [sourceBatch,~] = rayleigh_integral_frequency(...
            sourceBatch, sensor, medium);
       
    end
    
    % Add the results to the complete set of source or sensor points:
    if strcmp(medium.direction,'forward')
        sensor.pressures(I1:I2,:)  = sensorBatch.pressures;
    elseif strcmp(source.type,'pressure')
        source.pressures(I1:I2,:)  = sourceBatch.pressures;
    else
        source.velocities(I1:I2,:) = sourceBatch.velocities;
    end
        
end

fprintf('\nDone\n')

end
