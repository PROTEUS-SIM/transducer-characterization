function [FrequencyMHz1, InsertionLossdB, StandardDeviationdB, ...
    ReturnLossdB, FrequencyMHz2, GroupDelaynsec] = ...
    read_filter_data(filename)
%READ_FILTER_DATA reads the data from a Mini-Circuits analog filter data
%sheet.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [3, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = [...
    "FrequencyMHz1", ...
    "InsertionLossdB", ...
    "StandardDeviationdB", ...
    "ReturnLossdB", ...
    "FrequencyMHz2", ...
    "GroupDelaynsec"];

opts.VariableTypes = repmat("double",1,6);

% Specify file level properties
opts.ExtraColumnsRule = "ignore";

% Import the data
tbl = readtable(filename, opts);

% Convert to output type
FrequencyMHz1       = tbl.FrequencyMHz1;
InsertionLossdB     = tbl.InsertionLossdB;
StandardDeviationdB = tbl.StandardDeviationdB;
ReturnLossdB        = tbl.ReturnLossdB;
FrequencyMHz2       = tbl.FrequencyMHz2;
GroupDelaynsec      = tbl.GroupDelaynsec;

end