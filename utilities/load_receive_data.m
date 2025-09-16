function [y,Fs] = load_receive_data(filename)
% LOAD_RECEIVE_DATA reads the voltage receive data of the transducer from
% file filename.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

load(filename,'RcvData','Trans','Receive','TGC')

% Positive transmit pulse
startSample = Receive(1).startSample;
endSample   = Receive(1).endSample;

y1 = RcvData{1,1}(startSample:endSample,Trans.ConnectorES,:);
y1 = double(y1);
y1 = mean(y1,3);

% Negative transmit pulse
startSample = Receive(2).startSample;
endSample   = Receive(2).endSample;

y2 = RcvData{1,1}(startSample:endSample,Trans.ConnectorES,:);
y2 = double(y2);
y2 = mean(y2,3);

% Take the average of the positive and the negative transmit pulse
y = (y1 - y2)/2;

% Digitally reverse TGC
% From the Verasonics manual: The 40 dB TGC range of the Vantage system is
% roughly log-linear over this 0 to 1023 range of the control points.
TGC_range = 40; % dB

gain_dB = TGC.CntrlPts(1)/1023*TGC_range;
% gain_dB = gain_dB + RcvProfile.PgaGain + RcvProfile.LnaGain;
y = y*10^(-gain_dB/20);

Fs = Receive(1).decimSampleRate*1e6; % Sampling rate in Hz

end