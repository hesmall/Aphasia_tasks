function [ MIDI ] = freq2MIDIcont( f )
% from frequency in Hz to MIDI note number
    MIDI = 69 + (12*log2(f/440));

end

