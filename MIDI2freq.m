function [ f ] = MIDI2freq( MIDI )
%from MIDI number to frequency in Hz
f  =  2^((MIDI - 69)/12)*440;

end

