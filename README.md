
Matlab scripts that generate stimuli for pitch discrimination tasks for Aphasia study.
Both pure tones and complex tones (complex tone generation scripts form Malinda McPherson).

Main Matlab script to generate the stimuli is: GenerateStimuli.m
Can run sub-section by subsection given that all other scripts are in Matlab path.
Should change directory paths in the code.

Matlab structures:
dfs = [0.125 0.25 0.5 1 2]; are the frequency intervals in semitones, spaced logarithmically.
There are 5 dfs corresponding to 5 conditions.
nCond = 5; is the number of conditions.

f1s are the frequencies of first tones: a vector of nf1 values.
nf1 = 7; is the number of f1 values.
log-uniform distribution between 200 to 400 Hz

f2s are the frequencies of the second tones. A matrix of nf1 x nCond. Given that the first tone was f1 which appears in index i in f1s, 
and the condition of frequency interval is df which appears in column j in dfs, the frequency of the
second tone is - f2 = f2s(i,j)

.wav sound files are named due to this convention f1_i and f2_i_cj :
generally 
f1_'f1 index'_'pure' or 'complex'. wav
f2_'f1 index'_c 'condition index' _ 'pure' or 'complex' . wav

e.g.
f1_1_pure.wav
f1_2_pure.wav
...
and
f2_1_c1_pure.wav
...

and the same for complex





