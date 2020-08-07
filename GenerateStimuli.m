%% stimuli for pitch discrimination task:
%(similar to McPherson and McDermott 2018)
cd('/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Aphasia/2ToneTask/forHanna')
%% conditions: frequency intervals, in semitones
dS = [0.125 0.25 0.5 1 2];%difference in semitones
nCond = length(dS);
%% f1: the 1st out of the 2 tones
% logarithmically spaced between 200 and 400 Hz
nf1 = 7; % how many possible f1?
f1s = exp(linspace(log(200),log(400),nf1));
%% f2: all possible 2nd tones
% f2s is a matrix of nf1 x nCond. 
% The 2nd tone is selected with respect to the first tone (rows) and the condition (columns) 
f2s = nan(nf1,nCond);
for fi=1:nf1
    for cj=1:nCond
        f1=f1s(fi);M1=freq2MIDIcont(f1);M2=M1+dS(cj);
        f2s(fi,cj) = MIDI2freq(M2);
    end
end
%% synth pure tones
%params
    fs=44100;
    silence = [0 0];
    ramp = [30 30];%ms
    dur = 0.4;%s
    reduction = 0.5;
   
    %% f1
    soundwaves1=nan(nf1,fs*dur+1);
    for fi=1:nf1
        f = f1s(fi);
        [ data, t ] = SynthPureTone( f, dur, ramp, silence, reduction, fs, true );
        data=rmsnorm(data);
        soundwaves1(fi,:) = data;
        plot(t,data)
        soundsc(data,fs)  
        pause
        audiowrite(['f1pure' filesep 'f1_' num2str(fi) '.wav'],data,fs)
    end
    save(['f1pure' filesep 'f1all'],'soundwaves1')
    %% f2
    soundwaves2=nan(nf1,nCond,fs*dur+1);
    for fi=1:nf1
        hold off
        for cj=1:nCond
            f = f2s(fi,cj);
            [ data, t ] = SynthPureTone( f, dur, ramp, silence, reduction, fs, true );
            data=rmsnorm(data);
            plot(t,data)
            %pwelch(data,[],[],[],fs)
            soundsc(data,fs)
            pause
            audiowrite(['f2pure' filesep 'f2_' num2str(fi) '_c' num2str(cj) '.wav'],data,fs)
            soundwaves2(fi,cj,:) = data;
            hold all
        end
    end 
    save(['f2pure' filesep 'f2all'],'soundwaves2')
%% synth complex tones (Malinda script)
    %params
    harm_nums = 1:100;
    jitt_amt = 0;
    jitt = 0; % 2 to generate a single note with jitter (1 requires JitterString input)
    dist = 0;% 0 for exponential decay of spectral envelope
    JitterString = [];
    centroid = [];    
    %% f1
    
    soundwaves1=nan(nf1,fs*dur+1);
    for fi=1:nf1
        f = f1s(fi);
        data = zeros(1,size(soundwaves1,2));
        [signal, jitter_m] = generate_singlenote_vary_envelope_jitter_randphase_Tamar(f, harm_nums, jitt_amt,jitt, dur, fs, dist,JitterString, centroid);
        data(1:length(signal))=signal;
        t=0:1/fs:length(data)/fs-1/fs;
        data = hann(data,ramp(1),fs);
        data = rmsnorm(data);    
        soundwaves1(fi,:) = data;
        plot(t,data)
        soundsc(data,fs)  
        pause
        audiowrite(['f1complex' filesep 'f1_' num2str(fi) '.wav'],data,fs)
    end
    save(['f1complex' filesep 'f1all'],'soundwaves1')
    
    %% f2
    
    soundwaves2=nan(nf1,nCond,fs*dur+1);
    for fi=1:nf1
        hold off
        for cj=1:nCond
            f = f2s(fi,cj);
            data = zeros(1,size(soundwaves2,3));
            [signal, jitter_m] = generate_singlenote_vary_envelope_jitter_randphase_Tamar(f, harm_nums, jitt_amt,jitt, dur, fs, dist,JitterString, centroid);
            data(1:length(signal))=signal;
            t=0:1/fs:length(data)/fs-1/fs;
            data=rmsnorm(data);
            plot(t,data)
            %pwelch(data,[],[],[],fs)
            soundsc(data,fs)
            pause
            audiowrite(['f2complex' filesep 'f2_' num2str(fi) '_c' num2str(cj) '.wav'],data,fs)
            soundwaves2(fi,cj,:) = data;
            hold all
        end
    end 
    save(['f2complex' filesep 'f2all'],'soundwaves2')

    