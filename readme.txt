{\rtf1\ansi\ansicpg1252\cocoartf1504\cocoasubrtf840
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\b\fs24 \cf0 \ul \ulc0 Pitch discrimination task
\b0 \ulnone \
\
\ul Matlab structures:\
\ulnone \
dfs = [0.125 0.25 0.5 1 2]; are the frequency intervals in semitones, spaced logarithmically.\
There are 5 dfs corresponding to 5 conditions.\
nCond = 5; is the number of conditions.\
\
nf1 = 7; is the number of f1 values\
\
f1 are the first tones: a vector of nf1 values\
log-uniform distribution 200 to 400 Hz\
\
f2 are the second tones. A matrix of nf1xnCond. Given that the first one was f1 due to the row number, and the condition of frequency interval is c due to the columns, the frequency of the second tone is -  f2=f2s(f1,c)\
\
.wav sound files are named accordingly - \
f1_1\
f1_2\
\'85\
\
and \
f2_1_c1\
\'85\
generally f2_number of f1_condition number\
 }