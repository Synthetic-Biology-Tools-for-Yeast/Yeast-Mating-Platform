
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

raw = pd.read_csv('results9L1all.csv')

tresholdstart=50
treshold=10

df=raw.copy(deep=True)
df=df.rename(columns={"Frequency": "count"})
df=df[df['count']>=tresholdstart]

#Fg --- results

#first dataframe
Fg1Ago1 = pd.read_csv('results1FgAgo1.csv')
Fg1Ago1=Fg1Ago1.rename(columns={"Frequency": "Frequency1"})
Fg1Ago1=Fg1Ago1[Fg1Ago1['Frequency1']>=treshold]

#second dataframe
Fg2Ago2 = pd.read_csv('results2FgAgo2.csv')
Fg2Ago2=Fg2Ago2.rename(columns={"Frequency": "Frequency2"})
Fg2Ago2 = Fg2Ago2[Fg2Ago2['Frequency2']>=treshold]

#merge 1 and 2 on DNA sequence
df1DNA = pd.merge(Fg1Ago1, Fg2Ago2, how='inner', left_on='Sequence', right_on='Sequence')
#average
df1DNA["AvFrequency"]=(df1DNA["Frequency1"]+df1DNA["Frequency2"])/2
df1DNA.sort_values(by='AvFrequency', ascending=False, inplace=True)

#add translation
translationlist=[]
for i in df1DNA["Sequence"]:
    seq=Seq(i)
    translation=seq.translate()
    translationlist.append(str(translation))
df1DNA["translation"]=translationlist

#merge results with starting library on DNA level
FgDNA = pd.merge(df1DNA, df, how='inner', left_on='Sequence', right_on='Sequence')
FgDNA['normDNA']=(FgDNA['AvFrequency']-FgDNA['count'])/FgDNA['count']
FgDNA.sort_values(by='normDNA', ascending=False, inplace=True)

#summ coding sequenceses that have the same AA
colIwant=['AvFrequency','count','normDNA']
summed_df1 = FgDNA.groupby('translation')[colIwant].sum().reset_index()
summed_df1['norm']=(summed_df1['AvFrequency']-summed_df1['count'])/summed_df1['count']
summed_df1.sort_values(by='norm', ascending=False, inplace=True)


#Bb --- results

#first dataframe
Bb3Ago1 = pd.read_csv('results3BbAgo1.csv')
Bb3Ago1=Bb3Ago1.rename(columns={"Frequency": "Frequency1"})
Bb3Ago1=Bb3Ago1[Bb3Ago1['Frequency1']>=treshold]

#second dataframe
Bb4Ago2 = pd.read_csv('results4BbAgo2.csv')
Bb4Ago2=Bb4Ago2.rename(columns={"Frequency": "Frequency2"})
Bb4Ago2 = Bb4Ago2[Bb4Ago2['Frequency2']>=treshold]

#merge 1 and 2 on DNA sequence
df2DNA = pd.merge(Bb3Ago1, Bb4Ago2, how='inner', left_on='Sequence', right_on='Sequence')
#average
df2DNA["AvFrequency"]=(df2DNA["Frequency1"]+df2DNA["Frequency2"])/2
df2DNA.sort_values(by='AvFrequency', ascending=False, inplace=True)

#add translation
translationlist=[]
for i in df2DNA["Sequence"]:
    seq=Seq(i)
    translation=seq.translate()
    translationlist.append(str(translation))
df2DNA["translation"]=translationlist

#merge results with starting library on DNA level
BbDNA = pd.merge(df2DNA, df, how='inner', left_on='Sequence', right_on='Sequence')
BbDNA['normDNA']=(BbDNA['AvFrequency']-BbDNA['count'])/BbDNA['count']
BbDNA.sort_values(by='normDNA', ascending=False, inplace=True)

#summ coding sequenceses that have the same AA
colIwant=['AvFrequency','count','normDNA']
summed_df2 = BbDNA.groupby('translation')[colIwant].sum().reset_index()
summed_df2['norm']=(summed_df2['AvFrequency']-summed_df2['count'])/summed_df2['count']
summed_df2.sort_values(by='norm', ascending=False, inplace=True)


#Bc --- results

#first dataframe
Bc5Ago1 = pd.read_csv('results5BcAgo1.csv')
Bc5Ago1=Bc5Ago1.rename(columns={"Frequency": "Frequency1"})
Bc5Ago1=Bc5Ago1[Bc5Ago1['Frequency1']>=treshold]

#second dataframe
Bc6Ago2 = pd.read_csv('results6BcAgo2.csv')
Bc6Ago2=Bc6Ago2.rename(columns={"Frequency": "Frequency2"})
Bc6Ago2 = Bc6Ago2[Bc6Ago2['Frequency2']>=treshold]

#merge 1 and 2 on DNA sequence
df3DNA = pd.merge(Bc5Ago1, Bc6Ago2, how='inner', left_on='Sequence', right_on='Sequence')
#average
df3DNA["AvFrequency"]=(df3DNA["Frequency1"]+df3DNA["Frequency2"])/2
df3DNA.sort_values(by='AvFrequency', ascending=False, inplace=True)

#add translation
translationlist=[]
for i in df3DNA["Sequence"]:
    seq=Seq(i)
    translation=seq.translate()
    translationlist.append(str(translation))
df3DNA["translation"]=translationlist

#merge results with starting library on DNA level
BcDNA = pd.merge(df3DNA, df, how='inner', left_on='Sequence', right_on='Sequence')
BcDNA['normDNA']=(BcDNA['AvFrequency']-BcDNA['count'])/BcDNA['count']
BcDNA.sort_values(by='normDNA', ascending=False, inplace=True)

#summ coding sequenceses that have the same AA
colIwant=['AvFrequency','count','normDNA']
summed_df3 = BcDNA.groupby('translation')[colIwant].sum().reset_index()
summed_df3['norm']=(summed_df3['AvFrequency']-summed_df3['count'])/summed_df3['count']
summed_df3.sort_values(by='norm', ascending=False, inplace=True)


#Fo --- results

#first dataframe
Fo7Ago1 = pd.read_csv('results7FoAgo1.csv')
Fo7Ago1=Fo7Ago1.rename(columns={"Frequency": "Frequency1"})
Fo7Ago1=Fo7Ago1[Fo7Ago1['Frequency1']>=treshold]

#second dataframe
Fo8Ago2 = pd.read_csv('results8FoAgo2.csv')
Fo8Ago2=Fo8Ago2.rename(columns={"Frequency": "Frequency2"})
Fo8Ago2 = Fo8Ago2[Fo8Ago2['Frequency2']>=treshold]

#merge 1 and 2 on DNA sequence
df4DNA = pd.merge(Fo7Ago1, Fo8Ago2, how='inner', left_on='Sequence', right_on='Sequence')
#average
df4DNA["AvFrequency"]=(df4DNA["Frequency1"]+df4DNA["Frequency2"])/2
df4DNA.sort_values(by='AvFrequency', ascending=False, inplace=True)

#add translation
translationlist=[]
for i in df4DNA["Sequence"]:
    seq=Seq(i)
    translation=seq.translate()
    translationlist.append(str(translation))
df4DNA["translation"]=translationlist

#merge results with starting library on DNA level
FoDNA = pd.merge(df4DNA, df, how='inner', left_on='Sequence', right_on='Sequence')
FoDNA['normDNA']=(FoDNA['AvFrequency']-FoDNA['count'])/FoDNA['count']
FoDNA.sort_values(by='normDNA', ascending=False, inplace=True)

#summ coding sequenceses that have the same AA
colIwant=['AvFrequency','count','normDNA']
summed_df4 = FoDNA.groupby('translation')[colIwant].sum().reset_index()
summed_df4['norm']=(summed_df4['AvFrequency']-summed_df4['count'])/summed_df4['count']
summed_df4.sort_values(by='norm', ascending=False, inplace=True)

#Plot

yoFg=summed_df1.head(10)

#fig, ax = plt.subplots(figsize =(16, 40))
fig, ax = plt.subplots()
#ax.set_xlim(left=0, right=160)
ax.barh(yoFg["translation"], yoFg["norm"])
#plt.xscale("log")
plt.title('Fg')

yoFg=summed_df2.head(10)

#fig, ax = plt.subplots(figsize =(16, 40))
fig, ax = plt.subplots()
#ax.set_xlim(left=0, right=160)
ax.barh(yoFg["translation"], yoFg["norm"])
plt.title('Bb')

yoFg=summed_df3.head(10)

#fig, ax = plt.subplots(figsize =(16, 40))
fig, ax = plt.subplots()
#ax.set_xlim(left=0, right=160)
ax.barh(yoFg["translation"], yoFg["norm"])
plt.title('Bc')


yoFg=summed_df4.head(10)

#fig, ax = plt.subplots(figsize =(16, 40))
fig, ax = plt.subplots()
#ax.set_xlim(left=0, right=160)
ax.barh(yoFg["translation"], yoFg["norm"])
plt.title('Fo')
