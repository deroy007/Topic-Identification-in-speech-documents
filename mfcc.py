import os
import librosa
import numpy as np
count=0
l=[]
total=[]
file_list=[]
k=0
for filename in os.listdir(os.getcwd()):
	file_list.append(filename)
final_file_list=[]
for filename in file_list:
	if filename=="mfcc.py":
		k=1
	elif filename=="transfinal":
		k=1
	elif filename=="test.py":
		k=1
	elif filename=="testk.py":
		k=1
	elif filename=="MFCC":
		k=1
	elif filename=="posterior":
		k=1
	elif filename=="WAV":
		k=1
	else:
		final_file_list.append(filename)

for filename in final_file_list:
	y, sr = librosa.load(filename)
	mfcc = librosa.feature.mfcc(y=y, sr=sr,n_mfcc=13)
	mfcc_delta = librosa.feature.delta(mfcc)
	mfcc_delta_delta=librosa.feature.delta(mfcc,order=2)
	for j in range(0,len(mfcc)):
		#print len(mfcc)
		for i in range(0,len(mfcc[count])):
			l.append(mfcc[count][i])
#			print mfcc[count][i]
		
		for i in range(0,len(mfcc_delta[count])):
			l.append(mfcc[count])
		for i in range(0,len(mfcc_delta_delta[count])):
			l.append(mfcc[count])

		print len(l)
		count+=1
		l=[]
		total.append(l)
		#print len()
	
	#print count
	count=0
	
	#print mfcc_delta.shape
	#print mfcc_delta_delta.shape 
	
	#count+=1
	#print mfcc_delta
	#print mfcc_delta_delta
total=np.array(total)
print total
