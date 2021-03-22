# Voice Recognition based digit prediction
 Voice Recognition based digit prediction, developed using LPC analsis and HMM.
Deveopment:
Based on Cepstral coeffecients of phonemes codebook is generated and using it each word is divided into phonemes.
HMM is used to learn the underlying structure of each word and based on it model for each word is generated.
For each word in dataset a HMM model is there. 
For prediction, the model which gives maximum probability for a given observation sequence, corresponding word is predicted.
Accuracy achieved:
Prerecorded voice: 97%
Live recording: 90%
