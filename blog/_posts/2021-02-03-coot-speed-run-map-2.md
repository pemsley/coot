---
layout: post
title:  "Coot Speedrun Map 2"
date: Wed 3 Feb 14:16:44 GMT 2021
---

Coot Speedrun Map 2 is a step up from Map 1. The goal is to create, using de novo model-building, an accurate and full atomic
model for a cryo-EM reconstruction.

I have selected EMD-22898/7kjr which is SARS CoV-2 ORF 3a - an ion channel and potential therapeutic and vaccine target.
I think I recall seeing on twitter that the images for this reconstruction were collected in about 6 hours.

The starting point is

i) the map emd\_22898.map (get it from here: https://www.emdataresource.org/EMD-22898)

ii) the sequence

```
        10         20         30         40         50
                                          S LPFGWLIVGV 
        60         70         80         90        100
ALLAVFQSAS KIITLKKRWQ LALSKGVHFV CNLLLLFVTV YSHLLLVAAG 
       110        120        130        140        150
LEAPFLYLYA LVYFLQSINF VRIIMRLWLC WKCRSKNPLL YDANYFLCWH 
       160        170        180        190        200
TNCYDYCIPY NSVTSSIVIT SGDGTTSPIS EHDYQIGGYT EKWESGVKDC 
       210        220        230        240        250
VVLHSYFTSD YYQLYSTQLS TDTGVEHVTF FIYNKIVD
```

iii) the compound\_id for the ligand is PEE

iv) C2 symmetry has been applied along the Z axis to the reconstruction and hence should be reflected in the atomic model

The helices for the MSP are a bit marginal and will not be part of the scoring.

There will be a "Tool Assisted" version of this challenge - where one can start from the output of phenix.map\_to\_model. However,
can you complete the model with _Coot_ before phenix.map\_to\_model terminates?

I will need to write a progress-bar/timer widget to score the current model (I am not yet sure what the scoring system should be).

There is no closing date as yet - it will be several months from now.

Also note that this molecule is part of the 2021 Ligand Model Challenge (part of the EM Validation Challenges)

[EM Validation Challenges](https://challenges.emdataresource.org/?q=2021-launch-announce)

