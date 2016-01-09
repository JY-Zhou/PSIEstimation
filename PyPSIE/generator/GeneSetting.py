import json
import random

EXONNUM_L = 5
EXONNUM_U = 10

EXONLEN_L = 148
EXONLEN_U = 250

ISONUM_L = 1
ISONUM_U = 4

ISOEXON_L = 3

EXPLV_L = 1
EXPLV_U = 2

NG = 50
print('===NG===')
print(NG)

NE = []
for g in range(NG):
    eVal = random.randint(EXONNUM_L, EXONNUM_U)
    NE.append(eVal)
print('===NE===')
print(NE)

L = []
for g in range(NG):
    Ltemp = []
    for e in range(NE[g]):
        lVal = random.randint(EXONLEN_L, EXONLEN_U)
        Ltemp.append(lVal)
    L.append(Ltemp)
print('===L===')
for x in L:
    print(x)

IsoNum = []
for g in range(NG):
    isonumVal = random.randint(ISONUM_L, ISONUM_U)
    IsoNum.append(isonumVal)
print('===IsoNum===')
print(IsoNum)

Iso = []
for g in range(NG):
    Isotemp = []
    for i in range(IsoNum[g]):
        temp = random.sample(list(range(NE[g])), random.randint(ISOEXON_L, NE[g]))
        temp = sorted(temp)
        while temp in Isotemp:
            temp = random.sample(list(range(NE[g])), random.randint(ISOEXON_L, NE[g]))
            temp = sorted(temp)
        Isotemp.append(temp)
    Iso.append(Isotemp)
print('===Iso===')
for x in Iso:
    print(x)
    
expLv = []
for g in range(NG):
    expLvtemp = []
    for i in range(IsoNum[g]):
        expLvtemp.append(random.randint(EXPLV_L, EXPLV_U))
    expLv.append(expLvtemp)
print('===expLv===')
for x in expLv:
    print(x)
    
depth = 1
readLength = 75

config = {'NG' : NG,
          'NE' : NE,
          'L' : L,
          'Iso' : Iso,
          'expLv' : expLv,
          'depth' : depth,
          'readLength' : readLength}

json.dump(config, open('../kits/GenerationConfig.json', 'w'), indent = 4)