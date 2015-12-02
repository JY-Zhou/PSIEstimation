import json
import random

EXONNUM_L = 5
EXONNUM_U = 10

EXONLEN_L = 148
EXONLEN_U = 250

ISONUM_L = 1
ISONUM_U = 4

ISOEXON_L = 3

NG = 10
NE = []
for g in range(NG):
    eVal = random.randint(EXONNUM_L, EXONNUM_U)
    NE.append(eVal)
print(NE)

L = []
for g in range(NG):
    Ltemp = []
    for e in range(NE[g]):
        lVal = random.randint(EXONLEN_L, EXONLEN_U)
        Ltemp.append(lVal)
    L.append(Ltemp)
print(L)

IsoNum = []
for g in range(NG):
    isonumVal = random.randint(ISONUM_L, ISONUM_U)
    IsoNum.append(isonumVal)
print(IsoNum)

Iso = []
for g in range(NG):
    print('=====')
    Isotemp = []
    for i in range(IsoNum[g]):
        temp = random.sample(list(range(NE[g])), random.randint(3, NE[g]))
        temp = sorted(temp)
        print(temp)