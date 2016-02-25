import json
import math
import sys
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr

GroundTruthPath = sys.argv[1]
EstimationPath = sys.argv[2]

GroundTruthFile = open(GroundTruthPath, 'r')
EstimationFile = open(EstimationPath, 'r')
PsiGroundTruth = json.load(GroundTruthFile)
PsiEstimation = json.load(EstimationFile)

num = 0
tot = 0

prec = 0.01
avg = 0.0
maxv = 0.0
maxa = 0.0
maxb = 0.0
la = []
lb = []

for i in range(len(PsiGroundTruth)):
    #print('==' + str(i))
    for j in range(len(PsiGroundTruth[i])):
        a = PsiGroundTruth[i][j]
        b = PsiEstimation[i][j]
        la.append(a)
        lb.append(b)
        tot += 1
        #print('{0:.4f}'.format(a) + '\t' + '{0:.4f}'.format(b), end = '')
        if abs(a - b) < prec:
            num += 1
            #print('')
        else:
            #print('\t***')
            pass
        if abs(a - b) > maxv:
            maxv = abs(a - b)
            maxa = a
            maxb = b
        avg += (a - b) ** 2

avg = (avg / tot) ** 0.5

print('*** ' + '{0:.2f}'.format(num/tot*100) + '% PSI was correct under the percision ' + str(prec) + ' ***')
print('*** RSME is ' + '{0:.4f}'.format(avg)  + ' ***')
print('*** Pearson Correlation ratio is ' + '{0:.4f}'.format(pearsonr(la, lb)[0]) + ' ***')
print('*** Spearman Correlation ratio is ' + '{0:.4f}'.format(spearmanr(la, lb)[0]) + ' ***')
print('*** Maximum error is ' + '{0:.4f}'.format(maxv) + ' where '
        + '{0:.4f}'.format(maxa) + ' vs ' + '{0:.4f}'.format(maxb) + ' ***')
