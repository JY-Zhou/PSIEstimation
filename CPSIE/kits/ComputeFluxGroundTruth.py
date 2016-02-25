import json
import re

ProfileFile = open('./readsFlux.pro', 'r')
confFile = open('./GenerationConfig.json', 'r')

config = json.load(confFile)

NG = config['NG']
NE = config['NE']
iso = config['Iso']

Psi = []
tot = []

for g in range(NG):
    Psi.append([])
    tot.append(0.0)
    for e in range(NE[g]):
        Psi[g].append(0.0)

for line in ProfileFile:
    substr = line[:-1].split('\t')
    trxname = substr[1]
    isoinfo = re.findall(r'\d+', trxname)
    geneid = int(isoinfo[0])
    isoid = int(isoinfo[1])

    nmRNA = float(substr[5])
    ncDNA = float(substr[7])

    tot[geneid] += nmRNA
    for e in iso[geneid][isoid]:
        Psi[geneid][e] += nmRNA

for g in range(NG):
    for e in range(NE[g]):
        Psi[g][e] /= tot[g]

output = open('./PsiGroundTruthFlux.json', 'w')
json.dump(Psi, output, indent = 4)
