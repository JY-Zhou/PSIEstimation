import json
from generator.ReferenceGenerator import ReferenceGenerator
from generator.ReadGenerator import ReadGenerator

confFile = open('../tricky/groundTruth.json', 'r')

groundTruth = json.load(confFile)
referenceGen = ReferenceGenerator(groundTruth['NG'],
                                  groundTruth['NE'],
                                  groundTruth['L'],
                                  groundTruth['Iso'])
referenceGen.work()
readGen = ReadGenerator(groundTruth['expLv'],
                        groundTruth['depth'],
                        groundTruth['readLength'])
readGen.work()