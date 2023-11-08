import json
import argparse
import re

parser = argparse.ArgumentParser()

parser.add_argument("-n", "--nameFile", required=True)
parser.add_argument("-e", "--enrichmentFile", required=True)
parser.add_argument("-c", "--coverageCountsFile", required=True)
parser.add_argument("-w", "--coverageWindowsFile", required=True)
parser.add_argument("-s", "--saturationMetricsFile", required=True)
parser.add_argument("-d", "--dedupMetricsFile", required=True)
parser.add_argument("-u", "--summaryFile", required=True)
parser.add_argument("-a", "--alignmentFile", required=True)
parser.add_argument("-t", "--thaliaFile", required=True)


args = parser.parse_args()
data = {}
with open(args.nameFile, "r") as nameFile:
    name = (nameFile.read().split())

with open(args.enrichmentFile, "r") as enrichFile:
    enrich = (enrichFile.read().split())

with open(args.coverageCountsFile, "r") as coverageCountsFile:
    counts = (coverageCountsFile.read().split())

with open(args.coverageWindowsFile,"r") as coverageWindowsFile:
    windows = (coverageWindowsFile.read().split())

with open(args.saturationMetricsFile,"r") as saturationMetricsFile:
    saturation = (saturationMetricsFile.read().split())

with open(args.dedupMetricsFile, "r+") as dedupFile:
    new_f = dedupFile.readlines()
    dedupFile.seek(0)
    for line in new_f:
        if "#" not in line:
            if line != '\n':
                dedupFile.write(line)
    dedupFile.truncate()

with open(args.dedupMetricsFile, "r+") as dedupFileNew:
    dedup = re.split('\t|\n', dedupFileNew.read())


with open(args.summaryFile, "r+") as summaryFile:
    new_f = summaryFile.readlines()
    summaryFile.seek(0)
    for line in new_f:
        if "#" not in line:
            if line != '\n':
                summaryFile.write(line)
    summaryFile.truncate()
with open(args.summaryFile, "r+") as summaryFileNew:
    summary = re.split('\t|\n', summaryFileNew.read())

with open(args.alignmentFile, "r+") as alignmentFile:
    new_f = alignmentFile.readlines()
    alignmentFile.seek(0)
    for line in new_f:
        if "#" not in line:
            if line != '\n':
                alignmentFile.write(line)
    alignmentFile.truncate()
with open(args.alignmentFile, "r+") as alignmentFileNew:
    alignment = re.split('\t|\n', alignmentFileNew.read())

with open(args.thaliaFile, "r") as thaliaFile:
    thalia = (thaliaFile.read().split())


#Outputting the specific entries in the file into the JSON file

data[name[0]] = name[1]

for x in range(0,13):
    data[enrich[0+x]] = enrich[14+x]

for x in range(0,2):
    data[counts[0+x]] = counts[3+x]

for x in range(0,5):
    data[windows[1+x]] = windows[6+x]

for x in range(0,4):
    data[saturation[0+x]] = saturation[5+x]

if isinstance(dedup, list) and len(dedup) >= 20:
    dedup[0] = "LIBRARY_0"
    for x in range(0, 10):
        data[dedup[x]] = dedup[10 + x]
    if data[dedup[9]] == "set_size":
        data[dedup[9]] = ""
    data[dedup[0]] = data[dedup[0]].replace(' ', '_')

if isinstance(summary, list) and len(summary) >= 30:
    summary[12] = "SAMPLE_0"
    summary[13] = "LIBRARY_1"
    summary[14] = "READ_GROUP_0"
    for x in range(0, 15):
        data[summary[x]] = summary[15+x]
    data[summary[0]] = data[summary[0]].replace(' ', '_')

if isinstance(alignment, list) and len(alignment) >= 54:
    alignment[24] = "SAMPLE_1"
    alignment[25] = "LIBRARY_2"
    alignment[26] = "READ_GROUP_1"
    for x in range(0, 27):
        data[alignment[x]] = alignment[27+x]

for x in range(0, 6):
    data[thalia[0+x]] = thalia[6+x]

with open('qc_metrics.json', 'w') as f:
    json.dump(data, f, indent=4)
