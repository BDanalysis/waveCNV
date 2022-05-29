import sys
import numpy as np

result_start = []
result_end = []
result_type = []
file = sys.argv[1]
with open(file, 'r') as f:
    for line in f:
        linestr = line.strip()
        linestrlist = linestr.split('\t')
        result_start.append(int(linestrlist[0]))
        result_end.append(int(linestrlist[1]))
        if linestrlist[2] == 'gain':
            result_type.append("gain")
        else:
            result_type.append("loss")



truth_start = []
truth_end = []
truth_type = []
with open("GroundTruth", 'r') as f:
    line = f.readline()
    for line in f:
        linestr = line.strip('\n')
        linestrlist = linestr.split('\t')
        truth_start.append(int(linestrlist[0]))
        truth_end.append(int(linestrlist[1]))
        if linestrlist[3] == 'gain':
            truth_type.append("gain")
        else:
            truth_type.append("loss")


count = 0
for i in range(len(result_type)):
    for j in range(len(truth_type)):
        if truth_start[j] <= result_start[i] <= truth_end[j] and truth_type[j] == result_type[i]:
            if result_end[i] <= truth_end[j]:
                ch = (result_end[i] - result_start[i] + 1)
                tr = (truth_end[j]-truth_start[j] +1)
                if ch / tr > 0.5:
                    count += 1
                #count +=(result_end[i] - result_start[i] + 1)
                #print(count)
            elif result_end[i] >= truth_end[j]:
                ch = (truth_end[j] - result_start[i] + 1)
                tr = (truth_end[j] - truth_start[j] + 1)
                if ch / tr > 0.5:
                    count += 1
                #count += (truth_end[j] - result_start[i] + 1)
                #print(count)
            break
        elif truth_start[j] >= result_start[i] and truth_type[j] == result_type[i]:
            if truth_start[j] <= result_end[i] <= truth_end[j]:
                ch = (result_end[i] - truth_start[j] + 1)
                tr = (truth_end[j] - truth_start[j] + 1)
                if ch / tr > 0.5:
                    count += 1
                #count += (result_end[i] - truth_start[j] + 1)
                #print(count)
            elif result_end[i] >= truth_end[j]:
                ch = (truth_end[j] - truth_start[j] + 1)
                tr = (truth_end[j] - truth_start[j] + 1)
                if ch / tr > 0.5:
                    count += 1
                #count += (truth_end[j] - truth_start[j] + 1)
                #print(count)
            break

result_count = len(result_start)
#for i in range(len(result_start)):
#    result_count += (result_end[i] - result_start[i] + 1)

truth_count = len(truth_start)
#for i in range(len(truth_start)):
#    truth_count += (truth_end[i] - truth_start[i] + 1)

#print('count:',count, result_count, truth_count)

output = open("dp_score.txt", "a")
output.write(str(count/result_count) + '\t' + str(count/truth_count) + '\n')
#print(str(count/result_count) + '\t' + str(count/truth_count) + '\n')
