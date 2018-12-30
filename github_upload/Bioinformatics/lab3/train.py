
def DnaToInt (dna) :
    if dna == "A" :
        return 0
    elif dna == "C" :
        return 1
    elif dna == "G" :
        return 2
    else :
        return 3

# read files
f_genic = open("train_genic.fasta", 'r')
f_intergenic = open("train_intergenic.fasta", 'r')

# 3-d matrix
emission_genic = [[[0]*4 for i in range(4)] for j in range(4)]
emission_intergenic = [[[0]*4 for i in range(4)] for j in range(4)]
emission_start = [[[0]*4 for i in range(4)] for j in range(4)]
emission_stop = [[[0]*4 for i in range(4)] for j in range(4)]

# genic
while True:
    line = f_genic.readline()
    if not line: break
    line = f_genic.readline().rstrip()
    tmp = []
    for i in range(len(line)):
        tmp.append(DnaToInt(line[i]))
    for i in range(1,len(tmp)- 3) :
        emission_genic[tmp[i]][tmp[i+1]][tmp[i+2]]+= 1

    emission_start[tmp[0]][tmp[1]][tmp[2]] += 1
    emission_stop[tmp[-3]][tmp[-2]][tmp[-1]] += 1

f_genic.close()

# genic normalizing
for i in range(4):
    for j in range(4):
        sum = 0
        for k in range(4):
            sum += emission_genic[i][j][k]
        for k in range(4):
            emission_genic[i][j][k] = round(emission_genic[i][j][k]/sum, 4)

print("genic")
for i in range(4):
    for j in range(4):
        print(emission_genic[i][j])

# intergenic

while True:
    line = f_intergenic.readline()
    if not line: break
    line = f_intergenic.readline().rstrip()
    tmp = []
    for i in range(len(line)):
        tmp.append(DnaToInt(line[i]))
    for i in range(len(tmp) - 2):
        emission_intergenic[tmp[i]][tmp[i + 1]][tmp[i + 2]] += 1

f_intergenic.close()

# intergenic normalizing
for i in range(4):
    for j in range(4):
        sum = 0
        for k in range(4):
            sum += emission_intergenic[i][j][k]
        for k in range(4):
            emission_intergenic[i][j][k] = round(emission_intergenic[i][j][k] / sum, 4)

print("intergenic")
for i in range(4):
    for j in range(4):
        print(emission_intergenic[i][j])

# start normalizing
for i in range(4):
    for j in range(4):
        sum = 0
        for k in range(4):
            sum += emission_start[i][j][k]
        for k in range(4):
            try :
                emission_start[i][j][k] = round(emission_start[i][j][k] / sum, 4)
            except :
                emission_start[i][j][k] = 0
print("start")
for i in range(4):
    for j in range(4):
        print(emission_start[i][j])

# stop normalizing
for i in range(4):
    for j in range(4):
        sum = 0
        for k in range(4):
            sum += emission_stop[i][j][k]
        for k in range(4):
            try :
                emission_stop[i][j][k] = round(emission_stop[i][j][k] / sum, 4)
            except :
                emission_stop[i][j][k] = 0

print("stop")
for i in range(4):
    for j in range(4):
        print(emission_stop[i][j])

#write model.txt
f_model = open("model.txt", 'w')
#genic
f_model.write("emission_gene=\n")
for i in range(4):
    for j in range(4):
        data = str(emission_genic[i][j][0]) +" "+ str(emission_genic[i][j][1]) +" "+ \
        str(emission_genic[i][j][2]) + " " +str(emission_genic[i][j][3]) +"\n"
        f_model.write(data)
#intergenic
f_model.write("emission_intergene=\n")
for i in range(4):
    for j in range(4):
        data = str(emission_intergenic[i][j][0]) +" "+ str(emission_intergenic[i][j][1]) +" "+ \
        str(emission_intergenic[i][j][2]) + " " +str(emission_intergenic[i][j][3]) +"\n"
        f_model.write(data)

#start
f_model.write("emission_start=\n")
for i in range(4):
    for j in range(4):
        data = str(emission_start[i][j][0]) +" "+ str(emission_start[i][j][1]) +" "+ \
        str(emission_start[i][j][2]) + " " +str(emission_start[i][j][3]) +"\n"
        f_model.write(data)

#stop
f_model.write("emission_stop=\n")
for i in range(4):
    for j in range(4):
        data = str(emission_stop[i][j][0]) +" "+ str(emission_stop[i][j][1]) +" "+ \
        str(emission_stop[i][j][2]) + " " +str(emission_stop[i][j][3]) +"\n"
        f_model.write(data)
f_model.close()