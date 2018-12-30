def DnaToInt (dna) :
    if dna == "A" :
        return 0
    elif dna == "C" :
        return 1
    elif dna == "G" :
        return 2
    else :
        return 3


transition = [0.999, 0.001]

f_test = open("test.fasta", 'r')
f_model = open("model.txt", 'r')

#read test.fatsq
f_test.readline()
test_seq = f_test.readline().rstrip()
f_test.close()

#read model.txt

# 3-d matrix
emission_genic = [[[0]*4 for i in range(4)] for j in range(4)]
emission_intergenic = [[[0]*4 for i in range(4)] for j in range(4)]
emission_start = [[[0]*4 for i in range(4)] for j in range(4)]
emission_stop = [[[0]*4 for i in range(4)] for j in range(4)]

f_model.readline()
for i in range(4):
    for j in range(4):
        tmplst = f_model.readline().rstrip().split(" ")
        for k in range(4):
            emission_genic[i][j][k] = float(tmplst[k])
print("genic")
for i in range(4):
    for j in range(4):
        print(emission_genic[i][j])


f_model.readline()
for i in range(4):
    for j in range(4):
        tmplst = f_model.readline().rstrip().split(" ")
        for k in range(4):
            emission_intergenic[i][j][k] = float(tmplst[k])

print("intergenic")
for i in range(4):
    for j in range(4):
        print(emission_intergenic[i][j])

f_model.readline()
for i in range(4):
    for j in range(4):
        tmplst = f_model.readline().rstrip().split(" ")
        for k in range(4):
            emission_start[i][j][k] = float(tmplst[k])

print("start")
for i in range(4):
    for j in range(4):
        print(emission_start[i][j])


f_model.readline()
for i in range(4):
    for j in range(4):
        tmplst = f_model.readline().rstrip().split(" ")
        for k in range(4):
            emission_stop[i][j][k] = float(tmplst[k])

print("stop")
for i in range(4):
    for j in range(4):
        print(emission_stop[i][j])



f_model.close()


#viterbi algorithm
# first is gene, second is intergene
viterbi = [[None]*2 for i in range(len(test_seq))]
backtrack = [None] * len(test_seq)
viterbi[1][0] = 0.5
viterbi[1][1] = 0.5
backtrack[1] = True
for i in range(2, len(test_seq)) :
    #genic
    viterbi[i][0] = emission_genic[DnaToInt(test_seq[i - 2])][DnaToInt(test_seq[i - 1])][DnaToInt(test_seq[i])] * \
                    (viterbi[i - 1][0] * transition[0] + viterbi[i - 1][1] * transition[1])
    #intergenic
    viterbi[i][1] = emission_intergenic[DnaToInt(test_seq[i - 2])][DnaToInt(test_seq[i - 1])][DnaToInt(test_seq[i])] * \
                   (viterbi[i - 1][1] * transition[0] + viterbi[i - 1][0] * transition[1])

    sum = viterbi[i][0] + viterbi[i][1]
    for j in range(2):
        viterbi[i][j] /= sum



    if viterbi[i][0] < viterbi[i][1] :
        backtrack[i] =  True #genic
    else :
        backtrack[i] = False #not genic


# write result.txt
f_result = open("result.txt", "w")

startIdx = 0
num = 1
for i in range(3, len(test_seq)):
    if backtrack[i] and not backtrack[i-1] :
        startIdx = i

    elif not backtrack[i] and backtrack[i-1] :
        endIdx = i
        data = "genic_" + str(num) +" " + str(startIdx) + " " + str(endIdx) + "\n"
        f_result.write(data)
        num += 1

f_result.close()
