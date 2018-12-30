
def outputLCS (backtrack, str1, str2, len_str1 , len_str2) :
    stk = real_outputLCS(backtrack, len_str1, len_str2)

    numOfMatches = 0
    numOfMismatches = 0
    numOfGaps = 0
    result_str1 = str1[stk[0]]
    result_str2 = str2[stk[0]]

    if str1[stk[0]] != str2[stk[0]] :
        numOfMismatches += 1
    else :
        numOfMatches += 1
    for idx in range(1, int(len(stk)/2)) :
        i = 2 * idx
        j = 2 * idx + 1
        # gap of str1
        if stk[i] == stk[i-2] :
            result_str1 += "-"
            result_str2 += str2[stk[j]]
            numOfGaps += 1
        # gap of str2
        elif stk[j] == stk[j-2] :
            result_str1 += str1[stk[i]]
            result_str2 += "-"
            numOfGaps += 1
        # mismatch
        elif str1[stk[i]] != str2[stk[j]] :
            result_str1 += str1[stk[i]]
            result_str2 += str2[stk[j]]
            numOfMismatches += 1
        else :
            result_str1 += str1[stk[i]]
            result_str2 += str2[stk[j]]
            numOfMatches += 1

    print(result_str1)
    print(result_str2)
    return list([numOfMatches, numOfMismatches, numOfGaps])


'''
str1의 인덱스 i와 str2의 인덱스 j를 리스트로 반환하며 계속 붙여나간다.
'''
def real_outputLCS (backtrack, i , j) :
    if i == 0 or j == 0 :
        return list()
    if backtrack[i][j] == "down":
        tmp = real_outputLCS(backtrack, i - 1, j)
        return tmp + list([i,j])
    elif backtrack[i][j] == "right":
        tmp = real_outputLCS(backtrack, i, j - 1)
        return tmp + list([i, j])
    else :
        tmp = real_outputLCS(backtrack, i - 1, j -1)
        return tmp + list([i, j])



# read files
f_seq1 = open("seq1.txt", 'r')
f_seq2 = open("seq2.txt", 'r')
f_score = open("score.txt", 'r')

seq1_name = f_seq1.readline().rstrip()
seq1_str =  " "+f_seq1.readline().rstrip().upper()

seq2_name = f_seq2.readline().rstrip()
seq2_str =  " "+f_seq2.readline().rstrip().upper()

match = int(f_score.readline().split("=")[1].strip())
mismatch =  int(f_score.readline().split("=")[1].strip())
gap = int(f_score.readline().split("=")[1].strip())

# global alignment
# matrix initialization
matrix = []
for i in range(len(seq1_str)) :
    tmp = []
    for j in range(len(seq2_str)) :
        tmp.append(0)
    matrix.append(tmp)

for i in range(len(seq1_str)):
    matrix[i][0] = i * gap

for j in range(len(seq2_str)):
    matrix[0][j] = j * gap

#backtrack initialization
backtrack = []
for i in range(len(seq1_str)) :
    tmp = []
    for j in range(len(seq2_str)) :
        tmp.append("N")
    backtrack.append(tmp)


# LCS backtrack
for i in range(1, len(seq1_str)) :
    for j in range(1, len(seq2_str)) :
        if seq1_str[i] == seq2_str[j] :
            target = match
        else :
            target = mismatch

        matrix[i][j] = max([matrix[i-1][j] + gap , matrix[i][j-1] + gap, matrix[i-1][j-1] + target])
        if matrix[i][j] == matrix[i - 1][j] + gap:
            backtrack[i][j] = "down"
        elif matrix[i][j] == matrix[i][j - 1] + gap:
            backtrack[i][j] = "right"
        elif matrix[i][j] == matrix[i - 1][j - 1] + target:
            backtrack[i][j] = "diagonal"




print(seq1_name, seq1_str)
print(seq2_name, seq2_str)
print("============GLOBAL ALIGNMENT===========")
numOfMatches, numOfMismatches, numOfGaps = outputLCS(backtrack, seq1_str, seq2_str, len(seq1_str)-1, len(seq2_str)-1)
print("matches :", numOfMatches)
print("mismatches :", numOfMismatches)
print("gaps :", numOfGaps)
print("scores :", matrix[-1][-1])



# local alignment
# matrix initialization
matrix = []
for i in range(len(seq1_str)) :
    tmp = []
    for j in range(len(seq2_str)) :
        tmp.append(0)
    matrix.append(tmp)

#for i in range(len(seq1_str)):
#    matrix[i][0] = i * gap

#for j in range(len(seq2_str)):
#    matrix[0][j] = j * gap

#backtrack initialization
backtrack = []
for i in range(len(seq1_str)) :
    tmp = []
    for j in range(len(seq2_str)) :
        tmp.append("N")
    backtrack.append(tmp)

# LCS Backtrack
for i in range(1, len(seq1_str)) :
    for j in range(1, len(seq2_str)) :
        if seq1_str[i] == seq2_str[j]:
            target = match
        else:
            target = mismatch

        matrix[i][j] = max([0, matrix[i - 1][j] + gap, matrix[i][j - 1] + gap, matrix[i - 1][j - 1] + target])
        if matrix[i][j] == matrix[i - 1][j] + gap:
            backtrack[i][j] = "down"
        elif matrix[i][j] == matrix[i][j - 1] + gap:
            backtrack[i][j] = "right"
        elif matrix[i][j] == matrix[i - 1][j - 1] + target:
            backtrack[i][j] = "diagonal"


print("======LOCAL ALIGNMENT(not affine)======")
numOfMatches, numOfMismatches, numOfGaps = outputLCS(backtrack, seq1_str, seq2_str, len(seq1_str)-1, len(seq2_str)-1)
print("matches :", numOfMatches)
print("mismatches :", numOfMismatches)
print("gaps :", numOfGaps)
print("scores :", matrix[-1][-1])



# local alignment(affine)
# matrix initialization (0 : lower, 1: middle, 2:upper)

print('''
*****************Enter the penalties***************
(extension penalty is less than the open penalty)
example : open_penalty : 3, extension_penalty : 1
''')
open_penalty = int(input("open_penalty : "))
extension_penalty = int(input("extension_penalty : "))

matrix = []
for i in range(len(seq1_str)):
    for j in range(len(seq2_str)) :
        tmp = []
        for k in range(len(seq2_str)) :
            tmp.append(list([0,0,0]))
        matrix.append(tmp)


#for i in range(len(seq1_str)):
#    matrix[i][0] = i * gap

#for j in range(len(seq2_str)):
#    matrix[0][j] = j * gap

#backtrack initialization
backtrack = []
for i in range(len(seq1_str)) :
    tmp = []
    for j in range(len(seq2_str)) :
        tmp.append("N")
    backtrack.append(tmp)

# LCS Backtrack (diagonal)
for idx_i in range(1, len(seq1_str) + len(seq2_str)) :
    for idx_j in range(len(seq2_str)-1) :
        i = idx_i - idx_j
        j = idx_j + 1
        if i < 1:
            break
        if i >= len(seq1_str) :
            continue

        if seq1_str[i] == seq2_str[j]:
            target = match
        else:
            target = mismatch

        matrix[i][j][0] = max([matrix[i - 1][j][0] - extension_penalty, matrix[i - 1][j][1] - open_penalty])
        matrix[i][j][2] = max([matrix[i][j - 1][2] - extension_penalty, matrix[i][j - 1][1] - open_penalty])

        matrix[i][j][1] = max([matrix[i][j][0], matrix[i][j][2], matrix[i - 1][j - 1][1] + target])
        if matrix[i][j][1] == matrix[i][j][0]:
            backtrack[i][j] = "down"
        elif matrix[i][j][1] == matrix[i][j][2]:
            backtrack[i][j] = "right"
        elif matrix[i][j][1] == matrix[i - 1][j - 1][1] + target:
            backtrack[i][j] = "diagonal"



print("========LOCAL ALIGNMENT(affine)========")
numOfMatches, numOfMismatches, numOfGaps = outputLCS(backtrack, seq1_str, seq2_str, len(seq1_str)-1, len(seq2_str)-1)
print("matches :", numOfMatches)
print("mismatches :", numOfMismatches)
print("gaps :", numOfGaps)
print("scores :", matrix[len(seq1_str)-1][len(seq2_str)-1][1])

f_seq1.close()
f_seq2.close()
f_score.close()


'''
print(" ",end="")
for j in range(0, len(seq2_str)) :
    print("%5s" % seq2_str[j],end ="")
print("")

for i in range(len(seq1_str)) :
    print(seq1_str[i], end="")
    for j in range(len(seq2_str)) :
        print("%5d" % matrix[i][j], end="")
    print(" ")

for i in range(len(seq1_str)):
    for j in range(len(seq2_str)):
        print("%10s" %(backtrack[i][j]) ,end="")
    print(" ")
'''