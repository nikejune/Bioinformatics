import numpy as np

def DNAtoInt(str) :
    lst = []
    for i in str :
        if i == 'A' or i == 'a' :
            lst.append(0)
        elif i == 'C' or i == 'c':
            lst.append(1)
        elif i == 'G' or i == 'g':
            lst.append(2)
        elif i == 'T' or i == 't':
            lst.append(3)

    return lst



fp = open("test01.txt", 'r')

k = int(fp.readline().rstrip())
t = int(fp.readline().rstrip())
m = int(fp.readline().rstrip())
lst = []
lenofstringLst = []
max_pr = 0
answer_motif_lst = []
answer_pr_lst =[]
while True:
    line = fp.readline().rstrip()
    if not line: break
    lst.append(line)
    lenofstringLst.append(len(line))

for times in range(30) :
    print(str(times+1) + "번째 iter")
    #randomly initialize p
    p = np.random.randint(low = 1, high=100 ,size = (4, k))
    #normalize p
    col_sums = p.sum(axis=0)
    p = p / col_sums[ np.newaxis, :]

    #define z
    z =  np.zeros((t, m - k + 1), dtype=float)

    motif_complete_lst = []

    for iter in range(100):
        # set z
        for i in range(t) :
            for j in range(lenofstringLst[i] - k + 1) :
                pow = 1

                tmp = lst[i][j:j+k]

                tmp_int = DNAtoInt(tmp)
                for cur in range(k) :
                    pow *= p[tmp_int[cur]][cur]
                z[i][j] = pow

        #normalize z
        col_sums = z.sum(axis=1)
        z = z / col_sums[ :,  np.newaxis]

        #set p
        new_p = np.zeros((4,k), dtype=float)
        for i in range(t) :
            for j in range(lenofstringLst[i] - k + 1) :
                tmp = lst[i][j:j + k]
                tmp_int = DNAtoInt(tmp)
                for cur in range(k):
                    new_p[tmp_int[cur]][cur] += z[i][j]

        p = new_p
        for i in range(4) :
            for j in range(k) :
                p[i][j] = (p[i][j] + 1)/ (t+4)

        tmpstr = ""
        for i in range(k) :
            sum = 0
            max_j = 0
            max = 0
            for j in range(4):
                 if max < p[j][i] :
                    max = p[j][i]
                    max_j = j
            if max_j == 0:
                tmpstr += "A"
            elif max_j == 1:
                tmpstr += "C"
            elif max_j == 2:
                tmpstr += "G"
            elif max_j == 3:
                tmpstr += "T"
        motif_complete_lst.append(tmpstr)
        print(tmpstr)
        if len(motif_complete_lst) >= 2:
            if motif_complete_lst[-2] == motif_complete_lst[-1] :
                print("finish finding the motif")
                break

    motif_candi_lst = []

    for i in range(t) :
        tmpmotif = []
        min_score = 100000000
        for j in range(m - k + 1) :
            tmp = lst[i][j:j+k]
            score = 0
            for cur in range(k) :
                if tmp[cur] != motif_complete_lst[-1][cur] :
                    score +=1
            if min_score > score :
                min_score = score
                tmpmotif = tmp
        motif_candi_lst.append(tmpmotif)

    profile = np.zeros((4, k), dtype=int)
    pr_lst = []
    for i in motif_candi_lst :
        for j_idx , j in enumerate(i) :
            if j == 'A' or j == 'a':
                profile[0][j_idx] += 1
            elif j == 'C' or j == 'c':
                profile[1][j_idx] += 1
            elif j == 'G' or j == 'g':
                profile[2][j_idx] += 1
            elif j == 'T' or j == 't':
                profile[3][j_idx] += 1

    col_sums = profile.sum(axis=0)
    profile = profile / col_sums[ np.newaxis, :]

    for i in motif_candi_lst :
        mul = 1
        for j_idx , j in enumerate(i) :
            if j == 'A' or j == 'a':
                mul *= profile[0][j_idx]
            elif j == 'C' or j == 'c':
                mul *= profile[1][j_idx]
            elif j == 'G' or j == 'g':
                mul *= profile[2][j_idx]
            elif j == 'T' or j == 't':
                mul *= profile[3][j_idx]
        pr_lst.append(mul)

    sum_pr = 0
    for i in pr_lst :
        sum_pr += i
    if max_pr < sum_pr :
        max_pr = sum_pr
        answer_motif_lst = motif_candi_lst
        answer_pr_lst = pr_lst
        if max_pr > 1 :
            break

    for i in motif_candi_lst:
        print(i)

    for i in pr_lst:
        print(i)

    print("sumpr :", end='')
    print(sum_pr)


print("answer")
for i in answer_motif_lst :
    print(i)

for i in answer_pr_lst :
    print(i)

fp.close()