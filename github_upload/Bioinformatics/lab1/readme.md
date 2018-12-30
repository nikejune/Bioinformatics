# Find a motif from a set of strings using EM and gibbs sampling

```
pip install numpy
```

numpy 를 설치하여야 합니다.



```
fp = open("test01.txt", 'r')
```

19번째 라인에 test01.txt로 파일 입력을 받습니다. 



EM 알고리즘을 이용하였으며, p와 Z가 서로 학습을 하는 과정에서 p를 통해 알 수 있는consensus가 변하지 않을 때 EM알고리즘을 종료

그 다음 profile값을 저장한뒤 이 과정을 반복하여서 가장 높은 profile값이 나오는 것을 최종 답으로 채택합니다.

