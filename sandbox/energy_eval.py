import numpy as np
from numpy import linalg as LA

h2oSTO3g_F=[
[-32.255 ,  -2.7915  ,       0  ,0.008611   ,      0  , -0.1813  , -0.1813],
[-2.7915 ,  -8.2369  ,       0  ,-0.22829   ,      0  , -0.3858  , -0.3858],
[0       ,  0        , -7.5429  ,       0   ,      0  ,-0.11321  , 0.11321],
[0.008611,  -0.22829 ,        0 ,   -7.457  ,       0 , -0.11022 , -0.11022],
[0       ,  0        , 0        , 0         , -7.3471 ,        0 ,        0],
[-0.1813 ,  -0.3858  ,-0.11321  ,-0.11022   ,      0  ,  -4.033  ,  -0.044647],
[-0.1813 ,  -0.3858  , 0.11321  ,-0.11022   ,      0  , -0.044647,    -4.033]
]

density_mat=[
[   1.0650117 , -0.2852166,  -0.0000000,  -0.0195534,  -0.0000000,   0.0334496,   0.0334496],
[  -0.2852166 ,  1.2489657,   0.0000000,   0.1135594,   0.0000000,  -0.1442809,  -0.1442809],
[  -0.0000000 ,  0.0000000,   1.1258701,  -0.0000000,  -0.0000000,  -0.1461317,   0.1461317],
[  -0.0195534 ,  0.1135594,  -0.0000000,   1.0660638,   0.0000000,  -0.0993583,  -0.0993583],
[  -0.0000000 ,  0.0000000,  -0.0000000,   0.0000000,   1.0000000,  -0.0000000,  -0.0000000],
[   0.0334496 , -0.1442809,  -0.1461317,  -0.0993583,  -0.0000000,   0.0426802,   0.0047460],
[   0.0334496 , -0.1442809,   0.1461317,  -0.0993583,  -0.0000000,   0.0047460,   0.0426802]
]

core_ham=[
[ -32.5773954,  -7.5788328,   0.0000000,  -0.0144738,   0.0000000,  -1.2401023,  -1.2401023],
[  -7.5788328,  -9.2009433,   0.0000000,  -0.1768902,   0.0000000,  -2.9067098,  -2.9067098],
[   0.0000000,   0.0000000,  -7.4588193,   0.0000000,   0.0000000,  -1.6751501,   1.6751501],
[  -0.0144738,  -0.1768902,   0.0000000,  -7.4153118,   0.0000000,  -1.3568683,  -1.3568683],
[   0.0000000,   0.0000000,   0.0000000,   0.0000000,  -7.3471449,   0.0000000,   0.0000000],
[  -1.2401023,  -2.9067098,  -1.6751501,  -1.3568683,   0.0000000,  -4.5401711,  -1.0711459],
[  -1.2401023,  -2.9067098,   1.6751501,  -1.3568683,   0.0000000,  -1.0711459,  -4.5401711]]

result = np.zeros_like(core_ham)
for i in range(len(core_ham)):
    for j in range(len(core_ham[i])):
        result[i][j]+=core_ham[i][j]+h2oSTO3g_F[i][j]
print(result)
finmat=np.matmul(density_mat,result)
print(finmat)
print(np.trace(finmat))