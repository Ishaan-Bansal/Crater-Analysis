from sympy import *
import numpy as np
import pandas as pd

x = Symbol("x")
arr = np.zeros((6,7))
count1 = 0
count2 = 0

for r in range(20, 140, 20):
    count2 = 0
    for d in range(5, 40, 5):
        f = r**2 - (r**2/d**2) * x**2
        v1 = integrate(pi*f, (x, 0, d)).evalf()
        arr[count1, count2] = v1
        count2 += 1
    count1 += 1

DF = pd.DataFrame(arr)
DF.to_csv("data.csv")