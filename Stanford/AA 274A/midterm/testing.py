import numpy as np

a = np.array([8,2,6,1])
print(a[np.array([1,2,4,2])<3])

b = np.array([[2,4,1],[2,5,1]])
print(" ")
print(b[range(2),:][np.array([[True],[False]])])