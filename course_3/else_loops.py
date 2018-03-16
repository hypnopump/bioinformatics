# -*- coding: utf-8 -*-

n = 10 # magic number
for y in range(2, n):
	for x in range(2, y):
		if y%x == 0:
			print(y, "equals", x, "*", y//x)
			break
	else:
		# loop fell through without finding a factor
		print(y, "is a prime number")