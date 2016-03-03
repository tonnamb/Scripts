import os
print os.getcwd()
from itertools import islice

if not os.path.exists("data_slice"):
	os.makedirs("data_slice")

for i in range(1,51):
	with open('data/' + str(i)) as fin, open('data_slice/' + str(i), 'w') as fout:
		fout.writelines(islice(fin, None, None, 10))