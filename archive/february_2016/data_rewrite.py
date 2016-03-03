import os
print os.getcwd()

os.makedirs("data_3ns")

for i in range(1,102):
    with open("data/" + str(i)) as f:
		with open("data_3ns/" + str(i)) as out:
			for line in f[0:30000]:
				out.write(line)
			
			