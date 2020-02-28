import os

os.system('echo "====================== SERIAL TIME OF EXEC ======================"')

os.environ["input"] = "blocking"
os.environ["num_proc"] = "2"

for i in range(5,13):
	num = 2**i;
	os.environ["N"] = str(num)
	os.system('make compile')
	os.system('make run')


os.system('echo "================ BLOCKING P2P TIME OF EXEC (2 processes) ================"')

os.environ["input"] = "blocking"
os.environ["num_proc"] = "3"

for i in range(5,13):
	num = 2**i;
	os.environ["N"] = str(num)
	os.system('make compile')
	os.system('make run')


os.system('echo "================ NON-BLOCKING P2P TIME OF EXEC (2 processes) ================"')

os.environ["input"] = "noblocking"
os.environ["num_proc"] = "3"

for i in range(5,13):
	num = 2**i;
	os.environ["N"] = str(num)
	os.system('make compile')
	os.system('make run')


os.system('echo "=================== COLLEVTIVE MPI TIME OF EXEC (2 processes) ==================="')

os.environ["input"] = "collective"
os.environ["num_proc"] = "2"

for i in range(5,13):
	num = 2**i;
	os.environ["N"] = str(num)
	os.system('make compile')
	os.system('make run')


os.system('echo "=================== COLLEVTIVE MPI TIME OF EXEC (4 processes) ==================="')

os.environ["input"] = "collective"
os.environ["num_proc"] = "4"

for i in range(5,13):
	num = 2**i;
	os.environ["N"] = str(num)
	os.system('make compile')
	os.system('make run')
