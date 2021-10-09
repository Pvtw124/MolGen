from numba import jit, cuda

@jit
def test():
	print("hello")
	
if __name__ == "__main__":
	test()
