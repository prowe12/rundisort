
# test.so, which will be imported in python code
#test.so: test.f
#	poetry run f2py -c -m test test.f --fcompiler=gnu95

# test2.so, which will be imported in python code
test2.so: test2.f
	poetry run f2py -c -m test2 test2.f --fcompiler=gnu95
