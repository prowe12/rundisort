
# times3.so, which will be imported in python code
times3.so: times3.f
	poetry run f2py -c -m times3 times3.f --fcompiler=gnu95
