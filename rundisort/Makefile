
# Create the library disort_driver.so, which will be imported in python code

disort_driver_py.so: disort_driver_py.f disort_py.f errpack.f linpak.f rdi1mach.f bdref.f
	f2py -c -m --verbose --opt='-std=legacy' disort_driver_py disort_driver_py.f disort_py.f errpack.f linpak.f rdi1mach.f bdref.f --fcompiler=gnu95
