
# Create the library disort_driver.so, which will be imported in python code

disort_driver_py.so: disort_driver_py.f disort_py.f errpack.f linpak.f rdi1mach.f bdref.f
	f2py -c -m disort_driver_py disort_driver_py.f disort_py.f errpack.f linpak.f rdi1mach.f bdref.f --fcompiler=gnu95

# getmom.f 
# getmom.f 

# This is the legendre moments fortran library
#moms.so: getmompf.f getmom.f
#	f2py -c -m moms getmompf.f getmom.f  --fcompiler=gnu95

