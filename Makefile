.PHONY : build
build : 
	mkdir -p build && cd build && cmake .. && make -j12

.PHONY : run
run :
	cd bin &&./main

.PHONY : clean
clean : 
	rm -rf build