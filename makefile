target := wmsolve

ifndef HYWALL
$(error "Cannot find environment variable HYWALL")
endif

ifndef PTL
$(error "Cannot find environment variable PTL")
endif

main: setup
	mpicxx -I${HYWALL}/include -I${PTL}/include main.cc -o ${target} -L${HYWALL}/lib -lHyWall -L${PTL}/lib -lPTL
	cp -t bin ${target}

setup:
	mkdir -p bin

run: main
	./${target}

clean:
	rm -f ${target}
	rm -rf bin