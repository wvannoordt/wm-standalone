ifndef HYWALL
$(error "Cannot find environment variable HYWALL")
endif

ifndef PTL
$(error "Cannot find environment variable PTL")
endif

main:
	mpicxx -I${HYWALL}/include -I${PTL}/include main.cc -o program -L${HYWALL}/lib -lHyWall -L${PTL}/lib -lPTL

run: main
	./program

clean:
	rm -f program
