all: CarveItPosix CarveItPosixTimed CarveItPosixTimedOpt
tests: test testTimed testTimedOpt testTimedOpt1

CarveItPosix: CarveItPosixThreadsTimed.c
	gcc -o CarveItThreaded -O2 CarveItPosixThreadsTimed.c -lpthread -lm

CarveItPosixTimed: CarveItPosixThreadsTimed.c clock.c 
	gcc -o CarveItThreadedTimed -DTIMED CarveItPosixThreadsTimed.c clock.c -lpthread -lm

CarveItPosixTimedOpt: CarveItPosixThreadsTimed.c clock.c 
	gcc -o CarveItThreadedTimedOpt -g -O2 -DTIMED CarveItPosixThreadsTimed.c clock.c -lpthread -lm

test: carveInput.dat carveOutput.dat
	.\CarveItThreaded carveInput.dat carvePosixOutput.dat 4
	diff carveOutput.dat carvePosixOutput.dat

testTimed: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimed carveInput.dat carvePosixOutputTimed.dat 4 > timeData.txt
	diff carveOutput.dat carvePosixOutputTimed.dat

testTimedOpt: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt.dat 4 > timeDataOpt.txt
	diff carveOutput.dat carvePosixOutputTimedOpt.dat

testTimedOpt1: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt1.dat 1 > timeDataOpt1.txt
	diff carveOutput.dat carvePosixOutputTimedOpt1.dat

clean:
	erase carvePosixOutput.dat carvePosixOutputTimed.dat timeData.txt
	erase *.exe
	erase carvePosixOutputTimedOpt.dat timeDataOpt.txt
	erase carvePosixOutputTimedOpt1.dat timeDataOpt1.txt

