all: CarveItPosix CarveItPosixTimed CarveItPosixTimedOpt CarveItPosixTimedOptAll CarveAllMod CarveAllMod2

tests: test testTimed testTimedOpt16 testTimedOpt10 testTimedOpt12 testTimedOpt10 testTimedOpt8 testTimedOpt6 testTimedOpt4 testTimedOpt2 testTimedOpt1
opttests: testTimedOpt1 testTimedOpt2 testTimedOpt4 testTimedOpt6 testTimedOpt8 testTimedOpt10 testTimedOpt12 testTimedOpt14 testTimedOpt16 testTimedOpt18 testTimedOpt20 testTimedOpt36

GCC = C:\TDM-GCC-64\bin\gcc.exe

CarveItPosix: CarveItPosixThreadsTimed.c
	$(GCC) -o CarveItThreaded -O2 CarveItPosixThreadsTimed.c -lpthread -lm

CarveItPosixTimed: CarveItPosixThreadsTimed.c clock.c 
	$(GCC) -o CarveItThreadedTimed -DTIMED CarveItPosixThreadsTimed.c clock.c

CarveItPosixTimedOpt: CarveItPosixThreadsTimed.c clock.c 
	$(GCC) -o CarveItThreadedTimedOpt -g -O2 -DTIMED CarveItPosixThreadsTimed.c clock.c

CarveItPosixTimedOptAll: CarveItPosixThreadsTimedForAll.c clock.c 
	$(GCC) -o CarveItThreadedTimedForAll -O3 -DTIMED CarveItPosixThreadsTimedForAll.c clock.c

CarveAllMod: CarveAllMod.c clock.c
	$(GCC) -o CarveAllMod -O3 -DTIMED CarveAllMod.c clock.c

CarveAllMod2: CarveAllMod2.c clock.c
	$(GCC) -o CarveAllMod2 -O3 -DTIMED CarveAllMod2.c clock.c

test: carveInput.dat carveOutput.dat
	.\CarveItThreaded carveInput.dat carvePosixOutput.dat 4 ; diff carveOutput.dat carvePosixOutput.dat

testTimed: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimed carveInput.dat carvePosixOutputTimed.dat 4 > timeData4.txt ; diff carveOutput.dat carvePosixOutputTimed.dat

testTimedOpt36: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt36.dat 36 > timeDataOpt36.txt ; diff carveOutput.dat carvePosixOutputTimedOpt36.dat

testTimedOpt20: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt20.dat 20 > timeDataOpt20.txt ; diff carveOutput.dat carvePosixOutputTimedOpt20.dat

testTimedOpt18: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt18.dat 18 > timeDataOpt18.txt ; diff carveOutput.dat carvePosixOutputTimedOpt18.dat

testTimedOpt16: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt16.dat 16 > timeDataOpt16.txt ; diff carveOutput.dat carvePosixOutputTimedOpt16.dat

testTimedOpt14: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt14.dat 14 > timeDataOpt14.txt ; diff carveOutput.dat carvePosixOutputTimedOpt14.dat

testTimedOpt12: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt12.dat 12 > timeDataOpt12.txt ; diff carveOutput.dat carvePosixOutputTimedOpt12.dat

testTimedOpt10: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt10.dat 10 > timeDataOpt10.txt ; diff carveOutput.dat carvePosixOutputTimedOpt10.dat

testTimedOpt8: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt8.dat 8 > timeDataOpt8.txt ; diff carveOutput.dat carvePosixOutputTimedOpt8.dat

testTimedOpt6: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt6.dat 6 > timeDataOpt6.txt ; diff carveOutput.dat carvePosixOutputTimedOpt6.dat

testTimedOpt4: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt4.dat 4 > timeDataOpt4.txt ; diff carveOutput.dat carvePosixOutputTimedOpt4.dat

testTimedOpt2: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt2.dat 2 > timeDataOpt2.txt ; diff carveOutput.dat carvePosixOutputTimedOpt2.dat

testTimedOpt1: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedOpt carveInput.dat carvePosixOutputTimedOpt1.dat 1 > timeDataOpt1.txt ; diff carveOutput.dat carvePosixOutputTimedOpt1.dat

test40: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedForAll carveInput.dat carvePosixOutputTimed.dat 36 40 > out40.txt

test80: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedForAll carveInput.dat carvePosixOutputTimed.dat 36 80 > out80.txt

test100: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedForAll carveInput.dat carvePosixOutputTimed.dat 36 100 > out100ave.txt

test500: carveInput.dat carveOutput.dat
	.\CarveItThreadedTimedForAll carveInput.dat carvePosixOutputTimed.dat 36 500 > out500ave.txt

test500Mod: carveInput.dat carveOutput.dat
	.\CarveAllMod carveInput.dat carveOutputMod.dat 36 500 > out500mod.txt

test5Mod: carveInput.dat carveOutput.dat
	.\CarveAllMod carveInput.dat carveOutputMod.dat 36 5 > out5mod.txt

test500Mod2: carveInput.dat carveOutput.dat
	.\CarveAllMod2 carveInput.dat carveOutputMod.dat 36 500 > out500mod2.txt

test5Mod2: carveInput.dat carveOutput.dat
	.\CarveAllMod2 carveInput.dat carveOutputMod.dat 36 5 > out5mod2.txt

clean:
	erase carvePosixOutput.dat carvePosixOutputTimed.dat timeData.txt
	erase *.exe
	erase carvePosixOutputTimedOpt.dat timeDataOpt.txt
	erase carvePosixOutputTimedOpt1.dat timeDataOpt1.txt

