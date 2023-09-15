qsea_N2_finiteT: qsea_N2_finiteT.cpp
	g++ -fopenmp -o qsea_N2_finiteT qsea_N2_finiteT.cpp ../ProgressBar/progress.cpp

.PHONY: test clean

test: qsea_N2_finiteT
	./qsea_N2_finiteT

clean:
	rm -f qsea_N2_finiteT
