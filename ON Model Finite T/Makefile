qsea_ON_finiteT: qsea_ON_finiteT.cpp
	g++ -fopenmp -o qsea_ON_finiteT qsea_ON_finiteT.cpp ../ProgressBar/progress.cpp

.PHONY: test clean

test: qsea_ON_finiteT
	./qsea_ON_finiteT

clean:
	rm -f qsea_ON_finiteT
