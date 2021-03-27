NTL = -lntl -lgmp -L/usr/local/lib

test: test.o ECFactor.o EC_p.o
	g++ test.o ECFactor.o EC_p.o $(NTL)