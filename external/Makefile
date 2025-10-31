
BASE64LIB := WinBase64/libbase64.o
INCLUDES += -IWinBase64


.PHONY : all clean 

all: $(BASE64LIB) miniz.o    

miniz.o: miniz.c miniz.h
	$(CXX) -O3 -fPIC -c -o $@ $<

$(BASE64LIB): 
	make -C WinBase64


clean:
	make clean
