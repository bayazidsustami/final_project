FILE 		= COBA1

INCLUDE 	= 	-I./ \
				-I"/usr/include/opencv2" \
				-I"/usr/local/include"
				

LIB 		=		-L"/usr/include/opencv2" -lopencv_core -lopencv_highgui -lopencv_imgproc \
				-L"/usr/local/lib" -lfftw3f03
				
.PHONY: clean

%: %.cpp
	g++ -g $(INCLUDE) $(LIB) -o $@ $<

txt:$(FILE).exe
	./$(FILE)>output.txt
clean: 
	rm -f *.exe
