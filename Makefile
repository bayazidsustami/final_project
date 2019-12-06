FILE 		= COBA1

INC 	= 	-I./ \
				-I"/usr/include/opencv2" \
				-I"/usr/local/include"
				

LIB 		=		-L"/usr/include/opencv2" -lopencv_core -lopencv_highgui -lopencv_imgproc \
				-L "/usr/local/lib" -l fftw3f
				
.PHONY: clean

%: %.cpp
	g++ -o $@ $< $(INC) $(LIB) 

txt:$(FILE).exe
	./$(FILE)>output.txt
clean: 
	rm -f *.exe
