CC=g++
CFLAGS=-c -fPIC `root-config --cflags`
LDFLAGS=`root-config --libs` -lRooFit -lRooFitCore
SOURCES= SimpleKalmanCalculator.cc KalmanMuonCalibrator.cc calibrate.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=calibrate

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS)

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm ./*~ ./*.o ./calibrate


