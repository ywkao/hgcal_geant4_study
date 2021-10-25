ROOTLIBS = $(shell root-config --cflags --libs) -lMinuit -lMLP -lXMLIO -lTMVA -lGenVector

CC = g++
CFLAGS = -g -std=c++1y

TARGET = examine

all: $(TARGET)

examine: examine.cpp hits.C hits.h
	$(CC) $(CFLAGS) -I${ROOTSYS}/include -o examine examine.cpp $(ROOTLIBS)
