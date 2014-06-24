CPPFLAGS=`root-config --cflags`
LDFLAGS=-L${ROOTSYS}/lib -L${ROOFITSYS}/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lHtml -lMinuit -lRooFitCore -lRooFit -lTMVA -lRooStats -lFoam -lMathCore -lMathMore -g

all: RegVars RegVal AddRegVars

RegVars: src/RegVars.C
	g++ $(CPPFLAGS) $(LDFLAGS) -I/${ROOFITSYS}/include -o $@ $^

RegVal: src/RegVal.C
	g++ $(CPPFLAGS) $(LDFLAGS) -I/${ROOFITSYS}/include -o $@ $^

AddRegVars: src/AddRegVars.C
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

