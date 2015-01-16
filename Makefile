HOU_VERSION=$(HOUDINI_MAJOR_RELEASE).$(HOUDINI_MINOR_RELEASE)
HOU_DIR=$(HOME)/houdini$(HOU_VERSION)
OTL=otls/physhader.otl

all:
	hotl -c expanded-otl $(OTL)

install:
	cp -r otls vex gallery $(HOU_DIR)

clean:
	rm $(OTL)
