HOU_VERSION="$(HOUDINI_MAJOR_RELEASE).$(HOUDINI_MINOR_RELEASE)"
HOU_DIR="$(HOME)/houdini$(HOU_VERSION)"
DEST_OTLS="$(HOU_DIR)/otls"

OTL="physhader.otl"

all:
	hotl -c expanded-otl $(OTL)

install:
	cp $(OTL) $(DEST_OTLS)
	cp -r vex $(HOU_DIR)
	cp -r gallery $(HOU_DIR)

clean:
	rm $(OTL)
