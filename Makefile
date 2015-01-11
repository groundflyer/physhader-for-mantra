HOU_VERSION=`echo $(HOUDINI_MAJOR_RELEASE).$(HOUDINI_MINOR_RELEASE)`
HOU_DIR="$(HOME)/houdini$(HOU_VERSION)"
DEST_OTLS="$(HOU_DIR)/otls"

all:
	hotl -c expanded-otl physhader.otl

install:
	cp physhader.otl $(DEST_OTLS)
	cp -r vex $(HOU_DIR)

clean:
	rm physhader.otl
