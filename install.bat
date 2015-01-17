@set HOUDINI_MAJOR_RELEASE=14
@set HOUDINI_MINOR_RELEASE=0
@set HOUDINI_BUILD=201
@set HOUDINI_VERSION=%HOUDINI_MAJOR_RELEASE%.%HOUDINI_MINOR_RELEASE%.%HOUDINI_BUILD%.13
@set HOTL="C:\Program Files\Side Effects Software\Houdini %HOUDINI_VERSION%\bin\hotl.exe"
@set OTL=physhader.otl
@set HOUDINI_HOME_FOLDER="%USERPROFILE%\Documents\houdini%HOUDINI_MAJOR_RELEASE%.%HOUDINI_MINOR_RELEASE%"

%HOTL% -c expanded-otl %OTL%

move %OTL% "%HOUDINI_HOME_FOLDER%\otls\"
%systemroot%\System32\xcopy vex "%HOUDINI_HOME_FOLDER%\vex" /s/h/e/k/f/c
%systemroot%\System32\xcopy gallery "%HOUDINI_HOME_FOLDER%\gallery" /s/h/e/k/f/c

@pause
