if "%BUILD_VS_YEAR%"=="2019" set BUILD_VS_LIBDIRPOST=vc15
if "%BUILD_VS_YEAR%"=="2022" set BUILD_VS_LIBDIRPOST=vc15

set BUILD_VS_SVNDIR=win64_%BUILD_VS_LIBDIRPOST%
set BUILD_VS_LIBDIR=lib/windows_x64

if NOT "%verbose%" == "" (
	echo Library Directory = "%BUILD_VS_LIBDIR%"
)
if NOT EXIST %BUILD_VS_LIBDIR% (
	rem libs not found, but svn is on the system
	if not "%GIT%"=="" (
		echo.
		echo The required external libraries in %BUILD_VS_LIBDIR% are missing
				echo.
		echo Downloading %BUILD_VS_LIBDIR% libraries, please wait.
		echo.
		echo *********************************************************
		echo *                                                       *
		echo * Note: Once the initial download finishes and you see  *
		echo *       "Resolving deltas: 100%% (nnn/nnn) done"         *
		echo *       a second, much larger, update will occur with   *
		echo *       no visible updates. Please do not interrupt     *
		echo *       this process. It may take over an hour to       *
		echo *       complete depending on your internet connection. *
		echo *                                                       *
		echo *********************************************************
		"%GIT%" -C "%BLENDER_DIR%\" config --local "submodule.%BUILD_VS_LIBDIR%.update" "checkout"
		set GIT_LFS_SKIP_SMUDGE=1
		"%GIT%" -C "%BLENDER_DIR%\" submodule update --progress --init "%BUILD_VS_LIBDIR%"
		set GIT_LFS_SKIP_SMUDGE=
		"%GIT%" -C "./%BUILD_VS_LIBDIR%" lfs pull
		if errorlevel 1 (
			set /p LibRetry= "Error during download, retry? y/n"
			if /I "!LibRetry!"=="Y" (
				goto RETRY
			)
                        echo.
		        echo Error: Download of external libraries failed. 
		        echo Until this is resolved you CANNOT make a successful blender build.
		        echo.
		        exit /b 1
		)
	)
) else (
	if NOT EXIST %PYTHON% (
		if not "%SVN%"=="" (
			echo.
			echo Python not found in external libraries, updating to latest version
			echo.
			"%SVN%" update %BUILD_VS_LIBDIR%
		)
	)
)

if NOT EXIST %BUILD_VS_LIBDIR% (
	echo.
	echo Error: Required libraries not found at "%BUILD_VS_LIBDIR%"
	echo This is needed for building, aborting!
	echo.
	if "%SVN%"=="" (
		echo This is most likely caused by svn.exe not being available.
	)
	exit /b 1
)
