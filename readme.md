# Domain Transform Filter

Implementaiton of domain transform for edge-aware image and video 
processing [Gastal & Oliveira, SIGGRAPH 2011] in Nuke (13.1v5 et 15.1v5)

Only implement recursive filter (RF) at the moment.

# How to build

## Windows

Depending on Nuke version you will need different compiler:
 - Nuke 13 require [https://aka.ms/vs/15/release/vs_buildtools.exe](Visual Studio 15 2017).
 - Nuke 15 require [https://aka.ms/vs/17/release/vs_buildtools.exe](Visual Studio 17 2022).

Note: For Nuke 13.1v5 you may need to create a `cmake` directory in `Nuke13.1v5` installation root and move `NukeConfig.cmake` inside.
This issue may be related to our installation only, but I thought it was worth to note.


You can now build this plugin with: 

### Nuke 13.1v5
```cmd
cmake -DNuke_DIR="C:\PROGRA~1\Nuke13.1v5\cmake" -G "Visual Studio 15 2017" -A x64 -B </path/to/build/dir> -S <path/to/repository/root>
cmake --build </path/to/build/dir> --config Release
```

### Nuke 15.1v5
```cmd
cmake -DNuke_DIR="C:\PROGRA~1\Nuke15.1v5\cmake" -G "Visual Studio 17 2022" -A x64 -B </path/to/build/dir> -S <path/to/repository/root>
cmake --build </path/to/build/dir> --config Release
```

Once your `.dll` are built, you can just add to your plugin path the way you want.

# How to use it

This plugin implement a single node name `DTFilter`. There is a few parameters you can adjust:
 - Sigma Spatial: Affect filter based on pixel location on image.
 - Sigma Range: Affect filter based on pixel color.
 - Iteration: Number of iterations, 3 is a good starting point and higher value of sigma and higher resolution will require more iterations.

## Unix

WIP