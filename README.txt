MODELLING OF A SPACECRAFT HEAT SHIELD TILE

List of files:
RocketTileUI.mLapp
calctemp.m
plottempautoo.m
Analyser.m
tileThickness.m
temp468.jpg
temp480.jpg
temp502.jpg
temp590.jpg
temp597.jpg
temp711.jpg
temp730.jpg
temp850.jpg
Rocket.png

This code is used to model the flow of temperature through a space shuttle tile during the re-entry stage.
The primary way to interact with this software is using the graphical user interface RocketTileUI.mlapp.

Alternatively calctemp.m can be used to obtain data temperature again time and axial thickness,
plottempautoo.m automatically scans temperature data images,
Analyser.m can be used to analyse the accuracy of forward differencing, DuFort-Frankel,
    backward differencing and Crank-Nicolson numerical methods.
tileThickness.m can be used to determine the optimal thickness for a given temperature.
Note: In some cases, temperature data may need to be loaded into the workspace beforehand.