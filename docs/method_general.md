
## Basic Design of the Simulation

The code for the simulation is running in an animation loop. Each iteration of the loop will update and display the results of the simulation on a canvas displayed in a browser. Specifically, the animation loop will clear the canvas, draw the shapefile, draw food patches if added, compute the interindividual distances among all fish in the simulation, update their 3D positions according to the state-dependent forces acting on the fish and finally update the state of all fish. When the app begins the animation loop immediately begins to run using default settings. Essentially each iteration of the loop will check for a large set of parameters to determine the behaviour of the simulation. Users can update the parameters using a graphic user interface menu (built with the [lil-gui library](https://lil-gui.georgealways.com/)). For example, the animation loop keep track of all the fish in the animation, which is determined by the number of fish selected by the user. Once the user updates the number of fish, the animated fish will appear on the screen. Similarly, the parameters that control the state of the simulated fish can be updated in the gui, and will change the behaviour of simulation as it is running.  Below we describe in more detail how the state-switching process is implemented, and how the movement of the fish is determined. 

The positions of each fish are recorded, along with the fish's associated state and force-based parameters. Up to a maximum of 10 000 positions per fish are recorded before they begin to be overwritten. This labelled data can be exported in a csv file, by clicking download tracks in the gui menu. To simulate a more realistic acoustic telemetry sampling process it is possible to subsample the trajectory of the fish (e.g. a transmitter pings every 10 seconds), and to apply a random loss of detections with the "Detection Yield (%)" input. Finally, we can apply a positioning error, to simulate imprecise locations from the fish. In this way the labelled trajectories better matches the processes that generate empirical data and may be more informative for model evaluations and training. The positioning error is applied using a double-normal distribution, where there is a high chance (97%) that the narrow normal distribution ($\sigma = 2$) drawing the error component and leading to accurate positions, and a low chance (3%) of a long-tained distribution ($\sigma = 50$) drawing the error component leading the rare large errors.

The simulation begins with default shapefile, however users can input shapefiles to produce simulations within other waterbodies. Importantly, the simulation requires a .geojson file in coordinates that could be considered in units of metres. We recommend using a UTM coordinate system. This simulation scales the speed of the fish to the dimensions of the waterbody, so it will take a fish longer to swim across a larger area (i.e. it will appear to move slower in a larger shapefile). When a new shapefile is loaded the simulation is reset and will begin from the start. Our code initializes new fish and forces the initial locations of the simulated fish to be within the shape file. We have also implemented a shoreline repulsion force, which is explained in more detail later in the description. 

The code is written in a modular style allowing for future feature integration. The code is available here, and we welcome suggestions or collaboration for improving this simulation tool. 

<p align="center">
  <img src="AnimationLoopDiagram.png" alt="alt text">
  <br>
   <b>Figure 1.</b> A schematic diagram of the animation loop code, indicating the basic order of steps taken in each iteration of the simulation. </centre>
</p>




