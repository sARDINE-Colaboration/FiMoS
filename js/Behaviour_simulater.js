// Include the source file (not sure if the best approach, but it works)
const script_depth = document.createElement('script');
script_depth.src = './js/depthMapFunctions.js';
document.head.appendChild(script_depth);

// Get the canvas element
var canvas = document.getElementById("trailCanvas");
var ctx = canvas.getContext("2d");
ctx.font = "50px Arial";
ctx.fillStyle = "white"; // Set text color
ctx.textAlign = "center"; 

//second canvas for displaying time and scale
var HUDcanvas = document.getElementById("HUDCanvas");
var HUDctx = HUDcanvas.getContext("2d");
HUDctx.font = "50px Arial";
HUDctx.fillStyle = "white"; // Set text color
HUDctx.textAlign = "center"; 

// Set canvas size
canvas.width = window.innerWidth * 0.66;
canvas.height = window.innerHeight * 0.66;
console.log("window size:", window.innerWidth, window.innerHeight);
console.log("canvas size:", canvas.width, canvas.height);
HUDcanvas.width = window.innerWidth * 0.66;
HUDcanvas.height = window.innerHeight * 0.05;

// Array to store fish
var fishes = [];

//array to store food patches
var food_patches = [];
var n_food_patches = 6;
var food_patch_update_time = 1000;

// set simulation parameter
const dt = 0.05;
const dt_output = 1;
const sim_steps = Math.floor(dt_output / dt);

// variables to display simulation time
var time_passed = 0;

// temporary depth data (TODO: load maxDepth from data)
const maxDepth_ori = -8; // if no DepthMap is loaded, use this value
var maxDepth = -8;

//maximum trajectory length
//var max_iter = 10000;
var depthResolution = 4; // used in depthMapFunctions.js

// Variable to store loaded shape data
var geojson; // Declare the variable here
var geojson_xlimits_ori;
var geojson_ylimits_ori;
var geojson_xlimits;
var geojson_ylimits;
var geojson_scalebar = [];
var geojson_scalebar_length = 200;
var bad_geojson = false; //If true tell the user their geojson file's dimensions are too small
var pixel_per_meter; // unit = [pixel] / [meter] : scale to fit the geojson to the canvas AND keep the aspect ratio (important for simulation and rescaling at downloading tracks)
var x_as_max_extent = true; // if true, rescale x coordinates to fit the canvas width --> Y has offset to be in middle
var rescale_offset = 0; // offset to rescale the coordinates to fit the canvas
// would be more efficient to just save the bbox
var depth_map;
var depth_map_colours;
var renderedDepth;

// shared variables among all fish agents
var dist_matrix = [];

// TODO: transform these into a matrix
//       - # of rows: number of states
//parameter means between 0 and 1
// Define the parameter names
const parameterNames = ['beta', 'v0', 'D_phi', 'D_theta', 'D_v', 'patch_strength', 'strength_att', 'strength_align'];
//each vector of eight values corresponds to five parameters in this order: 
// beta, v0, D_phi, D_theta, D_v, patch_strength, strength_att, strength_align
const zero_vector =  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
const resting_vector =  [0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // 0.0, 0.2, 0.2];
const foraging_vector = [0.1, 0.4,  0.9, 0.9, 0.8, 0.0, 0.0, 0.0]; // 0.9, 0.2, 0.2];
const active_vector =   [0.8, 0.8,  0.1, 0.1, 0.2, 0.0, 0.0, 0.0]; // 0.5, 0.2, 0.2];
const resting_sojourn = [20]; // time in minutes (initially 100)
const foraging_sojourn = [20]; // initially 30 
const active_sojourn = [20]; // initially 30

// create a matrix with all available state vectors
let all_states = ['resting', 'foraging', 'active'];
let all_state_means = [resting_vector, foraging_vector, active_vector];
let all_sojourn_times = [resting_sojourn, foraging_sojourn, active_sojourn];

// create state_means_matrix that contains the active vectors resting_vector, ...
var state_means_matrix = [resting_vector];
var states_present = ['resting']; // if more added by gui looks like ['resting', 'active', 'foraging']
var substates = [[2]]; // number of initial substates for resting



// dwell time per state in minutes
var sojourn_times = [resting_sojourn];  // dwell times for resting, foraging, active
var transition_probs = [ //only initial
    [0, 0.5, 0.5],  // probabilities from resting to resting, foraging, active
    [0.5, 0, 0.5],  // probabilities from foraging to resting, foraging, active
    [0.5, 0.5, 0]   // probabilities from active to resting, foraging, active
];
function generateTransitionProbs(dimension) {
    var transition_probs = new Array(dimension).fill(0).map(() => new Array(dimension).fill(0));
    var prob = 1 / (dimension - 1);

    for (var i = 0; i < dimension; i++) {
        for (var j = 0; j < dimension; j++) {
            if (i !== j) {
                transition_probs[i][j] = prob;
            }
        }
    }

    return transition_probs;
}


function calculateTransitionMatrix(transition_probs, sojourn_times){
    var prob_staying = [];
    for (var state in transition_probs){
        prob_staying[state] = 1 - (1 / (sojourn_times[state] * 60 / dt_output)); // sojourn time in dt_output units
    }
    var transition_matrix = [];
    for (var state in transition_probs){
        var probs = transition_probs[state];
        var row = [];
        for (var next_state in probs){
            if (next_state == state){
                row.push(prob_staying[state]);
            }
            else{
                row.push((1-prob_staying[state])*probs[next_state]);
            }
        }
        transition_matrix.push(row);
    }
    return transition_matrix;
}
var transition_matrix = calculateTransitionMatrix(transition_probs, sojourn_times);

//see gui menu at bottom of script for global stdev variable

console.log("variables loaded");

class Vector extends Array {
    add(other) {
      return this.map((e, i) => e + other[i]);
    }
    sub(other) {
        return this.map((e, i) => e - other[i]);
    }
    mul(other) {
      return this.map((e, i) => e * other[i]);
    }
    div(other) {
      return this.map((e, i) => e / other[i]);
    }
    dot(other) {
        return this.reduce((total, current, i) => total + current * other[i], 0);
    }
    add_scalar(scalar) {
      return this.map(e => e + scalar);
    }
    mul_scalar(scalar) {
      return this.map(e => e * scalar);
    }
    div_scalar(scalar) {
      return this.map(e => e / scalar);
    }
    sum() {
        return this.reduce((total, current) => total + current, 0);
    }
    norm() {
        const sumOfSquares = this.reduce((total, current) => total + current ** 2, 0);
        return Math.sqrt(sumOfSquares);
    }
    normXY() {
        return Math.sqrt(this[0] ** 2 + this[1] ** 2);
    }
    phi() {
        return Math.atan2(this[1], this[0]);
    }
    theta() {
        return Math.atan2(this.normXY(), this[2]);
    }
    u_v() {
        return this.div_scalar(this.norm());
    }
    u_phi() {
        const phi = this.phi();
        const sin_theta = Math.sin(this.theta());
        return [-Math.sin(phi) * sin_theta, Math.cos(phi) * sin_theta, 0];
    }
    u_theta() {
        const phi = this.phi();
        const theta = this.theta();
        return [Math.cos(phi) * Math.cos(theta), Math.sin(phi) * Math.cos(theta), - Math.cos(theta)];
    }
  }

// lower and upper limits of state-parameters: beta, v0, D_phi, D_theta, D_v, patch_strength
var lower_limits = new Vector(0.1, 0.01, 0.01, 0.01, 0.01, 0, 0, 0);
var upper_limits = new Vector(3,   3,    1.,   1,    1, 10, 5, 5);


// state = 0; // first ticked state (initially resting)
// beta = 1; // speed relaxation coefficient (the higher, the quicker the preferred speed is reached)
// v0 = 1; // the preferred speed
// D_phi = 0.1; // angular noise intensity (xy-plane)
// D_theta = 0.1; // angular noise intensity (depth-plane)
// D_v = 0.4; // velocity noise intensity
// patch_strength = 0.2; // attraction to food patches
// strength_att = 0.2; // strength of attraction force
// strength_align = 0.2; // strength of alignment force



// Function to switch state
function state_switch(fish){
    if (states_present.length <= 1) {
        // Only one state, no transition needed
        // or no states selected --> no transition
        return;
    }
    var state = fish.state;
    var new_state = state;
    var probs = transition_matrix[state];
    var ran = Math.random();
    var cumulative_prob = 0;
    for (var i = 0; i < probs.length; i++) {
        cumulative_prob += probs[i];
        if (ran < cumulative_prob) {
            new_state = i;
            break;
        }
    }
    if (new_state != fish.state) {
    updateStateParameters(fish, new_state);
    }
}

// function to update the state-dependent parameters of a fish
function updateStateParameters(fish, new_state) {
    if (fish.parameter_array.length > 0) {
        var new_parameters = draw_state_parameter_from_array(fish.parameter_array, new_state);
        fish.state = new_state;
        fish.beta = new_parameters[0];
        fish.v0 = new_parameters[1];
        fish.D_phi = new_parameters[2];
        fish.D_theta = new_parameters[3];
        fish.D_v = new_parameters[4];
        fish.patch_strength = new_parameters[5];
        fish.strength_att = new_parameters[6];
        fish.strength_align = new_parameters[7];
    }
    else {
        fish.v0 = 0;
        fish.D_phi = 0;
        fish.D_theta = 0;
        fish.D_v = 0;
        fish.patch_strength = 0;
        fish.strength_att = 0;
        fish.strength_align = 0;
    }
}

// draw from normal dist using state_means and range(limits)*global_stdev for each fish
function draw_state_parameters(state){
    if (state in state_means_matrix) {
        var parameterVector = new Vector(...state_means_matrix[state].map(mean => {
            var result = randomNormal(mean, param.globalStateStdev);
            return result;
        }));
        parameterVector = parameterVector.map(value => Math.max(0, Math.min(1, value)));      
        parameterVector = parameterVector.mul(upper_limits.sub(lower_limits)).add(lower_limits);
        return parameterVector;
    }
    else {
        console.log("state not in dictionary");
        return false;
    }
}

// draw as many state vectors as there are substates for each (resting, foraging, active)
function draw_state_parameter_array(states_present){
    var state_param_array = [];
    for (var state in states_present){
        state_param_array[state] = [];
        for (let i = 0; i < substates[state]; i++) {
            // Call draw_state_parameters to generate a vector
            let state_param_vector = draw_state_parameters(state);
            // Add the generated vector to state_param_array
            state_param_array[state].push(state_param_vector);
        }
    }
    return state_param_array;
}

// draw one of the state vectors from the state_param_array
function draw_state_parameter_from_array(state_param_array, state){
    // prevent error if state_param_array is empty
    if (state_param_array.length == 0){
        return zero_vector;
    }
    // if not empty, draw a random vector from the array
    var index = Math.floor(Math.random() * substates[state]);
      return state_param_array[state][index];
}

// Fish class
// ATTENTION:
//      in PIXELS: [positions, velocity, position]
//      in METERS: [forces, speed_ms are in METER]
class Fish {
    constructor(id, geojson) {
        this.id = id;
        this.timestamp = [];
        this.timestamp_states = [];
        this.positions = []; // array of positions at each out_step in PIXELS
        this.states = []; //
        this.parameters = [];
        this.position = new Vector(0, 0, 0); // unit = PIXELS
        this.velocity = new Vector(0, 0, 0); // unit = PIXELS / sec
        this.speed_ms = 0.0;  // speed in METER / sec
        this.phi = 0.0;
        this.theta = Math.PI / 2;
        // state-independent parameters (turn-frictions)
        this.turn_phi = 2; // turn friction for the xy-plane angle
        this.turn_theta = 5; // turn friction for the depth angle
        // state-dependent parameters
        this.state = 0;
        // Initialize state-dependent parameters
        this.parameter_array = draw_state_parameter_array(states_present);
        const initial_parameters = draw_state_parameter_from_array(this.parameter_array, this.state);
        // const initial_parameters = draw_state_parameters(this.state);
        this.beta = initial_parameters[0]; // speed relaxation coefficient (the higher, the quicker the preferred speed is reached)
        this.v0 = initial_parameters[1]; // the preferred speed
        this.D_phi = initial_parameters[2]; // angular noise intensity (xy-plane)
        this.D_theta = initial_parameters[3]; // angular noise intensity (depth-plane)
        this.D_v = initial_parameters[4]; // velocity noise intensity
        // social force parameters
        this.strength_att = initial_parameters[6]; // strength of attraction force
        this.strength_align = initial_parameters[7]; // strength of alignment force
        this.social_nn = []; // array of ids of nearest (interacting) neighbours
        // environmental force parameters
        this.patch_strength = initial_parameters[5]; //update based on state
        this.patch_sensing_length = 10; // sensing range for food patches in METER
        this.shore_time = 10;
        this.shore_strength = 5;
        this.depth_time = 2;
        this.depth_stength = 2;
        this.shore_avoidance_sign = 0; // to avoid direction switching during shore-avoidance
        if (geojson) {
            this.createFish(id, geojson);
        }
    }

    createFish(id, geojson) {
        if (geojson) {
            const polygonCoordinates = geojson.features[0].geometry.coordinates;
            let tempx, tempy;
            do {
                tempx = Math.random() * canvas.width;
                tempy = Math.random() * canvas.height;
            } while (!pointInsidePolygon([tempx, tempy], polygonCoordinates));
            this.position[0] = tempx;
            this.position[1] = tempy;
            this.position[2] = 0;
        };
        this.phi = Math.random() * Math.PI * 2;
        this.theta = Math.PI/2 + Math.random() * Math.PI / 8;
        this.speed_ms = Math.max(0, randomNormal(this.v0, this.v0 / 4)); // speed in METER / sec
        this._computeVelocity_inPixelPerSec()
        this.positions.push(this.position);
        this.timestamp.push(0);
        this.timestamp_states.push(0);
        this.states.push(this.state);
        this.parameters.push([this.beta, this.v0, this.D_phi, this.D_theta, this.D_v, this.patch_strength,this.patch_sensing_length, this.strength_att, this.strength_align]);
    }

    updatePosition(polygonCoordinates, dt) {
        // forces
        // ATTENTION: forces are in units of [m/s^2]
        //            --> if your force depends on [distance, speed]
        //                ensure to convert them from pixel to m via pixel_per_meter
        var force = new Vector(0, 0, 0);
        const force_social = socialForce(this, this.strength_att, this.strength_align);
        const force_patch_attraction = patchAttractionForce(this, this.patch_sensing_length, this.patch_strength);
        const force_shore_repulsion = shoreRepulsionForce(this, polygonCoordinates, this.shore_time, this.shore_strength);
        const force_surface_repulsion = surfaceRepulsionForce(this, this.depth_time, this.depth_stength);
        const force_ground_repulsion = groundRepulsionForce(this, this.depth_time, this.depth_stength);
        force = force.add(force_social
            ).add(force_patch_attraction
            ).add(force_shore_repulsion
            ).add(force_surface_repulsion
            ).add(force_ground_repulsion);
        // define unit vectors
        var speed_ms = this.speed_ms, phi = this.phi, theta = this.theta;
        const cos_phi = Math.cos(phi), sin_phi = Math.sin(phi);
        const cos_theta = Math.cos(theta), sin_theta = Math.sin(theta);
        const u_v = new Vector(cos_phi * sin_theta, sin_phi * sin_theta, cos_theta); 
        const u_phi = new Vector(-sin_phi * sin_theta, cos_phi * sin_theta, 0); 
        const u_theta = new Vector(cos_phi * cos_theta, sin_phi * cos_theta, -sin_theta); 
        // modify the speed_ms
        const force_v = force.dot(u_v);
        speed_ms += (this.beta * (this.v0 - speed_ms) + force_v) * dt; //update speed_ms
        speed_ms += Math.sqrt(this.D_v * dt) * randomNormal(0, 1); //add noise to speed_ms
        speed_ms = Math.max(0, Math.min(speed_ms, 20 * this.v0)); //restrict speed_ms between min and max - max is 20*preferred speed_ms
        // modify the xy-plane angle
        const force_phi = force.dot(u_phi);
        const rand_phi = Math.sqrt(this.D_phi * dt) * randomNormal(0, 1)
        phi += (force_phi * dt + rand_phi) / (speed_ms + this.turn_phi);
        // modify the depth angle
        const force_theta = force.dot(u_theta);
        const rand_theta = Math.sqrt(this.D_phi * dt) * randomNormal(0, 1)
        theta += (force_theta * dt + rand_theta) / (speed_ms + this.turn_theta);
        // if the fish makes a looping
        if (theta < 0){
            phi += Math.PI;
            theta *= -1;
        }
        if (theta > Math.PI) {
            phi += Math.PI;
            theta = Math.PI - (theta - Math.PI);
        }
        this.phi = phi, this.theta = theta, this.speed_ms = speed_ms; //apply newly calculated parameters this fish
        // update velocity vector
        this._computeVelocity_inPixelPerSec() // use parameters to update vector
        // Update position based on current velocity
        this.position = this.position.add(this.velocity.mul_scalar(dt));
        this.stayInPolygon(polygonCoordinates);
    }
    
    stayInPolygon(polygonCoordinates) {
        // Check if fish is outside the polygon
        if (!pointInsidePolygon([this.position[0], this.position[1]], polygonCoordinates)) {
            // If outside, reflect back by changing angle
            this.phi += Math.PI; // Reverse direction
            this.velocity = this.velocity.mul_scalar(-1); // Reverse direction for this update
            this.velocity[2] = -this.velocity[2]; // keep depth-direction
            // undo the last position change
            this.position = this.position.add(this.velocity.mul_scalar(3*dt));
            // now the fish is inside the polygon and points in the opposite direction
        }
    }
    _computeVelocity_inPixelPerSec() {
        this.velocity[0] = Math.cos(this.phi) * Math.sin(this.theta);
        this.velocity[1] = Math.sin(this.phi) * Math.sin(this.theta);
        this.velocity[2] = Math.cos(this.theta);
        this.velocity = this.velocity.mul_scalar(this.speed_ms * pixel_per_meter);
    }
    draw_and_save_position(max_iter) {
        const timestamp_now = this.timestamp[this.timestamp.length-1] + 1; 

        // save the position
        this.positions.push(this.position);
        this.timestamp.push(timestamp_now);
        
        if (this.positions.length > max_iter) {
            this.positions.shift();
            this.timestamp.shift();
        }

        // save the state if state-switch happened
        if (this.state != this.states[this.states.length-1]) {
            this.states.push(this.state);
            this.timestamp_states.push(timestamp_now);
            this.parameters.push([this.beta, this.v0, this.D_phi, this.D_theta, this.D_v, this.patch_strength,this.patch_sensing_length, this.strength_att, this.strength_align]);

            if ( (timestamp_now - this.timestamp_states[0]) > max_iter) {
                this.states.shift();
                this.timestamp_states.shift();
                this.parameters.shift();
            }
        }

        // Calculate hue value based on fish state
        // let hue = (this.speed / 10) * 240; // Hue ranges from 0 to 240 (blue to red)
        let hue = this.state / (states_present.length) * 240; // Hue ranges from 0 to 240 (blue to red)  //CTM - removed states_present.length - 1   to avoid the fish changing color when only one state is present
        // Ensure hue stays within the range [0, 360]
        hue = hue % 360;

        // Set color based on fish state 
        ctx.fillStyle = "hsl(" + hue + ", 100%, 50%)";
        // ctx.fillStyle = "white";

        // Calculate the depth ratio
        const depthRatio = Math.abs(this.position[2] / maxDepth); // Adjust maxDepth according to your global variable

        // Adjust the size of the triangle based on the depth
        const triangleSize = 15 - depthRatio * 10; // Decrease triangle size as the depth increases

        // Draw a triangle representing the fish
        ctx.beginPath();
        const x = this.position[0], y = this.position[1];
        ctx.moveTo(x, y);
        const angleA = this.phi + Math.PI * (15 / 18); // 120 degrees
        const angleB = this.phi + Math.PI * (21  / 18); // 240 degrees
        const angleC = this.phi + Math.PI; // 180 degrees
        ctx.lineTo(x + Math.cos(angleA) * triangleSize, y + Math.sin(angleA) * triangleSize);
        ctx.lineTo(x + Math.cos(angleC) * triangleSize*0.5, y + Math.sin(angleC) * triangleSize*0.5);
        ctx.lineTo(x + Math.cos(angleB) * triangleSize, y + Math.sin(angleB) * triangleSize);
        ctx.fill();
        
        // Draw lines to social neighbours
        for (let i = 0; i < this.social_nn.length; i++) {
            const neighbour = fishes[this.social_nn[i]];
            ctx.beginPath();
            ctx.moveTo(this.position[0], this.position[1]);
            ctx.lineTo(neighbour.position[0], neighbour.position[1]);
            ctx.strokeStyle = "rgba(0, 0, 0, 0.5)"; // Set the alpha value to 0.5
            ctx.stroke();
            ctx.fill();
        }
    }
}

// this function computes the distance matrix between all fish in a sparse way
function compute_dist_matrix_sparse(fishes){
    dist_matrix = [];
    var row_first = [];
    for (let i = 0; i < fishes.length; i++) {
        const row = [];
        for (let j = i+1; j < fishes.length; j++) {
            const distance = fishes[i].position.sub(fishes[j].position).norm();
            row.push(distance);
        }
        // fill up the matrix with the already computed distances, and zero for the diagonal
        let row_first = [];
        for (let k = 0; k < i; k++) {
            row_first.push(dist_matrix[k][i]);
        }
        row_first = row_first.concat(0);
        dist_matrix.push(row_first.concat(row));
    };
}

class FoodPatch{
    constructor(id, geojson) {
        this.id = id;
        this.position = new Vector(0, 0, 0);
        if (geojson) {
            this.createPatch(id, geojson);
        }
    }
        
    createPatch(id, geojson) {
        if (geojson) {
            const polygonCoordinates = geojson.features[0].geometry.coordinates;
            let tempx, tempy;
            do {
                tempx = Math.random() * canvas.width;
                tempy = Math.random() * canvas.height;
            } while (!pointInsidePolygon([tempx, tempy], polygonCoordinates));
            this.position[0] = tempx;
            this.position[1] = tempy;
            this.position[2] = 0;
        };
    }

    updatePatch() {
        if (geojson) {
            const polygonCoordinates = geojson.features[0].geometry.coordinates;
    
            // Determine the chance of updating the position
            const updateChance = 1 / food_patch_update_time; // Example: if updateTime = 100, the chance is 1%
    
            // Check if the patch should update based on the update chance
            if (Math.random() < updateChance) {
                let tempx, tempy;
                do {
                    tempx = Math.random() * canvas.width;
                    tempy = Math.random() * canvas.height;
                } while (!pointInsidePolygon([tempx, tempy], polygonCoordinates));
    
                // Update the patch position
                this.position[0] = tempx;
                this.position[1] = tempy;
                this.position[2] = 0;
            }
        }
    }

    draw() {
        
        // Set color based on fish speed
        ctx.fillStyle = "white";

        // Calculate the depth ratio
        //const depthRatio = Math.abs(this.position[2] / maxDepth); // Adjust maxDepth according to your global variable

        // Adjust the size of the triangle based on the depth
        const patchSize = 5 ;
        // Draw a triangle representing the fish
        ctx.beginPath();
        const x = this.position[0], y = this.position[1];
        ctx.moveTo(x, y);
        ctx.arc(x, y, patchSize, 0, Math.PI * 2);
        ctx.fill();
    }
}



// Function to animate the fish
function animateFish() {
    // Clear canvas
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    HUDctx.clearRect(0, 0, canvas.width, canvas.height);

    console.log("animation loop running, sim_steps", sim_steps, ", N_fishes", fishes.length);

    // Draw loaded shape data if available
    if (geojson) { // Check if geoJSONCoordinates is defined
        drawGeoJSON(geojson);
    }

    //draw the food patches
    food_patches.forEach(function (food_patch) {
        food_patch.updatePatch();
        food_patch.draw();
    });

    // Update the fish position in dt-time-steps
    for (let i = 0; i < sim_steps; i++) {
        // compute distance matrix
        compute_dist_matrix_sparse(fishes)
        // update the position of each fish
        fishes.forEach(function (fish) {
            fish.updatePosition(geojson.features[0].geometry.coordinates, dt);
        });
    }

    // console.log("dist_matrix.length", dist_matrix.length);
    // if (FileSystemWritableFileStream.length > 1) {
    //     console.log("dist_matrix[0].length", dist_matrix[1].length);
    // }

    // draw each fish (at output time step dt_output = 1s)
    fishes.forEach(function (fish) {
        fish.draw_and_save_position(param.maxIter); // important 1. record the position and state according to its past movement 2. switch state
        state_switch(fish); // switch state before updating position
    });
    time_passed += dt_output;
    time_passed_stamp();
    drawScaleBar();

     //send error if geojson file is too small in dimensions
     if(bad_geojson){
        //console.log("running code in here");
        ctx.fillStyle = "white";
        ctx.fillText("Geojson dimensions are too small. Please use UTM coordinates, or a coordinate system where the x and y dimensions are in meters", 
                canvas.width/4, canvas.height/8); 
    }

    // Request next frame
    requestAnimationFrame(animateFish);

    
}

function time_passed_stamp(){
    // draw a timer on the canvas
    // Calculate the time in days, hours, minutes, and seconds
    const secs_passed = time_passed * dt_output;
    const days = Math.floor(secs_passed / (24 * 60 * 60));
    const hours = Math.floor((secs_passed % (24 * 60 * 60)) / (60 * 60));
    const minutes = Math.floor((secs_passed % (60 * 60)) / 60);

    // Format the time as "Days:Hours:Minutes:Seconds"
    const formattedHours = hours.toString().padStart(2, '0');
    const formattedMinutes = minutes.toString().padStart(2, '0');
    const timerText = `${days} Days  ${formattedHours}:${formattedMinutes} (h:m)`;

    // Draw the timer on the canvas
    HUDctx.fillStyle = "white";
    HUDctx.font = '20px Arial'; // Set a clear and readable font size
    HUDctx.fillText(timerText, HUDcanvas.width * 0.15 / 5, HUDcanvas.height*0.5);
} 

function rescaleToGeoJSON(x, y, z) {
    // Scale the coordinates
    if (x_as_max_extent) {
        const scaledX = (x)                  / pixel_per_meter + geojson_xlimits_ori[0];
        const scaledY = (y - rescale_offset) / pixel_per_meter + geojson_ylimits_ori[0];
    } else {
        const scaledX = (x - rescale_offset) / pixel_per_meter + geojson_xlimits_ori[0];
        const scaledY = (y)                  / pixel_per_meter + geojson_ylimits_ori[0];
    }
    const scaledZ = z / pixel_per_meter;
    return [scaledX, scaledY, scaledZ];
}

// Function to rescale coordinates
function rescaleCoordinates2D(coordinates) {
    const minX = geojson_xlimits_ori[0];
    const maxX = geojson_xlimits_ori[1];
    const minY = geojson_ylimits_ori[0];
    const maxY = geojson_ylimits_ori[1];
    
    // first adjust the width such that it fits to the canvas
    pixel_per_meter = canvas.width / (maxX - minX);
    x_as_max_extent = true;
    rescale_offset = ( canvas.height - (maxY-minY)*pixel_per_meter) / 2;

    // now check if the rescaled coordinates fit the height
    // if True: keep the scale
    // if False: recompute the scale using height
    if ((maxY - minY) * pixel_per_meter > canvas.height) {
        x_as_max_extent = false; // rescale y coordinates to fit the canvas height
        pixel_per_meter = canvas.height / (maxY - minY);
        rescale_offset = ( canvas.width - (maxX-minX)*pixel_per_meter ) /2;
        
        //shift x positions to centre
        var rescaled_coords = coordinates.map(coord => [
            (coord[0] - minX) * pixel_per_meter + rescale_offset,
            canvas.height - (coord[1] - minY) * pixel_per_meter]);
        
    }else{
        //otherwise shift Y positions to centre
        var rescaled_coords = coordinates.map(coord => [
            (coord[0] - minX) * pixel_per_meter,
            canvas.height - ( (coord[1] - minY) * pixel_per_meter + rescale_offset ) ] );
    }

    geojson_scalebar_length = Math.round(((maxX-minX)*0.15)/10)*10;

    
    console.log("the canvas scale is (pixel_per_meter, xlen, ylen)",
        pixel_per_meter, (maxX-minX), (maxY-minY));
    return rescaled_coords;
}

// Function to rescale x,y coordinates but keep z untouched
function rescaleCoordinates3D(coordinates) {
    const minX = geojson_xlimits_ori[0];
    const maxX = geojson_xlimits_ori[1];
    const minY = geojson_ylimits_ori[0];
    const maxY = geojson_ylimits_ori[1];
    
    // ATTENTION: this function assumes that rescaleCoordinates2D was called before
    // --> pixel_per_meter, x_as_max_extent, rescale_offset are already defined

    if (x_as_max_extent) {
        //shift y positions to centre
        var rescaled_coords = coordinates.map(coord => [
            (coord[0] - minX) * pixel_per_meter,
            canvas.height - ( (coord[1] - minY) * pixel_per_meter + rescale_offset ),
            coord[2] * pixel_per_meter ] );
    } else {
        //shift x positions to centre
        var rescaled_coords = coordinates.map(coord => [
            (coord[0] - minX) * pixel_per_meter + rescale_offset,
            canvas.height - (coord[1] - minY) * pixel_per_meter,
            coord[2] * pixel_per_meter ] );
    }

    return rescaled_coords;
}

function drawGeoJSON(geojson) {
    const features = geojson.features;
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Render depth_map_colours

    ctx.putImageData(renderedDepth, 0, 0);

    ctx.strokeStyle = 'white'; // Set stroke color to white
    features.forEach(feature => {
        const geometry = feature.geometry;
        const type = geometry.type;
        const coordinates = geometry.coordinates;

        if (type === 'LineString') {
            drawLineString(coordinates);
        }
    });

}

// Function to draw a line shape
function drawLineString(coordinates) {
    ctx.beginPath();
    ctx.lineWidth = 2;
    ctx.moveTo(coordinates[0][0], coordinates[0][1]);
    for (let i = 1; i < coordinates.length; i++) {
        ctx.lineTo(coordinates[i][0], coordinates[i][1]);
    }
    ctx.stroke();
}
function depthMapFromData(geojson) {
    // get a grid detailed depth map with depthResolution as grid-distance
    depth_map = sample_depth_map(geojson)
    const flatDepths = depth_map.flat().filter(d => !isNaN(d));
    maxDepth = Math.min(...flatDepths);
    depth_map_colours = getDepthMapColors(depth_map);
    renderedDepth = makeDepthMapImageData(depth_map_colours);
}

function depthMapFromShoreDistance(geojson) {
    maxDepth = 1 * maxDepth_ori; // reset the maxDepth to the original value
    // depth map
    const depthMatrix = calculateDistanceToPolygon(geojson);
    depth_map = mapDepth(depthMatrix);
    depth_map_colours = getDepthMapColors(depth_map);
    renderedDepth = makeDepthMapImageData(depth_map_colours);
}


function resetAndCreateFish(num, geojson) {
    if (!isNaN(num)) {
        // Clear existing fish
        fishes = [];
        time_passed = 0;
        // Create specified number of fish
        for (var i = 0; i < num; i++) {
            fishes.push(new Fish(i, geojson));
        }
    }
}

// Function to handle file selection

function handleGeoJSONFile(event) {
    const file = event.target.files[0];
    const reader = new FileReader();

    reader.onload = function (e) {
        const result = e.target.result;
        try {
            geojson = JSON.parse(result);
            IntegrateGeoJSONShapeAndDepth();

            //reset and create new fish to keep them inside the new map
            var num = parseInt(param.numberOfFish);
            resetAndCreateFish(num, geojson);
            
        } catch (error) {
            console.error("Error parsing GeoJSON:", error);
        }
    };

    reader.readAsText(file);
}


function TestGeoJSONFile() {
    // the former version of this function only checked the limits of the first feature
    bad_geojson = false; // Reset to false
    const minX = geojson_xlimits_ori[0];
    const maxX = geojson_xlimits_ori[1];
    const minY = geojson_ylimits_ori[0];
    const maxY = geojson_ylimits_ori[1];

    console.log(minX);
    console.log("The x dimension is", (maxX - minX));
    if ((maxX - minX) < 5 || (maxY - minY) < 5) {
        bad_geojson = true;
        console.log("Bad GeoJSON");
    }
    return bad_geojson;
}


function IntegrateGeoJSONShapeAndDepth() {
    // Convert geometry to LineString if it's not already
    if (geojson.features.length > 0 && geojson.features[0].geometry.type !== 'LineString') {
        const coordinates = geojson.features[0].geometry.coordinates[0];
        geojson.features[0].geometry = {
            type: 'LineString',
            coordinates: coordinates
        };
    }

    // get limits of the original GeoJSON
    const bbox = calculateBbox(geojson)
    geojson_xlimits_ori = [bbox[0], bbox[2]];
    geojson_ylimits_ori = [bbox[1], bbox[3]];

    // Testing if the GeoJSON is valid
    bad_geojson = TestGeoJSONFile();

    // Rescale coordinates after loading GeoJSON
    canvas.width = window.innerWidth * 0.66; // Reset canvas dimensions
    canvas.height = window.innerHeight * 0.66;

    // Rescale coordinates if geometry is LineString
    const feature0 = geojson.features[0];
    if (feature0.geometry.type === 'LineString') {
        feature0.geometry.coordinates = rescaleCoordinates2D(feature0.geometry.coordinates);
    }

    // Compute limits
    const limits = calculateMinMax_2D_matrix(feature0.geometry.coordinates);
    geojson_xlimits = [limits[0], limits[2]];
    geojson_ylimits = [limits[1], limits[3]];
    console.log(geojson_xlimits);

    //setup scale bar
    calculateScaleBar(geojson_scalebar_length);

    // depth map (if available use it, otherwise compute from shore distance)
    if (geojson.features.length > 1) {
        const feature1 = geojson.features[1];
        if (feature1.geometry.type === 'MultiPoint') {
            // depth map from MultiPoint geometry
            feature1.geometry.coordinates = rescaleCoordinates3D(feature1.geometry.coordinates);
            depthMapFromData(geojson);
            // const depthPoints = feature1.geometry.coordinates;
        } else {
            console.warn("Expected MultiPoint geometry for depth map, but found:", feature1.geometry.type);
        }
    } else {
        depthMapFromShoreDistance(geojson); 
    }

    //make new food patches
    makeFoodPatches();

    // Draw the GeoJSON
    drawGeoJSON(geojson);
}


function InitialGeoJSONFile(filePath) {
    // Fetch the GeoJSON file from the provided filePath
    fetch(filePath)
    .then(response => {
        if (!response.ok) {
            throw new Error('Network response was not ok');
        }
        // Correctly parse the response as JSON
        return response.json(); // Parse as JSON directly
    })
        .then(the_map => {
            geojson = the_map;
            IntegrateGeoJSONShapeAndDepth();
        })
        .catch(error => {
            console.error("Error fetching and parsing GeoJSON:", error);
        });
}

function calculateScaleBar(n_meters){
    geojson_scalebar[0] = canvas.width*0.8 //x1
    geojson_scalebar[1] = HUDcanvas.height*0.5 //y1
    geojson_scalebar[2] = canvas.width*0.8 + pixel_per_meter*n_meters //x2
    geojson_scalebar[3] = HUDcanvas.height*0.5 //y1
}

//function to draw scale bar in "m"
function drawScaleBar(){
    HUDctx.beginPath();                     
    HUDctx.moveTo(geojson_scalebar[0], geojson_scalebar[1]);
    HUDctx.lineTo(geojson_scalebar[2], geojson_scalebar[3]); 
    HUDctx.lineWidth = 2;                   
    HUDctx.strokeStyle = 'white';           
    HUDctx.stroke(); 
    
    const text = `${geojson_scalebar_length} Meters`; 
    const textX = geojson_scalebar[2] + 10; // Position text slightly to the right of the bar (adjust as needed)
    const textY = geojson_scalebar[3]; // Align text vertically with the scale bar

    HUDctx.font = '16px Arial'; // Set font style and size
    HUDctx.fillStyle = 'white';   // Set text color
    HUDctx.textBaseline = 'middle'; // Align text vertically to the middle of the scale bar
    HUDctx.fillText(text, textX, textY); // Draw the text
}



// Function to compute the min and max values of a polygon
function calculateMinMax_2D_matrix(polygon) {
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    polygon.forEach(coord => {
        minX = Math.min(minX, coord[0]);
        minY = Math.min(minY, coord[1]);
        maxX = Math.max(maxX, coord[0]);
        maxY = Math.max(maxY, coord[1]);
    });
    return [minX, minY, maxX, maxY];
}

// function to get a specific column from a 2D array
function getColumn(matrix, col) {
    var coll = matrix.map(d => d[col]);
    return coll;
}


function calculateBbox(thegeojson) {
    let bbox = [Infinity, Infinity, -Infinity, -Infinity];
    if (thegeojson.type === "FeatureCollection") {
        thegeojson.features.forEach(feature => {
            const coords = feature.geometry.coordinates;
            if (feature.geometry.type === "Point") {
                bbox[0] = Math.min(bbox[0], coords[0]);
                bbox[1] = Math.min(bbox[1], coords[1]);
                bbox[2] = Math.max(bbox[2], coords[0]);
                bbox[3] = Math.max(bbox[3], coords[1]);
            } else if (feature.geometry.type === "LineString" ||
                       feature.geometry.type === "Polygon" ||
                       feature.geometry.type === "MultiPoint") {
                console.log("calculating bbox for", feature.geometry.type);
                coords.forEach(coord => {
                    bbox[0] = Math.min(bbox[0], coord[0]);
                    bbox[1] = Math.min(bbox[1], coord[1]);
                    bbox[2] = Math.max(bbox[2], coord[0]);
                    bbox[3] = Math.max(bbox[3], coord[1]);
                });
            }
        });
    }
    return bbox;
}

// function to calculate the intersection of two lines
// base on https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
// section "Given two points on each line segment"
// idea: p1 = position, p2 = position + velocity * time_threshold
//       p3 = polygon[i], p4 = polygon[i+1]
function lineSegmentIntersection(p1, p2, p3, p4) {
    const x1 = p1[0], y1 = p1[1];
    const x2 = p2[0], y2 = p2[1];
    const x3 = p3[0], y3 = p3[1];
    const x4 = p4[0], y4 = p4[1];
    var nominator = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4);
    const denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    const t = nominator / denominator;
    if ((t >= 0) && (t <= 1)) {
        nominator = (x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3);
        const u = - nominator / denominator;
        if ((u >= 0) && (u <= 1)) {
            return [t, x1 + t * (x2 - x1), y1 + t * (y2 - y1)];
        }
    }
    return false
}

function surfaceRepulsionForce(fish, response_time, strength) {
    const p = fish.position, v = fish.velocity, phi = fish.phi, theta = fish.theta;
    const p_after_respTime = p.add(v.mul_scalar(response_time));
    var force = new Vector(0, 0, 0);
    const depth = p_after_respTime[2]; // depth > 0: fish in the air
    if (depth > 0) {
        force[2] = -strength * depth;
    }
    return force;
}

// Dummy function that assumes fixed depth irrespective of the position 
function groundRepulsionForce(fish, response_time, strength) {
    const p = fish.position, v = fish.velocity, phi = fish.phi, theta = fish.theta;
    const p_after_respTime = p.add(v.mul_scalar(response_time));
    var force = new Vector(0, 0, 0);
    const depth = p_after_respTime[2] - maxDepth; // depth < 0: fish in the earth 
    if (depth < 0) {
        force[2] = - strength * depth;
    }
    return force;
}

function shoreRepulsionForce(fish, polygon, response_time, strength) {
    const p = fish.position, v = fish.velocity, phi = fish.phi, theta = fish.theta;
    const p_after_respTime = p.add(v.mul_scalar(response_time));
    var force = new Vector(0, 0, 0);
    var intersection_distance_in_t = 1.1; // if t>1 the lines do not intersect
    var index_line_segment = null;
    // get the closest intersection
    for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
        let intersection = lineSegmentIntersection(p, p_after_respTime, polygon[i], polygon[j]); 
        if (intersection != false) {
            if (intersection[0] < intersection_distance_in_t){
                intersection_distance_in_t = intersection[0];
                index_line_segment = j;
            }
        }
    }
    // if there was an intersection --> compute the vector pointing away from the shore-line-segment
    if (index_line_segment !== null){
        const j = index_line_segment;
        const p1 = polygon[j], p2 = polygon[j+1];
        const vel_seg = new Vector(p1[0] - p2[0], p1[1]- p2[1]);
        var phi_force = vel_seg.phi() + Math.PI/2; // meant to point away from the shore + towards the fish
        // if norm points not towards the fish, flip the direction
        if (Math.cos(phi_force) * Math.cos(phi) + Math.sin(phi_force) * Math.sin(phi) > 0){
            phi_force += Math.PI;
        }
        if (fish.shore_avoidance_sign != 0){
            // point perpendicular to fish-swimming direction + away from shore
            var phi_force2 = fish.phi + fish.shore_avoidance_sign * Math.PI/2; 
        }
        else{
            var phi_force2 = fish.phi + Math.PI/2; 
            fish.shore_avoidance_sign = 1;
            // if phi_force2 points not away from shore (=in same direction as phi_force), flip the direction 
            if (Math.cos(phi_force) * Math.cos(phi_force2) + Math.sin(phi_force) * Math.sin(phi_force2) < 0){
                phi_force2 += Math.PI;
                fish.shore_avoidance_sign = -1;
            }
        }
        // force[0] = Math.cos(phi_force) + 2 * Math.cos(phi_force2); 
        // force[1] = Math.sin(phi_force) + 2 * Math.sin(phi_force2); 
        force[0] = Math.cos(phi_force2); 
        force[1] = Math.sin(phi_force2); 
        force = force.mul_scalar( strength * (1 - intersection_distance_in_t) / force.norm());
    }
    else {
        fish.shore_avoidance_sign = 0;
    }
    return force;
}

function patchAttractionForce(fish, patch_sensing_length, strength) {
    const p = fish.position;
    let force = new Vector(0, 0, 0);
    let dist_threshold = patch_sensing_length * pixel_per_meter; // if t > 1, the lines do not intersect
    let index_patch = null;

    // Get the closest patch
    food_patches.forEach((patch, index) => {
        let patch_dist = p.sub(patch.position).norm(); 
        if (patch_dist < dist_threshold) {
            dist_threshold = patch_dist;
            index_patch = index;
        }
    });

    // If there is a closest patch, compute the vector pointing toward it
    if (index_patch !== null) {
        const closest_patch = food_patches[index_patch];
        var direction_to_patch = closest_patch.position.sub(p);
        const magnitude = direction_to_patch.norm()
        if (magnitude > 0) {
            direction_to_patch = direction_to_patch.div_scalar(magnitude);
        }

        // Scale the force based on strength and proximity
        force = direction_to_patch.mul_scalar(strength); // Adjust the scaling factor as needed
    }
    return force;
}

function socialForce(fish, strength_att, strength_align) {
    // the sensory range depends on the current speed (slow fish have a smaller range)
    const p = fish.position, v_meter = fish.velocity.div_scalar(pixel_per_meter);
    const t_min = 2;
    const t_max = 10;
    const d_min = t_min * fish.speed_ms * pixel_per_meter;
    const d_max = t_max * fish.speed_ms * pixel_per_meter;
    const dist_to_others = dist_matrix[fish.id];
    var force_rep = new Vector(0, 0, 0);
    var force_att = new Vector(0, 0, 0);
    var force_ali = new Vector(0, 0, 0);
    var counter_rep = 0;
    fish.social_nn = [];
    for (let i = 0; i < dist_to_others.length; i++) {
        let d = dist_to_others[i];
        if ( d < d_max && i != fish.id) {
            let dist_weight = 1 - (d - d_min) / (d_max - d_min)
            fish.social_nn.push(i);
            // ATTRACTION or REPULSION FORCE
            let direction_to_other = fishes[i].position.sub(p);
            if (d > 0) {
                direction_to_other = direction_to_other.div_scalar(d);
            }
            direction_to_other = force_rep.add(direction_to_other.mul_scalar(strength_att * dist_weight));
            if (d < d_min) { // REPULSION FORCE
                counter_rep += 1;
                force_rep = force_rep.add(direction_to_other.mul_scalar(- 1));
            }
            else {
                // ATTRACTION FORCE
                force_att = force_att.add(direction_to_other);
                // ALIGNMENT FORCE
                v_meter_other = fishes[i].velocity.div_scalar(pixel_per_meter);
                v_diff = v_meter_other.sub(v_meter);
                force_ali = force_ali.add(v_diff.mul_scalar(strength_align * dist_weight));
            }
        }
    }
    if (counter_rep > 0) {
        return force_rep.div_scalar(counter_rep);
    }
    else {  // ATTRACTION + ALIGNMENT FORCE
        return force_att.add(force_ali);
    }
}



// Function to calculate distance between two points
function calculateDistance(x1, y1, x2, y2) {
    return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
}

function randomNormal(mean, sd) {
    const u = 1 - Math.random(); //Converting [0,1) to (0,1]
    const v = Math.random();
    const num = Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
    return num * sd + mean;
}

function makeFoodPatches(){
         if (geojson) {
             // Clear existing patches 
             food_patches = [];
             // Create specified number of patches 
             for (var i = 0; i < (n_food_patches-1); i++) {
                 food_patches.push(new FoodPatch(i + 1, geojson));
             }
         } else {
             console.log("no geojson file");
         }
}


// Function to check if a point lies inside a polygon using ray casting algorithm
// note: this function could be generalized to return the distance to the shore
//       (in the second step the x-distance is already computed)
function pointInsidePolygon(point, polygon) {
    const x = point[0];
    const y = point[1];
    let inside = false;
    for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
        const xi = polygon[i][0], yi = polygon[i][1];
        const xj = polygon[j][0], yj = polygon[j][1];
        const intersect = ((yi > y) !== (yj > y)) &&
            (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }
    return inside;
}

// Function to subsample fish tracks
function subsampleTracks(tracks, ssTrackpar) {
    // PPK-Note (ToDo): it should be (ssTrackpar <= dt_output), currently it is ok since dt_output = 1....
    if (ssTrackpar <= 1) {
        return tracks;
    } else {
        const indices = Array.from({ length: Math.ceil(tracks.length / ssTrackpar) }, (_, i) => i * ssTrackpar);
        return indices.map(index => tracks[index]);
    }
}

// Function to randomly remove rows based on detection yield
function applyMissedDetections(tracks, detYieldpar) {
    if (detYieldpar >= 100) {
        return tracks;
    } else {
        const indices = tracks.map((_, index) => Math.random() * 100 <= detYieldpar ? index : null).filter(index => index !== null);
        return indices.map(index => tracks[index]);
    }
}

// Function to apply positioning error
// double normal distribution - where there is a high chance that the narrow normal distribution
// is used, but rare events lead to a long tailed distribution. Not perfectly natural
// as the long tail errors are often autocorrelated in time and space, but close for now. 
// error can be toggled off or on with a parameter - toggling occurs in the downloadTracks function
function applyPositioningError(tracks, errorProbHigh, errorSDHigh, errorSDLow) {
    // Implement positioning error logic here
    const withError = tracks.map(track => {
        const trackWithError = track.map(point => {
            if (Math.random() < errorProbHigh) {
                point[0] += randomNormal(0, errorSDHigh);
                point[1] += randomNormal(0, errorSDHigh);
            }else{
                point[0] += randomNormal(0, errorSDLow);
                point[1] += randomNormal(0, errorSDLow);
            }
            return point;
        });
        return trackWithError;
    });
    return withError;
}

function downloadData() {
    downloadTracks();
    downloadStates();
}

function downloadStates() {

    const stateSwitches = fishes.map(fish => {
        return fish.states.map((state, index) => {
            return [state, fish.timestamp_states[index],
                    fish.parameters[index][0], fish.parameters[index][1], fish.parameters[index][2], fish.parameters[index][3], 
                    fish.parameters[index][4], fish.parameters[index][5], fish.parameters[index][6], fish.parameters[index][7],
                    fish.parameters[index][8]];
        });
    });

    // Create CSV content
    let csvContent = "data:text/csv;charset=utf-8,";
    csvContent += "Fish ID,state,Timestamp,beta,v0,D_phi,D_theta,D_v,patch_strength,patch_dist,social_strength,social_align\n";
    stateSwitches.forEach((stateRecord, fishIndex) => {
        stateRecord.forEach(row => {
            csvContent += (fishIndex + 1) + "," + row[0] + "," + row[1] + "," + row[2] + "," + row[3] + "," + row[4] + "," + row[5] + "," + row[6] + "," + row[7] + "," + row[8] + "," + row[9] + "," + row[10] + "\n";
        });
    });

    // Create a link element and trigger download
    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", "fish_states.csv");
    document.body.appendChild(link);
    link.click();

}


//needs to be refactored for fish position update
function downloadTracks() {
    // Rescale fish tracks to original coordinates
    //const ssTrack = parseInt(document.getElementById("ssTrack").value);
    //const detYield = parseInt(document.getElementById("detYield").value);
    //const posErr = document.getElementById("posErr").checked;

    const ssTrack = param.ssTrack;
    const detYield = param.detYield;
    const posErr = param.applyPosErr;

    console.log(posErr);

    const rescaledTracks = fishes.map(fish => {
        return fish.positions.map((position, index) => {
            const originalCoord = rescaleToGeoJSON(position[0], position[1], position[2]);
            return [originalCoord[0], originalCoord[1], originalCoord[2], fish.timestamp[index]];
        });
    });

    // Process tracks: Subsample, apply missed detections, and apply positioning error
    let processedTracks = rescaledTracks.map(track => {
        let subsampledTrack = subsampleTracks(track, ssTrack);
        subsampledTrack = applyMissedDetections(subsampledTrack, detYield);
        if(posErr){
            console.log("checked!");
            subsampledTrack = applyPositioningError(subsampledTrack, 0.97, 2, 50);
        }
        return subsampledTrack;
    });

    // Remove undefined tracks
    processedTracks = processedTracks.filter(track => track && track.length);

    // Create CSV content
    let csvContent = "data:text/csv;charset=utf-8,";
    csvContent += "Fish ID,X,Y,Z,Timestamp\n";
    processedTracks.forEach((track, fishIndex) => {
        track.forEach(row => {
            csvContent += (fishIndex + 1) + "," + row[0] + "," + row[1] + "," + row[2] + "," + row[3] + "\n";
        });
    });

    // Create a link element and trigger download
    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", "fish_tracks.csv");
    document.body.appendChild(link);
    link.click();
}

//initialize Menu
const gui = new GUI();
//gui.add( document, 'title' );
const param = {
  geojsonFile: function() {geoJSONFileInput.click()},
  numberOfFish: 0,
  restingState:true,
  _SubstatesResting: 1, // Use a private variable
  get SubstatesResting() {
    return this._SubstatesResting;
  },
  set SubstatesResting(value) {
    this._SubstatesResting = Math.max(1, value); // Ensure value cannot go below 1
  },
  activeState: false,
  _SubstatesActive: 1, // Use a private variable
  get SubstatesActive() {
    return this._SubstatesActive;
  },
  set SubstatesActive(value) {
    this._SubstatesActive = Math.max(1, value); // Ensure value cannot go below 1
  },
  foragingState: false,
  _SubstatesForaging: 1, // Use a private variable
  get SubstatesForaging() {
    return this._SubstatesForaging;
  },
  set SubstatesForaging(value) {
    this._SubstatesForaging = Math.max(1, value); // Ensure value cannot go below 1
  },
  addCustomState: false,
  globalStateStdev: 0.05,
  ssTrack: 1,
  detYield: 100,
  maxIter: 43200,  // = 12 hours * 60 minutes * 60 seconds
  applyPosErr: false,
  download: function() {downloadData()},
};


//set up gui menu
const geoJSONFileInput = document.createElement('input');
geoJSONFileInput.setAttribute('type', 'file');
geoJSONFileInput.setAttribute('accept', '.geojson');
geoJSONFileInput.style.display = 'none';
document.body.appendChild(geoJSONFileInput);

//allow buttons to flash
function flashButton() {
    this.flash();
}

gui.add(param, 'geojsonFile').name('Select GeoJSON file');
gui.add(param, 'numberOfFish').name('Number of fish').onChange(value => {
    // Handle change
    var num = parseInt(value);
    if (geojson) {
        resetAndCreateFish(num, geojson);
    } else {
        console.log("no geojson file");
    }
});

// Function to add or remove a state
function manageStateChecking(stateName, vector, sojourn, substate, include) {
    
    const N_states_present_before_switch = states_present.length
    if (include) {
        state_means_matrix.push(vector);
        states_present.push(stateName);
        sojourn_times.push(sojourn);
        substates.push(substate);
    } else {
        state_means_matrix = state_means_matrix.filter((_, index) => states_present[index] !== stateName);
        sojourn_times = sojourn_times.filter((_, index) => states_present[index] !== stateName);
        substates = substates.filter((_, index) => states_present[index] !== stateName);
        states_present = states_present.filter(state => state !== stateName);
    }

    // Update transition_probs and transition_matrix
    var dimension = states_present.length;
    transition_probs = generateTransitionProbs(dimension);
    transition_matrix = calculateTransitionMatrix(transition_probs, sojourn_times);

    // Update fish.parameter_array with the new states_present
        fishes.forEach(function (fish) {
            fish.parameter_array = draw_state_parameter_array(states_present);
        }); 

    if ( (!include) || (N_states_present_before_switch == 0)) {
      // update fish to be in state 0 if a state is removed
      // or if no state was present before
        fishes.forEach(function (fish) {
           updateStateParameters(fish, 0);
           });
         }; 
        
}

function updateSubstateEntry(stateName, newNumberSubstates) {
        const state = states_present.indexOf(stateName);
        if (state !== -1) { // Check if the 'resting' exists in states_present
            substates[state] = newNumberSubstates; // Update the corresponding entry in substates

        // Update fish.parameter_array with the new states_present
        fishes.forEach(function (fish) {
            fish.parameter_array[state] = [];
            for (let i = 0; i < substates[state]; i++) {
                // Call draw_state_parameters to generate a vector
                let state_param_vector = draw_state_parameters(state);
                // Add the generated vector to parameter_array
                fish.parameter_array[state].push(state_param_vector);
            }
        }); 
        } else {
            console.log(stateName, "not found in states_present");
        }
}

// Function to add input fields for each entry in the state vector with the specified names
// and to update parameter array for each fish on change and let it restart in state 0
function adjustStateParameters(folder, state) {
  let adjustedState = all_state_means[state];
    adjustedState.forEach((value, index) => {
        folder.add({ [parameterNames[index]]: value }, parameterNames[index])
            .name(parameterNames[index])
            .onChange((newValue) => {
                const parsedValue = parseFloat(newValue);
                if (!isNaN(parsedValue)) {
                    all_state_means[state][index] = parsedValue;
                    adjustedState[index] = parsedValue;
                }
                fishes.forEach(function (fish) {
                    fish.parameter_array = draw_state_parameter_array(states_present);
                    updateStateParameters(fish, 0);
                });
            });
    });


    let adjustedSojourn = all_sojourn_times[state][0];
    console.log("Adding sojournTime to folder:", adjustedSojourn);
    folder.add({ sojournTime: adjustedSojourn }, 'sojournTime')
        .name('Sojourn Time (Min)')
        .onChange((newValue) => {
            const parsedValue = parseFloat(newValue);
            if (!isNaN(parsedValue)) {
                all_sojourn_times[state] = parsedValue;
                adjustedSojourn = parsedValue;
                // replace the updated sojourn time in sojourn_times
                sojourn_times[states_present.indexOf(all_states[state])] = parsedValue;
            }
            // Update transition_probs and transition_matrix
            var dimension = states_present.length;
            transition_probs = generateTransitionProbs(dimension);
            transition_matrix = calculateTransitionMatrix(transition_probs, sojourn_times);
        });
}

gui.add(param, 'restingState').name('Include resting state (y/n)').onChange(value => {
    manageStateChecking(all_states[0], all_state_means[0], all_sojourn_times[0], param.SubstatesResting, value);
});
gui.add(param, 'SubstatesResting').name('Select number of resting states').onChange(value => {
    updateSubstateEntry('resting', value)
});
const restingfolder = gui.addFolder('Adjust resting state');
restingfolder.close(); // Close the folder by default
adjustStateParameters(restingfolder, 0)

gui.add(param, 'activeState').name('Include active state (y/n)').onChange(value => {
    manageStateChecking(all_states[2], all_state_means[2], all_sojourn_times[2], param.SubstatesActive, value);
});
gui.add(param, 'SubstatesActive').name('Select number of active states').onChange(value => {
    updateSubstateEntry('active', value)
});
const activefolder = gui.addFolder('Adjust active state');
activefolder.close(); // Close the folder by default
adjustStateParameters(activefolder, 2)

gui.add(param, 'foragingState').name('Include foraging state (y/n)').onChange(value => {
    manageStateChecking(all_states[1], all_state_means[1], all_sojourn_times[1], param.SubstatesForaging, value);
});
gui.add(param, 'SubstatesForaging').name('Select number of foraging states').onChange(value => {
    updateSubstateEntry('foraging', value)
});
const foragingfolder = gui.addFolder('Adjust foraging state');
foragingfolder.close(); // Close the folder by default
adjustStateParameters(foragingfolder, 1)

// Function to handle adding a custom state
function addCustomState() {
    const newStateIndex = all_states.length;
    // Prompt the user for a state name
    let newStateName = prompt("Please enter a name for the new state:", `state_${newStateIndex}`);
   
    // Add the new state to the states array of all available state
    all_states.push(newStateName);
    all_state_means.push(all_state_means[0]);
    all_sojourn_times.push(all_sojourn_times[0]);

    param[newStateName] = false;
    let newStateSubstatesName = `Substates${newStateName}`;
    const newStateSubstates = 1;
    param[newStateSubstatesName] = newStateSubstates;


    manageStateChecking(all_states[newStateIndex], 
        all_state_means[newStateIndex], 
        all_sojourn_times[newStateIndex], param[newStateSubstatesName], false);
    // then decide whether to include the state or not
    gui.add(param, newStateName).name(`Include ${newStateName} (y/n)`).onChange(value => {
        manageStateChecking(all_states[newStateIndex], 
            all_state_means[newStateIndex], 
            all_sojourn_times[newStateIndex], param[newStateSubstatesName], value);
    });
    gui.add(param, newStateSubstatesName).name(`Select number of ${newStateName} substates`).onChange(value => {
        param[newStateSubstatesName] = Math.max(1, value); // Ensure value cannot go below 1
        updateSubstateEntry(newStateName, param[newStateSubstatesName])
    });
    
    stateFolders[newStateName] = gui.addFolder(`Adjust ${newStateName}`);
    stateFolders[newStateName].close(); // Close the folder by default
    adjustStateParameters(stateFolders[newStateName], newStateIndex)

    // Logic to add a custom state
    console.log("Custom state added");
}

// Add custom states
// Button to add a new state
const stateFolders = {};
gui.add({ addCustomState: addCustomState }, 'addCustomState').name('Add custom state');

// gui.add(param, 'addForagingState').name('add additional foraging state');
gui.add(param, 'globalStateStdev').name('Behavioural variation');
gui.add(param, 'ssTrack').name('Subsample track (s)');
gui.add(param, 'detYield').name('Detection yield (%)');
gui.add(param, 'maxIter').name(`Max. track length (${dt_output} s)`);
gui.add(param, 'applyPosErr').name('Apply position error (y/n)');
gui.add(param, "download").name('Download tracks');
geoJSONFileInput.addEventListener('change', handleGeoJSONFile);

/* console.log("Resting state parameters:", draw_state_parameters(0));
// Test state_switch function
const testFish = new Fish(1, geojson);
console.log("Initial state:", testFish.state);
console.log("Initial parameters:", testFish.beta, testFish.v0, testFish.D_phi, testFish.D_theta, testFish.D_v);

state_switch(testFish);
console.log("After state switch:");
console.log("New state:", testFish.state);
console.log("New parameters:", testFish.beta, testFish.v0, testFish.D_phi, testFish.D_theta, testFish.D_v);


var arraytest = draw_state_parameter_array(states_present);
// Call the draw_state_parameter_array function with the states_present array
console.log("Result of draw_state_parameter_array:", arraytest);
// call result of function draw_state_parameter_from_array
console.log("Result of draw_state_parameter_from_array:", draw_state_parameter_from_array(arraytest, 0));
 */
// Start animation loop for fish
InitialGeoJSONFile("../FiMoS/data/Most_shoreline_polygon_UTM33.geojson");
animateFish();