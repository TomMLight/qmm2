// Code throughout adapts, translates, and expands upon original Quantum Marble Maze by Dr. Crispin Cooper
// Taken from "fiftysevendegreesofrad" on GitHub, v1.0 formalized 17-05-2024
// Accessed 25-01-2025
// https://github.com/fiftysevendegreesofrad/quantum

"use strict"; // Directive to indicate code should be executed in strict mode.

//FUNCTION DECLARATIONS --------------------------------------------------------------------------------------------------------------------------

// Unpacks the RGBA data spat out by the hacked GameRender.update() and AmpColourMap.process() methods and repacks it as an ImageData object.
function setRGB(startX, startY, width, height, rgbArray, offset, scansize) {
  const dataArray = new Uint8ClampedArray(width * height * 4); // Create an array of unsigned 8 bit values, with 4 values (R, G, B, and A) for each pixel.
  for(let i = 0; i < width * height; i++) { //Iterate over each value in the incoming array (1 for each pixel):
    //For each pixel, unpack each of the 4 values from the 2nd dimension of input array and put in one large array.
    dataArray[(i * 4) + 0] = rgbArray[i][0];
    dataArray[(i * 4) + 1] = rgbArray[i][1];
    dataArray[(i * 4) + 2] = rgbArray[i][2];
    dataArray[(i * 4) + 3] = rgbArray[i][3];
  }
  return new ImageData(dataArray, width); // Create a new ImageData object using our unpacked array, and return it.
}

// Emulates the System.nanoTime() method from Java.
function nanoTime() {
  return performance.now() * 1_000_000 // Converting milliseconds to nanoseconds.
}

// Partial port of UpdateTask.run() from GameWindow.java in original, hacked to get requestAnimationFrame() to work. Async so sleep functionality works.
async function graphicsLoop() {
  // In the original, there is a main loop that breaks when the goal is sounded, and calls a repaint when a certain amount of time has passed.
  // I have struggled greatly to directly port this to JavaScript due to requestAnimationFrame, which as far as I can tell is the only way to achieve runtime animation.
  // If this is not the case, I will fix this asap.
  // Instead, this loop just executes until it is time to draw a new gfx frame, when it exits and recursively calls requestAnimationFrame again. 
  while(true) {
    quantumframes_this_frame++;
    if(quantumframes_this_frame < quantum_frames_per_gfx_frame) {
      manager.getQD().step();
    } else {
      const timesincelastframe = nanoTime() - lastframetime;
      const sleeptime = gfxframetime - timesincelastframe;
      if(sleeptime > 0) {
        await new Promise(r => setTimeout(r, sleeptime/1000000)); // Sleeps for a number of milliseconds equal to argument provided. This code was taken from Stack Overflow post made by user "Dan Dascalescu" on 18-07-2018, edited by "blackgreen♦" on 03-04-2024. Accessed 27-02-2025. https://stackoverflow.com/questions/951021/what-is-the-javascript-version-of-sleep
      }
    }
    const currenttime = nanoTime();
    if(currenttime - lastframetime > gfxframetime) { //"30fps"
      quantumframes_this_frame = 0;
      lastframetime = currenttime;
      manager.updateGraphics();
      break; // New graphics frame needed, exit the loop and request new frame.
    }
  }
  requestAnimationFrame(graphicsLoop); // Recursively call requestAnimationFrame to queue next frame.
}
//------------------------------------------------------------------------------------------------------------------------------------------------


//CLASS DECLARATIONS -----------------------------------------------------------------------------------------------------------------------------

// # prefix is used to denote private properties/methods in JavaScript
// this. prefix required for all property references in JavaScript, even within same class.
// Can't really figure out how to do final without const, which doesn't work for class properties? If there is some easy way to emulate this, let me know - but I don't think it is actually necessary or even desirable.

// Partial port of class from QuantumData.java in original
class Complex {
  #real; #imag; // final float
  constructor(r = 0, i = 0) {this.#real = r; this.#imag = i;} // Preserving default values of 0 from original.
  mod2() {return this.#real*this.#real+this.#imag*this.#imag} //NOTE TO SELF: MOD2 STANDS FOR MODULUS SQUARED!!!
}

// Partial port of class from QuantumData.java in original
class FloatArray2d {
  #data; // float[]
  #width; // int
  constructor(width, height) {
    this.#data = new Float32Array(width * height).fill(0); //Float32Array closest to float[] in Java. In Java arrays are initialized with 0/False/Null, whereas they are left undefined in JavaScript. fill() fixes this.
    this.#width = width;
  }
  setEqualTo(other) {
    this.#width = other.getWidth();  // getWidth() necessary as in JavaScript private properties cannot be read even by other objects of same class, unlike Java.
    this.#data = [...other.getData()]; // Spread operator used in this way equivalent to .clone() in Java. getData() necessary as in JavaScript private properties cannot be read even by other objects of same class, unlike Java.
  }
  get(x, y) {return this.#data[x+this.#width*y];}
  set(x, y, v) {this.#data[x+this.#width*y]=v;}
  getWidth() {return this.#width;} // Necessary as in JavaScript private properties cannot be read even by other objects of same class, unlike Java.
  getData() {return this.#data;} // Necessary as in JavaScript private properties cannot be read even by other objects of same class, unlike Java.
  setData(v) {this.#data = v;} // Used to copy results of kernel into array. Maybe temp if I use textures?
}

// Partial port of class from QuantumData.java in original
class BoolArray2d {
  #data; // final boolean[]
  #width; // int
  constructor(width, height) {
    this.#data = new Array(width * height).fill(false); //In Java arrays are initialized with 0/False/Null, whereas they are left undefined in JavaScript. fill() fixes this.
    this.#width = width;
  }
  width() {return this.#width;}
  get(x,y) {return this.#data[x+this.#width*y];}
  set(x,y,v) {this.#data[x+this.#width*y]=v;}
  getData() {return this.#data;} // Need to grab entirety of data to load into GPU.
}

// Partial port of class from QuantumData.java in original
class QuantumData {
  #levelDesignPotential; #pot_cache; // final FloatArray2d
  sink_mult; // final FloatArray2d
  #real; #imag; #init_real; #init_imag; // final FloatArray2d
  #walls; #sink; // BoolArray2d
  #delta_t; #maxtilt; // float
  running = false; // boolean

  width; height; // final int
  #controlstate; // ControlState

  #gpu; // To save calling the constructor every time QuantumData.step() is called, a new gpu property is added to the class for future reference.
  #updateCompontent; // Kernel used to update real and imaginary components of simulation in step().
  #addGaussComp; // Kernel used to add gaussian to real and imaginary components of simulation in addGaussian().
  #addWallsComp; // Kernel used to set real and imaginary components of simulation to 0 when a wall is present.

  constructor(width, height, ks) {
    this.#controlstate = ks;
    this.width = width;
    this.height = height;
    this.#real = new FloatArray2d(width, height);
    this.#imag = new FloatArray2d(width, height);
    this.#init_real = new FloatArray2d(width, height);
    this.#init_imag = new FloatArray2d(width, height);
    this.#walls = new BoolArray2d(width, height);
    this.#sink = new BoolArray2d(width, height);
    this.sink_mult = new FloatArray2d(width, height);
    this.#levelDesignPotential = new FloatArray2d(width, height);
    this.#pot_cache = new FloatArray2d(width, height);
    this.#setupKernels(); // Kernels are several lines long, so in my opinion having them stored seperately makes code more readable.
  }
  getReal() {
    return this.#real.getData();
  }
  getImag() {
    return this.#imag.getData();
  }
  #setupKernels() {
    this.#gpu = new GPU.GPU(); // Create GPU object.

    // Next, define kernels:
    
    // updateComponent - Handles updating real and imaginary (two seperate calls) components of simulation in step().
    // Takes in various properties of QuantumData, with FloatArray2d objects being passed as their raw 1D "data" property.
    // Returns a 1D array that can then be copied into the "data" property of that component's FloatArray2d object.
    // "update" and "ref" are the component being updated and merely referenced respectively.
    // "sign" allows a toggle between addition and subtraction so kernel can be reused for both real and imaginary component updates.
    // Simply flip which component is "update" vs "ref" and "sign" from 1 to -1.
    this.#updateCompontent = this.#gpu.createKernel(function(w, h, delta_t, update, ref, walls, sink_mult, pot_cache, sign) {
        // Calculate x and y co-ordinates from index in 1D array.
        const x = this.thread.x % w;
        const y = Math.floor(this.thread.x / w);
    
        // If the co-ordinates are on the border or inside a wall they should not be simulated - return 0.
        if(x < 1 || x >= w - 1 || y < 1 || y >= h - 1 || walls[x+w*y] == 1) {
            return 0;
        }
    
        // Else, return updated value of component.
        return sink_mult[x+w*y] * (update[x+w*y] + sign * delta_t * (-0.5 * (ref[x+w*(y-1)]+ref[x+w*(y+1)]+ref[(x-1)+w*y]+ref[(x+1)+w*y]-4*ref[x+w*y]) + pot_cache[x+w*y]*ref[x+w*y]));
    }).setOutput([this.width * this.height]); // Set output size to be equal to "data" property of FloatArray2d objects.
    
    // addGaussComp - Handles adding real and imaginary components (two seperate calls) of a gaussian to the simulation data.
    // Takes in width and height of simulation, various pre-computed statistics of gaussian, a phase shift, and the component to be updated.
    // Phase shift allows elegant toggle between using sin and cos - allowing the same kernel to be reused for both components.
    // comp = real and phase = pi/2, or comp = imag and phase = 0.
    this.#addGaussComp = this.#gpu.createKernel(function(width, height, xc, yc, a, d, omegax, omegay, phase, comp) {
      // Calculate x and y co-ordinates from index in 1D array.
      const x = this.thread.x % width;
      const y = Math.floor(this.thread.x / width);

      // If the co-ordinates are on the border they should not be simulated - return 0.
      if(x < 1 || x >= width - 1 || y < 1 || y >= height - 1) {
          return 0;
      }

      // Else, update component.
      const r2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);
      const v = a * Math.exp(-r2/d)
          * Math.sin((omegax*x/width) + phase)
          * Math.sin((omegay*y/height) + phase);
      return comp[this.thread.x] + v;

    }).setOutput([this.width * this.height]); // Set output size to be equal to "data" property of FloatArray2d objects.

    // addWallsComp - Handles setting real and imaginary components (two seperate calls) of simulation to 0 when a wall is present in that cell.
    // Takes in width and height of simulation, an array of walls, and the component to be updated.
    // To update both components, simply call the function twice, supplying comp as real component and then imaginary component.
    this.#addWallsComp = this.#gpu.createKernel(function(width, height, walls, comp) {
      // Calculate x and y co-ordinates from index in 1D array.
      const x = this.thread.x % width;
      const y = Math.floor(this.thread.x / width);

      // MAY BE UNECESSARY?
      // If the co-ordinates are on the border they should not be simulated - return 0.
      if(x < 1 || x >= width - 1 || y < 1 || y >= height - 1) {
          return 0;
      }

      // Else check for presence of walls:
      if(walls[this.thread.x] == 1) {
        return 0
      }
      return comp[this.thread.x]; // No walls, return existing value.

    }).setOutput([this.width * this.height]); // Set output size to be equal to "data" property of FloatArray2d objects.
  }
  #saveInitialState() {
    this.#init_real.setEqualTo(this.#real);
    this.#init_imag.setEqualTo(this.#imag);
  }
  resetInitialState() {
    this.#real.setEqualTo(this.#init_real);
    this.#imag.setEqualTo(this.#init_imag);
  }
  setDeltaT(dt) {this.#delta_t = dt;}
  setMaxTilt(mt) {this.#maxtilt = mt;}
  addGaussian(xc, yc, sigma, fx, fy, ascale) {
    const a = ascale * Math.pow(2*Math.PI*sigma*sigma,-0.25);
    const d = 4*sigma*sigma;
    const omegax = 2*Math.PI*fx; // "fixme this seems wrong" (Crispin has confirmed since this is not the case)
    const omegay = 2*Math.PI*fy; // "fixme this seems wrong" (Crispin has confirmed since this is not the case)

    // Run kernels to add Gaussian to real and imaginary components, copying the results into the "data" property of that component's FloatArray2d object.
    this.#real.setData(this.#addGaussComp(this.width, this.height, xc, yc, a, d, omegax, omegay, Math.PI/2, this.#real.getData()));
    this.#imag.setData(this.#addGaussComp(this.width, this.height, xc, yc, a, d, omegax, omegay, 0, this.#imag.getData()));
  }
  #reset_potential_cache() {
    // "potentials >0 are problematic"
    // "pixel wide band with potential +1 above background - tunnelling"
    // "potential of -5 over width of universe - good for steering"

    // "if tilting 2 directions at once reduce tilt to compensate"
    const totalslope = Math.abs(this.#controlstate.getXSlope())+Math.abs(this.#controlstate.getYSlope());
    let tilt = 0;
    if(totalslope <= 1) {
      tilt = this.#maxtilt;
    } else {
      tilt = this.#maxtilt/totalslope;
    }

    // "compute desired relative potentials of corners"
    const biggerdim = Math.max(this.width, this.height);
    const right_change = -this.#controlstate.getXSlope()*tilt*this.width/biggerdim;
    const bottom_change = -this.#controlstate.getYSlope()*tilt*this.height/biggerdim;
    const topleft = -right_change - bottom_change;
    const topright = right_change - bottom_change;
    const bottomleft = -right_change + bottom_change;
    const bottomright = right_change + bottom_change;

    // "adjust all potentials to be <0"
    const max = Math.max(Math.max(Math.max(topleft,topright), bottomleft), bottomright); //In JavaScript Math.max() can actually take an arbitrary number of parameters instead of just 2 in Java, but I have nested several as in the Java original to avoid confusion for now, will change later.
    const newtopleft = topleft - max;

    // "compute per-simulation-element steps in potential to efficiently compute"
    const x_pot_step = right_change/this.width;
    const y_pot_step = bottom_change/this.width;
    let left_edge_pot = newtopleft;
    for(let y=1;y<this.height-1;y++) {
      left_edge_pot += y_pot_step;
      let current_pot = left_edge_pot;
      for(let x=1;x<this.width-1;x++) {
        current_pot += x_pot_step;
        this.#pot_cache.set(x,y,current_pot + this.#levelDesignPotential.get(x,y));
      }
    }
  }
  #ensure_no_positive_potential() {
    let maxpot = Number.NEGATIVE_INFINITY;
		for(let x=0;x<this.width;x++) {
      for(let y=0;y<this.height;y++) {
        const pot = this.#levelDesignPotential.get(x,y)
        if(pot>maxpot) {maxpot = pot};
      }
    }
		for(let x=0;x<this.width;x++) {
      for(let y=0;y<this.height;y++) {
        this.#levelDesignPotential.set(x, y,
            this.#levelDesignPotential.get(x,y)-maxpot);
      }
    }
  }
  #add_walls() {
    // Run kernels to add set real and imaginary components to 0 if a wall is present.
    this.#real.setData(this.#addWallsComp(this.width, this.height, this.#walls.getData(), this.#real.getData()));
    this.#imag.setData(this.#addWallsComp(this.width, this.height, this.#walls.getData(), this.#imag.getData()));
  }
  step() {
    if(!this.running) {
      this.running = true;
      this.#setupSinkMult();
      this.#add_walls(); // "must be done before saveInitialState"
      this.#saveInitialState();
      this.#ensure_no_positive_potential();
    }
    this.#controlstate.step();
    this.#reset_potential_cache();

    // Run kernels to update real and imaginary components, copying the results into the "data" property of that component's FloatArray2d object.
    this.#real.setData(this.#updateCompontent(this.width, this.height, this.#delta_t, this.#real.getData(), this.#imag.getData(), this.#walls.getData(), this.sink_mult.getData(), this.#pot_cache.getData(), 1));
    this.#imag.setData(this.#updateCompontent(this.width, this.height, this.#delta_t, this.#imag.getData(), this.#real.getData(), this.#walls.getData(), this.sink_mult.getData(), this.#pot_cache.getData(), -1)); 

  }
  #setupSinkMult() {
		// "flood fill sink_mult with 0 where not a sink; otherwise distance in pixels from non-sink"
		// "...basically a mini-Dijkstra"
    class Pixel {
      x;y;d;
      constructor(xx, yy, dd){this.x=xx;this.y=yy;this.d=dd;}
      // Instead of Pixel inheriting from comparable and defining a compareTo function as in the original, an equivalent comparator is supplied as a parameter when creating the queue instead.
    }
    const queue = new PriorityQueue((a, b) => a.d < b.d); // Comparator equivalent to Pixel.compareTo in original. I have confirmed with tests that items are pulled in the same order from this and original queue implementation.
    for(let y=0;y<this.height;y++) {
      for(let x=0;x<this.width;x++) {
        this.sink_mult.set(x,y,Number.POSITIVE_INFINITY);
        if(!this.#sink.get(x,y) && !this.#walls.get(x,y)) {
          queue.push(new Pixel(x, y, 0))
        }
      }      
    }
    while(!queue.isEmpty()) {
      const p = queue.pop(); // pop() equivalent to Java poll()
      if(this.sink_mult.get(p.x,p.y) > p.d) {
        this.sink_mult.set(p.x,p.y,p.d);
        for(let dx=-1;dx<=1;dx+=2) {
          for(let dy=-1;dy<=1;dy+=2) {
            const q = new Pixel(p.x+dx,p.y+dy,p.d+1);
            if(q.x>=0 && q.x<this.width && q.y>=0 && q.y<this.height) {
              queue.push(q); // push() equivalent to Java add()
            }
          }          
        }
      }
    }
    // "now convert these to actual sink_mults"
    const suddenness = 0.005;
    for(let y=0;y<this.height;y++) {
      for(let x=0;x<this.width;x++) {
        const dist = this.sink_mult.get(x,y);
        this.sink_mult.set(x,y,Math.exp(-Math.pow(dist/2,2)*suddenness));
      }      
    }
  }
  // Slightly hacky work around to get control state from main function in order to bind event listener.
  getCS() {
    return this.#controlstate;
  }
}

// Partial port of class from QuantumData.java in original
// Apparently interfaces aren't really used in JavaScript and are sort of against the spirit of the language?
// Only using AmpColourMap right now so may add rainbow later, but this will suffice for now.
class AmpColourMap {
  #gamma = 0.7; // float
  #maxindex = 255; // int
  #lookup; // int[]
  #max = 0; // float
  #gain = 0; // float

  #gpu; // To save calling the constructor every time process() is called, a new gpu property is added to the class for future reference.
  // Similarly, kernels are defined once and saved for future use.
  #mod2;
  #renderFrame;
  constructor(size) {
    this.#lookup = new Int32Array(this.#maxindex + 1); //Int32Array closest to int[] in Java.
    for (let i=0;i<this.#maxindex+1;i++) {
      this.#lookup[i] = Math.trunc(255*Math.pow(i/this.#maxindex,this.#gamma)) //Math.trunc() closest to how casting to int works in Java.
    }
    this.#setupGPU(size);
  }
  #setupGPU(size) {
    this.#gpu = new GPU.GPU(); // Create GPU object.
    
    // Next, create kernels:
    // mod2() - Takes in array of complex numbers (one array for each component), and returns the modulus squared of each complex number.
    this.#mod2 = this.#gpu.createKernel(function(real, imag) {
      return real[this.thread.x] * real[this.thread.x] + imag[this.thread.x] * imag[this.thread.x];
    }).setOutput([size]);  // One value for each cell.

    // renderFrame() - Takes in arrays of source strengths, a gain value, an alpha value lookup array, and a max index for that array.
    // Returns an array of RGBA values for the sources.
    this.#renderFrame = this.#gpu.createKernel(function(source, gain, maxindex, lookup) {
      // R, G, and B values should all be 255 - only every 4th value (Alpha) needs to be calculated.
      if(this.thread.x % 4 != 3) {
        return 255;
      }

      // Calculate relative co-ordinate for non-RGBA array.
      const pos = Math.floor(this.thread.x / 4);

      // Calculate index, cap at "maxIndex", then return corresponding value in "lookup".
      const index = (Math.min(Math.trunc(source[pos] * gain), maxindex));
      return lookup[index];
    }).setOutput([size * 4]); // RGBA channel for each cell = number of cells * 4 channels.
  }
  process(real, imag) {
    const source = this.#mod2(real, imag); // Call mod2 on values to get source strengths.
    this.#max = Math.max(...source); // Grab max source strength. THIS WILL BE INEFFICIENT WITH TEXTURES, REWRITE?!
    return Uint8ClampedArray.from(this.#renderFrame(source, this.#gain, this.#maxindex, this.#lookup)); // Render and return frame.
  }
  resetGain() {
    this.#gain = this.#maxindex / this.#max;
    this.#max = 0;
  }
}

// Partial port of class from QuantumData.java in original
class GameRender {
  data; //int[]
  qd; // QuantumData
  colourmap; // AmpColourMap (for now)
  #image; // ImageData
  width() {return this.qd.width;}
  height() {return this.qd.height;}
  constructor(source) {
    this.qd = source;
    this.#image = new ImageData(new Uint8ClampedArray(this.width() * this.height() * 4), this.width()); // ImageData replaced BufferedImage
    this.data = new Uint8ClampedArray(this.width() * this.height() * 4); // Stores RGBA values used to create ImageData.
    this.colourmap = new AmpColourMap(this.width() * this.height()); // Supply size of array to colour map so kernel output can be set.
  }
  update() {
    this.data = this.colourmap.process(this.qd.getReal(), this.qd.getImag()); // Push complex numbers to kernel and get RGBA values back.
    this.colourmap.resetGain();
    this.#image = new ImageData(this.data, this.width()); // Create new ImageData object using values retrieved from kernel.
  }
  getImage() {
    return this.#image;
  }
}

// Partial port of class from LevelManager.java in original
class LevelManger {
  qd; // QuantumData
  gr; // GameRender
  controlstate; // ControlState
  #quantumframetimenanos; // long
  #scale; // float
  constructor() {
    this.controlstate = new ControlState();
  }
  init(scale, dt, maxtilt, thousanditertimesecs) {
    this.#scale = scale;
    this.qd = new QuantumData(Math.trunc(canvas.width/scale), Math.trunc(canvas.height/scale), this.controlstate);
    this.qd.setDeltaT(dt);
    this.qd.setMaxTilt(maxtilt);
    this.gr = new GameRender(this.qd);
    this.#quantumframetimenanos = Math.trunc(thousanditertimesecs/1000*1000000000); //Math.trunc closest to casting to long.
  }
  addGaussian(x, y, sigma, px, py, a) {
    this.qd.addGaussian(Math.trunc(x/this.#scale), Math.trunc(y/this.#scale), sigma, Math.trunc(px/this.#scale), Math.trunc(py/this.#scale), a);
  }
  getQD() {
    return this.qd;
  }
  updateGraphics() {
    if(this.controlstate.resetRequested()) {
      this.resetInitialState();
    }
    this.gr.update();
    // Following two lines replace lc.repaint in the original: clear canvas and then draw data in gr.#image to it.
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    createImageBitmap(this.gr.getImage()).then(renderer => ctx.drawImage(renderer, 0, 0, canvas.width, canvas.height)); // Creates new bitmap image using imageData and then scales it. This code was taken from Stack Overflow post made by user "Kaiido" on 18-07-2018. Accessed 27-02-2025. https://stackoverflow.com/questions/51387989/change-image-size-with-ctx-putimagedata
  }
  resetInitialState() {
    this.qd.resetInitialState();
  }
  quantumFrameTimeNanos() {
    return this.#quantumframetimenanos;
  }
}

//Basic priority queue in JavaScript
//This code was taken from Stack Overflow post made by user "gyer" on 21-03-2017
//Accessed 04-02-2025
//https://stackoverflow.com/questions/42919469/efficient-way-to-implement-priority-queue-in-javascript
{
  const top = 0;
  const parent = i => ((i + 1) >>> 1) - 1;
  const left = i => (i << 1) + 1;
  const right = i => (i + 1) << 1;
  
  class PriorityQueue {
    constructor(comparator = (a, b) => a > b) {
      this._heap = [];
      this._comparator = comparator;
    }
    size() {
      return this._heap.length;
    }
    isEmpty() {
      return this.size() == 0;
    }
    peek() {
      return this._heap[top];
    }
    push(...values) {
      values.forEach(value => {
        this._heap.push(value);
        this._siftUp();
      });
      return this.size();
    }
    pop() {
      const poppedValue = this.peek();
      const bottom = this.size() - 1;
      if (bottom > top) {
        this._swap(top, bottom);
      }
      this._heap.pop();
      this._siftDown();
      return poppedValue;
    }
    replace(value) {
      const replacedValue = this.peek();
      this._heap[top] = value;
      this._siftDown();
      return replacedValue;
    }
    _greater(i, j) {
      return this._comparator(this._heap[i], this._heap[j]);
    }
    _swap(i, j) {
      [this._heap[i], this._heap[j]] = [this._heap[j], this._heap[i]];
    }
    _siftUp() {
      let node = this.size() - 1;
      while (node > top && this._greater(node, parent(node))) {
        this._swap(node, parent(node));
        node = parent(node);
      }
    }
    _siftDown() {
      let node = top;
      while (
        (left(node) < this.size() && this._greater(left(node), node)) ||
        (right(node) < this.size() && this._greater(right(node), node))
      ) {
        let maxChild = (right(node) < this.size() && this._greater(right(node), left(node))) ? right(node) : left(node);
        this._swap(node, maxChild);
        node = maxChild;
      }
    }
  }
  window.PriorityQueue=PriorityQueue;
}
// End of referenced code

// Partial port of class from LevelManager.java in original
// No idea about the inheritance works in the original. For now I don't need any of the key or action listening functions so I will leave as is for now.
class ControlState {
  xslope; yslope; // float
  cursordisabled; // boolean
  speed = 0.08; // float
  constructor() {
    this.xslope = 0;
    this.yslope = 0;
    this.cursordisabled = false;
  }
  resetRequested() {
    // "only called once so can return current value and reset"
    if(this.#reset) {
      this.#reset = false;
      return true;
    } else {
      return false;
    }
  }
  #up = false; #down = false; #left = false; #right = false; #reset = false;
  #queueReset() {
    this.#reset = true;
  }
  set(key, value) {
    switch(key) {
      case "ArrowRight":
        this.#right = value;
        break;
      case "ArrowLeft":
        this.#left = value;
        break;
      case "ArrowUp":
        this.#up = value;
        break;
      case "ArrowDown":
        this.#down = value;
        break;
      case "KeyR":
        this.#queueReset();
        break;
    }
  }
  #getTargetXSlope() {
    if(this.#left && this.#right) {
      return 0;
    } else if(this.#left) {
      return -1;
    } else if (this.#right){
      return 1;
    } else {
      return 0;
    }
  }
  #getTargetYSlope() {
    if(this.#up && this.#down) {
      return 0;
    } else if(this.#up) {
      return -1;
    } else if (this.#down){
      return 1;
    } else {
      return 0;
    }
  }
  getXSlope() {
    if(!this.cursordisabled) {
      return this.xslope;
    } else {
      return 0
    }
  }
  getYSlope() {
    if(!this.cursordisabled) {
      return this.yslope;
    } else {
      return 0
    }
  }
  #getNewSlope(oldslope, target) {
    if(Math.abs(oldslope-target)<this.speed) {
      return target;
    } else {
      const movedir = Math.sign(target-oldslope); // Math.sign() in JavaScript and Math.signum() in Java are equivalent, I checked.
      let newslope = oldslope+this.speed*movedir;
      if(newslope<-1) {newslope=-1;}
      if(newslope>1) {newslope=1;}
      return newslope;
    }
  }
  step() {
    this.xslope = this.#getNewSlope(this.xslope,this.#getTargetXSlope());
    this.yslope = this.#getNewSlope(this.yslope,this.#getTargetYSlope());
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------------


// MAIN METHOD -----------------------------------------------------------------------------------------------------------------------------------
// Retrieving data from HTML.
const canvas = document.getElementById("game"); // Grab canvas from HTML
const ctx = canvas.getContext("2d"); // Create 2D context.

// Setting up objects.
const manager = new LevelManger(); // Create new LevelManager (top level object at the moment)

manager.init(2.5, 0.1, 10, 5/1.5); // manager.init() will have been called elsewhere by the time UpdateTask.run() is executed, so do it now. scale = 2.5, dt = 0.1, maxtilt = 10, thousanditertimesecs = 5/1.5 - all values taken from toofast.xml
manager.addGaussian(153, 263, 10, 0, 0, 1); // Also no point running a simulation if nothing to simulate, so add a gaussian. Values again taken from toofast.xml.

// Adding and binding key listeners.
document.addEventListener('keydown', keyDownEvent);
document.addEventListener('keyup', keyUpEvent);

function keyDownEvent(e) {
  manager.getQD().getCS().set(e.code, true);
}

function keyUpEvent(e) {
  manager.getQD().getCS().set(e.code, false);
}

// Begin implementation of UpdateTask.run()
const gfxframetime = 33000000; // "30 fps"
const quantum_frames_per_gfx_frame = gfxframetime / manager.quantumFrameTimeNanos();
let lastframetime = nanoTime();
let quantumframes_this_frame = 0;

requestAnimationFrame(graphicsLoop); // Start recursive calling of requestAnimationFrame, using graphicsLoop() as the function.