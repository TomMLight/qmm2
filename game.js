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

// Partial port of UpdateTask.run() from GameWindow.java in original, hacked to get requestAnimationFrame() to work.
function graphicsLoop() {
  // In the original, there is a main loop that breaks when the goal is sounded, and calls a repaint when a certain amount of time has passed.
  // I have struggled greatly to directly port this to JavaScript due to requestAnimationFrame, which as far as I can tell is the only way to achieve runtime animation.
  // If this is not the case, I will fix this asap.
  // Instead, this loop just executes until it is time to draw a new gfx frame, when it exits and recursively calls requestAnimationFrame again. 
  let loop = true;
  while(loop) {
    quantumframes_this_frame++;
    if(quantumframes_this_frame < quantum_frames_per_gfx_frame) {
      manager.getQD().step();
    } else {
      console.log("Should sleep!"); // This is never *ever* called in my tests, so I haven't actually implemented the sleeping yet!
    }
    const currenttime = nanoTime();
    if(currenttime - lastframetime > gfxframetime) {
      quantumframes_this_frame = 0;
      lastframetime = currenttime;
      manager.updateGraphics();
      loop = false; // New graphics frame needed, exit the loop and request new frame.
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

  #gpu // To save calling the constructor every time QuantumData.step() is called, a new gpu property is added to the class for future reference.

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
    this.#gpu = new GPU.GPU();
  }
  #saveInitialState() {
    this.#init_real.setEqualTo(this.#real);
    this.#init_imag.setEqualTo(this.#imag);
  }
  setDeltaT(dt) {this.#delta_t = dt;}
  setMaxTilt(mt) {this.#maxtilt = mt;}
  addGaussian(xc, yc, sigma, fx, fy, ascale) {
    const a = ascale * Math.pow(2*Math.PI*sigma*sigma,-0.25);
    const d = 4*sigma*sigma;
    const omegax = 2*Math.PI*fx; // "fixme this seems wrong" (Crispin has confirmed since this is not the case)
    const omegay = 2*Math.PI*fy; // "fixme this seems wrong" (Crispin has confirmed since this is not the case)
    for(let x=1;x<this.width-1;x++) {
      for(let y=1;y<this.height-1;y++) {
        const r2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);
        const vr = a * Math.exp(-r2/d)
            * Math.cos(omegax*x/this.width)
            * Math.cos(omegay*y/this.height);
        const vi = a * Math.exp(-r2/d)
            * Math.sin(omegax*x/this.width)
            * Math.sin(omegay*y/this.height);
        this.#real.set(x,y,this.#real.get(x,y)+vr);
        this.#imag.set(x,y,this.#imag.get(x,y)+vi);
      }
    }
  }
  get(x, y) {return new Complex(this.#real.get(x,y),this.#imag.get(x,y));}
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
    for(let x=1;x<this.width-1;x++) {
      for(let y=1;y<this.height-1;y++) {
        if(this.#walls.get(x,y)) {
          this.#real.set(x,y,0);
          this.#imag.set(x,y,0);
        }
      }      
    }
  }
  #unpack(arrayProperty) {
    const unpackedData = [];
    for(let x = 0; x < this.width; x++) {
      const thisRow = [];
      for(let y = 0; y < this.height; y++) {
        thisRow.push(arrayProperty.get(x, y));
      }
      unpackedData.push(thisRow);
    }
    return unpackedData
  }
  #repack(unpackedData, width, height) {
    const repackedData = new FloatArray2d(width, height);
    for(let x = 0; x < this.width; x++) {
      for(let y = 0; y < this.height; y++) {
        repackedData.set(x, y, unpackedData[y][x]); // Data appears to be transposed at some point during unpacking and operating on it. Switch it back now.
      }
    }
    return repackedData; 
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
    // "boundaries are never computed, hence left at 0"

    // No data (including class properties) can be read outside of kernel, so they must be unpacked into a normal 2D array now and passed in when the kernel is called.
    // In future the clear solution is to store properties in this form throughout instead of converting back and forth, but this will suffice for now as a test.
    let unpackedReal = this.#unpack(this.#real); // Let rather than const as will be replaced with result of computeReal.
    let unpackedImag = this.#unpack(this.#imag); // Let rather than const as will be replaced with result of computeImag.
    const unpackedWalls = this.#unpack(this.#walls);
    const unpackedSinkMult = this.#unpack(this.sink_mult);
    const unpackedPotCache = this.#unpack(this.#pot_cache);

    // Create a new "kernel" (GPU-executed function), stored in "computeReal", defined as follows:
    const computeReal = this.#gpu.createKernel(function(width, height, delta_t, real, imag, walls, sink_mult, pot_cache) {
      // Grab thread x and y now for easy reference later.
      const x = this.thread.x;
      const y = this.thread.y;

      // Clumsily skipping boundaries. Check if a better alternative exists later.
      // If not this can probably be combined with the wall check to make something more elegant.
      if(x < 1 || x >= width - 1 || y < 1 || y >= height - 1) {
        return real[x][y]; // If at a boundry no update should occur, so return existing value. (Maybe 0 would be better?)
      }

      // Contents of original double-nested loop. 
      if(walls[x][y] == 0) { // Booleans passed to/read by kernel as 1 or 0 instead of true or false.
        // Return the new value of real[x][y]. Calculation is perfect copy of original, just referencing normal 2D array instead of custom object.
        return sink_mult[x][y] * 
          (real[x][y] - delta_t * (-0.5 *
            (imag[x][y-1]+imag[x][y+1]+imag[x-1][y]+imag[x+1][y]-4*imag[x][y])
          + pot_cache[x][y] * imag[x][y]));
      } else {
        return real [x][y]; // If a wall is present no update should occur, so return existing value. (Maybe 0 would be better?)
      }
    }).setOutput([this.height, this.width]); // Set size of kernel output. Counter-intuitively, the y-dimension must be supplied before the x-dimension.

    // Call computeReal with unpacked parameters, and save the result in unpackedReal for repacking later (and use in computeImag).
    unpackedReal = computeReal(this.width, this.height, this.#delta_t, unpackedReal, unpackedImag, unpackedWalls, unpackedSinkMult, unpackedPotCache);
    this.#real.setEqualTo(this.#repack(unpackedReal, this.width, this.height)); // Unpack result into #real.

    /* Old non-GPU updating of real component commented out:
    
    for(let y=1;y<this.height-1;y++) {
      for(let x=1;x<this.width-1;x++) {
        if(!this.#walls.get(x,y)) { //Real component "fixed" at last!!!
          this.#real.set(x,y,
            this.sink_mult.get(x,y)*
          (this.#real.get(x,y) + this.#delta_t * (-0.5 *
            (this.#imag.get(x,y-1)+this.#imag.get(x,y+1)+this.#imag.get(x-1,y)+this.#imag.get(x+1,y)-4*this.#imag.get(x,y))
          + this.#pot_cache.get(x,y)*this.#imag.get(x,y))));
        }
      }      
    }
    */

    // "I have inlined del2, it does make it faster"
		// "Inlining could happen automatically with vm options -XX:FreqInlineSize=50 -XX:MaxInlineSize=50"
		// "But these are not universally supported or guaranteed not to change in future"
    for(let y=1;y<this.height-1;y++) {
      for(let x=1;x<this.width-1;x++) {
        if(!this.#walls.get(x,y)) { // The following calculation is EXACTLY equivalent to original after removing: whitespace, "this." from properties, "#" for private properties, and an "f" in the Java to denote a float.
          this.#imag.set(x,y,
              this.sink_mult.get(x,y)*
            (this.#imag.get(x,y) - this.#delta_t * (-0.5 *
              (this.#real.get(x,y-1)+this.#real.get(x,y+1)+this.#real.get(x-1,y)+this.#real.get(x+1,y)-4*this.#real.get(x,y))
            + this.#pot_cache.get(x,y)*this.#real.get(x,y))));
        }
      }      
    }
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
  // VERY temp
  getCS() {
    return this.#controlstate;
  }
  // ALSO VERY temp
  getWall(x, y) {
    return this.#walls.get(x, y);
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
  constructor() {
    this.#lookup = new Int32Array(this.#maxindex + 1); //Int32Array closest to int[] in Java.
    for (let i=0;i<this.#maxindex+1;i++) {
      this.#lookup[i] = Math.trunc(255*Math.pow(i/this.#maxindex,this.#gamma)) //Math.trunc() closest to how casting to int works in Java.
    }
  }
  process(c) {
    const source = c.mod2();
    if (source>this.#max) this.#max = source;
    let index = Math.trunc(source*this.#gain);
    if (index>this.#maxindex) index=this.#maxindex;
    const alpha = this.#lookup[index];
    const red = 255;
    const green = 255;
    const blue = 255;
    return[red, green, blue, alpha]; // In original bit shifting is done here to pack all 4 values into one int. In the original this is a requirement due to how setRGB, but this version has no such constraints. So, I've just returned it as a length 4 array to make my life easier. If this proves to break later I will change it.
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
  width() {return Math.trunc(this.qd.width);}
  height() {return Math.trunc(this.qd.height);}
  constructor(source) {
    this.qd = source;
    this.#image = new ImageData(new Uint8ClampedArray(this.width() * this.height() * 4), this.width()); //ImageData replaced BufferedImage
    this.data = new Array(this.width()*this.height()); // To avoid bit-shifting back and forth data is a 2D array instead of a 1D array of RGBA values packed into a single int. May attempt bit shifting later to quadruple check this isn't causing issues.
    this.colourmap = new AmpColourMap();
  }
  update() {
    const showpotential = false; // "for debugging"
    // "it may seem perverse to calculate 'data' only to copy it into 'image'"
    // "rather than just calculate image.  but i profiled and it's faster."
    for (let y = 0; y < this.qd.height; y++) {
      for (let x = 0; x < this.qd.width; x++) {
        const point = showpotential ? new Complex(0, 0) : this.qd.get(x,y) // No level potentials present at the moment, so imaginary component always 0.
        this.data[x + this.qd.width * y] = this.colourmap.process(point);
        if(this.qd.getWall(x, y)) {
          this.data[x + this.qd.width * y] = [100, 100, 100, 255]; // Shade walls in uniform grey. Functionality ok but feels very hacky maybe fix later.
        }
      }
    }
    
    this.colourmap.resetGain();
    this.#image = setRGB(0, 0, this.width(), this.height(), this.data, 0, this.width()); // In the original this is a method that updates the RGB values of 'image'. As imageData cannot be edited once created in JavaScript, this instead returns an ImageData object which replaces the existing one.
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
  constructor() {
    this.controlstate = new ControlState();
  }
  init(dt, maxtilt, thousanditertimesecs) {
    this.qd = new QuantumData(canvas.width, canvas.height, this.controlstate); // I don't really understand how scale works in the original? As I'm not yet parsing levels there's no real need, so just pass canvas dimensions instead.
    this.qd.setDeltaT(dt);
    this.qd.setMaxTilt(maxtilt);
    this.gr = new GameRender(this.qd);
    this.#quantumframetimenanos = Math.trunc(thousanditertimesecs/1000*1000000000); //Math.trunc closest to casting to long.
  }
  addGaussianQUnits(x, y, sigma, px, py, a) {
    this.qd.addGaussian(x, y, sigma, px, py, a);
  }
  getQD() {
    return this.qd;
  }
  updateGraphics() {
    this.gr.update();
    ctx.putImageData(this.gr.getImage(), 0, 0); // Draw data in gr.#image to canvas. Replaces lc.repaint in original.
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
  up = false; down = false; left = false; right = false;
  set(key, value) {
    switch(key) {
      case "ArrowRight":
        this.right = value;
        break;
      case "ArrowLeft":
        this.left = value;
        break;
      case "ArrowUp":
        this.up = value;
        break;
      case "ArrowDown":
        this.down = value;
        break;
    }
  }
  #getTargetXSlope() {
    if(this.left && this.right) {
      return 0;
    } else if(this.left) {
      return -1;
    } else if (this.right){
      return 1;
    } else {
      return 0;
    }
  }
  #getTargetYSlope() {
    if(this.up && this.down) {
      return 0;
    } else if(this.up) {
      return -1;
    } else if (this.down){
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

manager.init(0.1, 2.5, 1.6/1.5); // manager.init() will have been called elsewhere by the time UpdateTask.run() is executed, so do it now. dt = 0.1, maxtilt = 2.5, thousanditertimesecs = 1.6/1.5 - all values taken from c-bounce.xml
manager.addGaussianQUnits(200, 109, 2, 0, 0, 1); // Also no point running a simulation if nothing to simulate, so add a gaussian. Values taken from c-bounce.xml and then tweaked.

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