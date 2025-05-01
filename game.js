// Code throughout adapts, translates, and expands upon original Quantum Marble Maze by Dr. Crispin Cooper
// Taken from "fiftysevendegreesofrad" on GitHub, v1.0 formalized 17-05-2024
// Accessed 25-01-2025
// https://github.com/fiftysevendegreesofrad/quantum

"use strict"; // Directive to indicate code should be executed in strict mode.

//FUNCTION DECLARATIONS --------------------------------------------------------------------------------------------------------------------------

// graphicsLoop - A partial port of the run() method of the UpdateTask runnable object in QMM1.
// This function is exclusively called by the browser after earlier being passed as a parameter to a call of requestAnimationFrame().
// After requestAnimationFrame() is called, when the browser next decides it is time to render a graphics frame (30 fps if tab is in focus),
// it first calls the function that was passed to requestAnimationFrame().
// In this case, the function contains a loop which renders quantum simulation frames until A) a target number of quantum frames are reached (30),
// B) a time limit is reached (1/30th of a second), or C) the next level has been queued.
// In the case of A and B a new graphics frame is rendered and then requestAnimationFrame() is called again, restarting the progress.
// In C's case getNextLevel() is called, and after the next level is loaded requestAnimationFrame() is called with graphicsLoop() to restart the
// simulation.
async function graphicsLoop() {
  let nextLevel = false;                                            // Need to know when loop ends if this is because a new frame has been rendered or we need to load the next level.
  while(true) {                                                     // Until break called by frame time limit or level ending:
    quantumframes_this_frame++;                                       // Increment the quantum frame count for this graphics frame stored in global scope (needed to get around requestAnimationFrame).
    if(quantumframes_this_frame < quantum_frames_per_gfx_frame) {     // If we haven't hit our target number of qFrames yet:
      totalQFrames++;                                                   // Increment the number of total qFrames for testing
      //console.time("qFrame " + totalQFrames);                         // Start timing
      manager.getQD().step();                                           // Advance the simulation one step.
      //console.timeEnd("qFrame " + totalQFrames);                      // Finish timing and display.
    } else {                                                          // Else if we have hit our target:
      const timesincelastframe = performance.now() - lastframetime;     // Calculate how much time has passed since last graphics frame.
      const sleeptime = gfxframetime - timesincelastframe;              // Calculate the amount of time until the end of the frame
      if(sleeptime > 0) {                                               // If this time is greater than zero:
        await new Promise(r => setTimeout(r, sleeptime));                 // Sleep for the requisite amount of time.
      }
    }
    const currenttime = performance.now();                            // Take a time reading
    if(currenttime - lastframetime > gfxframetime) {                  // Calculate the amount of time since last frame. If it is greater than the target time:        
      quantumframes_this_frame = 0;                                     // Reset the qFrame count to zero.
      lastframetime = currenttime;                                      // Note the new time since last frame.
      totalGfxFrames++;                                                 // Increment the graphics frame counter for testing.
      //console.time("gfxFrame " + totalGfxFrames);                     // Start timing
      manager.updateGraphics();                                         // Render and display a new graphics frame.
      //console.timeEnd("gfxFrame " + totalGfxFrames);                  // Finish timing and display.
      break;                                                            // Need to pass control back to browser to await next graphics frame, so break.
    }
    if(manager.shouldTerminate()) {                                   // If terminate level flag has been raised:
      nextLevel = true;                                                 // Amend boolean flag so we remember to call getNextLevel()
      break                                                             // Break so we can call it.
    }
  }
  if(nextLevel) {                                                   // If we broke loop to load new level:
    getNextLevel();                                                   // Do so
  } else {                                                          // Otherwise:
    requestAnimationFrame(graphicsLoop);                              // Schedule the next animation frame and pass control back to browser.
  }
}

function loadLevel(baseurl, levelfile, lm) {
  try {
    fetch(baseurl + levelfile + ".xml")
      .then((response) => response.text())
      .then((xmlString) => {
          const parser = new DOMParser();
          const dom = parser.parseFromString(xmlString, "text/xml");
          const root = dom.documentElement;
          lm.init(root.getAttribute("scale"),
              Number(root.getAttribute("dt")),
              root.getAttribute("maxtilt"),
              root.getAttribute("qft")/1.5);
          
          const nodes = root.childNodes;

          // "get mask first" (done differently here because async)
          const img = document.createElement('img');
          img.crossOrigin = "Anonymous";
          img.src = baseurl + dom.getElementsByTagName("mask")[0].childNodes[0].nodeValue;
          img.decode().then(() => {
            ctx.drawImage(img, 0, 0);
            lm.mask = ctx.getImageData(0, 0, img.width, img.height);

            for(let i = 0; i < nodes.length; i++) {
              const node = nodes[i];
              //element?
              switch(node.nodeName) {
                case "synth":
                  //dummy
                  break;
                case "audioloop":
                  //dummy
                  break;
                case "background":
                  let bg = new Image();
                  bg.src = baseurl + node.childNodes[0].nodeValue;
                  lm.setBackground(bg);
                  break;
                case "text":
                  //dummy
                  break;
                case "potentialplane":
                  lm.addPotentialPlane(node.getAttribute("tl"),
                      node.getAttribute("tr"),
                      node.getAttribute("bl"),
                      node.getAttribute("br"),
                      node.getAttribute("mask"));
                  break;
                case "potentialwell":
                  //dummy
                  break;
                case "potentialcone":
                  //dummy
                  break;
                case "gaussian":
                  lm.addGaussian(node.getAttribute("x"),
                      node.getAttribute("y"),
                      node.getAttribute("sigma"),
                      node.getAttribute("px"),
                      node.getAttribute("py"),
                      node.getAttribute("a"));
                  break;
                case "delta":
                  //dummy
                  break;
                case "goal":
                  lm.addGoal(node.getAttribute("mask"), node.getAttribute("target"));
                  break;
                case "reward":
                  //dummy
                  break;
                case "walls":
                  lm.setWalls(node.getAttribute("mask"));
                  break;
                case "sink":
                  //dummy
                  break;
                case "collapse":
                  lm.addCollapse(node.getAttribute("mask"), node.getAttribute("target"), node.getAttribute("sigma"));
                  break;
                case "trap":
                  //dummy
                  break;
                case "steepeningvalley":
                  //dummy
                  break;
              }
            }
          
            quantum_frames_per_gfx_frame = gfxframetime / lm.quantumFrameTimeNanos();
            lastframetime = performance.now();
            quantumframes_this_frame = 0;
            totalQFrames = 0;
            totalGfxFrames = 0;
            requestAnimationFrame(graphicsLoop);
          })
      });
  } catch(error) {
    //dummy
  }
}

function getNextLevel() {
  manager = new LevelManger("./levels/", levelNames.shift());
}
//------------------------------------------------------------------------------------------------------------------------------------------------


//CLASS DECLARATIONS -----------------------------------------------------------------------------------------------------------------------------

// # prefix is used to denote private properties/methods in JavaScript
// this. prefix required for all property references in JavaScript, even within same class.
// Can't really figure out how to do final without const, which doesn't work for class properties? If there is some easy way to emulate this, let me know - but I don't think it is actually necessary or even desirable.

// Partial port of class from QuantumData.java in original
class QuantumData {
  #levelDesignPotential; #pot_cache;
  sink_mult;
  #real; #imag; #init_real; #init_imag;
  #walls; #sink;
  #counters;
  #delta_t; #maxtilt; // float
  running = false; // boolean

  width; height; // final int
  #controlstate; // ControlState

  #gpu;
  #stepComponent;
  #toTexture;
  #resetPotentialCacheKernel;
  #loadWalls;

  constructor(width, height, ks) {
    this.#controlstate = ks;
    this.width = width;
    this.height = height;
    this.#real = new Float32Array(width * height).fill(0);
    this.#imag = new Float32Array(width * height).fill(0);
    this.#init_real = new Float32Array(width * height).fill(0);
    this.#init_imag = new Float32Array(width * height).fill(0);
    this.#walls = new Array(width * height).fill(false);
    this.#sink = new Array(width * height).fill(false);
    this.sink_mult = new Float32Array(width * height).fill(0);
    this.#levelDesignPotential = new Float32Array(width * height).fill(0);
    this.#pot_cache = new Float32Array(width * height).fill(0);
    this.#counters = [];
    this.#gpu = new GPU.GPU();
    this.#setupKernels();
  }
  getReal() {
    return this.#real;
  }
  getImag() {
    return this.#imag;
  }
  getGPU() {
    return this.#gpu;
  }
  #setupKernels() {
    this.#stepComponent = this.#gpu.createKernel(function(w, h, delta_t, update, ref, walls, sink_mult, pot_cache, sign) {
      // Calculate x and y co-ordinates from index in 1D array.
      const x = this.thread.x % w;
      const y = Math.floor(this.thread.x / w);
  
      // If the co-ordinates are on the border or inside a wall they should not be simulated - return 0.
      if(x < 1 || x >= w - 1 || y < 1 || y >= h - 1 || walls[x+w*y] == 1) {
          return 0;
      }
  
      // Else, return the updated value of component using the same formula as QMM1 (with the addition of "sign" for kernel reuse).
      return sink_mult[x+w*y] * (update[x+w*y] + sign * delta_t * (-0.5 * (ref[x+w*(y-1)]+ref[x+w*(y+1)]+ref[(x-1)+w*y]+ref[(x+1)+w*y]-4*ref[x+w*y]) + pot_cache[x+w*y]*ref[x+w*y]));
    }).setOutput([this.width * this.height]).setPipeline(true).setImmutable(true);

    this.#toTexture = this.#gpu.createKernel(function(values){
      return values[this.thread.x];
    }).setOutput([this.width * this.height]).setPipeline(true).setImmutable(true);

    this.#loadWalls = this.#gpu.createKernel(function(values){
      return values[this.thread.x];
    }).setOutput([this.width * this.height]).setPipeline(true).setImmutable(true);

    this.#resetPotentialCacheKernel = this.#gpu.createKernel(function(w, h, newtopleft, x_pot_step, y_pot_step, levelDesignPotential) {
      // Calculate x and y co-ordinates from index in 1D array.
      const x = this.thread.x % w;
      const y = Math.floor(this.thread.x / w);

      // If the co-ordinates are on the border they should not be simulated - return 0.
      if(x < 1 || x >= w - 1 || y < 1 || y >= h - 1) {
          return 0;
      }

      // Else calculate now potential.
      return newtopleft + (x_pot_step * x) + (y_pot_step * y) + levelDesignPotential[x+w*y];
    }).setOutput([this.width * this.height]).setPipeline(true).setImmutable(true);
  }
  #saveInitialState() {
    this.#init_real = this.#real;
    this.#init_imag = this.#imag;
  }
  resetInitialState() {
    this.#real = this.#toTexture(this.#init_real);
    this.#imag = this.#toTexture(this.#init_imag);
  }
  texturize() {
    this.#real = this.#toTexture(this.#real);
    this.#imag = this.#toTexture(this.#imag);
  }
  clearWaveFunction() {
    this.#real = new Float32Array(this.width * this.height).fill(0);
    this.#imag = new Float32Array(this.width * this.height).fill(0);
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
        this.#real[x+this.width*y] = this.#real[x+this.width*y] + vr;
        this.#imag[x+this.width*y] = this.#imag[x+this.width*y] + vi;
      }
    }
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

    const oldCache = this.#pot_cache;
    this.#pot_cache = this.#resetPotentialCacheKernel(this.width, this.height, newtopleft, x_pot_step, y_pot_step, this.#levelDesignPotential);
    oldCache.delete();
  }
  #ensure_no_positive_potential() {
    let maxpot = Number.NEGATIVE_INFINITY;
		for(let x=0;x<this.width;x++) {
      for(let y=0;y<this.height;y++) {
        const pot = this.#levelDesignPotential[x+this.width*y]
        if(pot>maxpot) {maxpot = pot};
      }
    }
		for(let x=0;x<this.width;x++) {
      for(let y=0;y<this.height;y++) {
        this.#levelDesignPotential.set(x, y,
            this.#levelDesignPotential[x+this.width*y]-maxpot);
      }
    }
  }
  #add_walls() {
    for(let x=1;x<this.width-1;x++) {
      for(let y=1;y<this.height-1;y++) {
        if(this.#walls[x+this.width*y]) {
          this.#real[x+this.width*y] = 0;
          this.#imag[x+this.width*y] = 0;
        }
      }      
    }
  }
  step() {
    if(!this.running) {
      this.running = true;
      this.#setupSinkMult();
      this.#add_walls(); // "must be done before saveInitialState"
      this.#saveInitialState();
      this.#ensure_no_positive_potential();
      this.#walls = this.#loadWalls(this.#walls);
      this.#real = this.#toTexture(this.#real);
      this.#imag = this.#toTexture(this.#imag);
      this.#pot_cache = this.#toTexture(this.#pot_cache);
    }
    this.#controlstate.step();
    this.#reset_potential_cache();

    const oldReal = this.#real;
    this.#real = this.#stepComponent(this.width, this.height, this.#delta_t, this.#real, this.#imag, this.#walls, this.sink_mult, this.#pot_cache, 1);
    oldReal.delete();
    const oldImag = this.#imag;
    this.#imag = this.#stepComponent(this.width, this.height, this.#delta_t, this.#imag, this.#real, this.#walls, this.sink_mult, this.#pot_cache, -1);
    oldImag.delete();
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
        this.sink_mult[x+this.width*y] = Number.POSITIVE_INFINITY;
        if(!this.#sink[x+this.width*y] && !this.#walls[x+this.width*y]) {
          queue.push(new Pixel(x, y, 0))
        }
      }      
    }
    while(!queue.isEmpty()) {
      const p = queue.pop(); // pop() equivalent to Java poll()
      if(this.sink_mult[p.x+this.width*p.y] > p.d) {
        this.sink_mult[p.x+this.width*p.y] = p.d;
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
        const dist = this.sink_mult[x+this.width*y];
        this.sink_mult[x+this.width*y] = Math.exp(-Math.pow(dist/2,2)*suddenness);
      }      
    }
  }
  addPotential(x, y, p) {
    this.#levelDesignPotential[x + this.width * y] = p + this.#levelDesignPotential[x + this.width * y];
  }
  setWalls(submask) {
    this.#walls = submask;
  }
  addCounter(submask) {
    this.#counters.push(submask);
  }
  getCounterScores(counters2) {
    if(this.#counters.length != counters2.length) {
      console.log("Error in getCounterScores(): this.#counters.length != counters2.length")
    }
    const tempReal = this.#real.toArray();
    const tempImag = this.#real.toArray();

    //"find total probability for renormalization"
    let totalprob = 0;
    for (let i = 0; i < this.width * this.height; i++) {
      totalprob += tempReal[i] * tempReal[i] + tempImag[i] * tempImag[i];
    }

    for (let counterid = 0; counterid < this.#counters.length; counterid++) {
      let score = 0;
      const currentcounter = this.#counters[counterid];
      for (let i = 0; i < this.width * this.height; i++) {
        if(currentcounter[i]) {
          score += tempReal[i] * tempReal[i] + tempImag[i] * tempImag[i];
        }
      }
      counters2[counterid].setValue(Math.floor(score/totalprob * 100) || 0); // converts "falsy" values apparently, avoids divide by 0 after collapse
    }

  }
}

class PosGoalCounter {
  target;
  value = 0;
  textLocation;
  lm;
  constructor(lm, target, textLocation) {
    this.target = target;
    this.textLocation = textLocation;
    this.lm = lm;
  }
  setValue(v) {this.value = v;}
  getTextLocation() {return this.textLocation;}
  getColour() {return "#ffffff";}
  getText() {return "GOAL+\n" + this.value + "/" + this.target + "%"}
  check() {
    if(this.value < this.target) this.lm.reportUnsatisfiedGoal();
  }
  reset() {this.value = 0;}
}

class CollapseCounter {
  target;
  value = 0;
  textLocation;
  lm;
  sigma;
  active = true;
  constructor(lm, target, textLocation, sigma) {
    this.target = target;
    this.textLocation = textLocation;
    this.sigma = sigma;
    this.lm = lm;
  }
  setValue(v) {this.value = v;}
  getTextLocation() {return this.textLocation;}
  getColour() {return "#ffff00"}
  getText() {
    if(this.active) {
      return "COLLAPSE+\n" + this.value + "/" + this.target + "%";
    } else {
      return "";
    }
  }
  check() {
    if(this.active && this.value >= this.target) {
      this.lm.clearWaveFunction();
      this.lm.addGaussianQUnits(this.textLocation[0], this.textLocation[1], this.sigma, 0, 0, 1);
      this.active = false;
    }
  }
  reset() {this.value = 0;}
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
  #gpu;
  #mod2;
  #render;
  constructor(width, height, gpu) {
    this.#lookup = new Int32Array(this.#maxindex + 1); //Int32Array closest to int[] in Java.
    for (let i=0;i<this.#maxindex+1;i++) {
      this.#lookup[i] = Math.trunc(255*Math.pow(i/this.#maxindex,this.#gamma)) //Math.trunc() closest to how casting to int works in Java.
    }
    this.#setupGPU(width, height, gpu);
  }
  #setupGPU(width, height, gpu) {
    this.#gpu = gpu;
    this.#mod2 = this.#gpu.createKernel(function(real, imag) {
      return real[this.thread.x]*real[this.thread.x]+imag[this.thread.x]*imag[this.thread.x];
    }).setOutput([width * height]).setPipeline(true).setImmutable(true);

    this.#render = this.#gpu.createKernel(function(source, gain, maxindex, lookup) {
      if(this.thread.x % 4 != 3) {
        return 255;
      }
      const index = Math.min(Math.trunc(source[this.thread.x / 4]*gain), maxindex);
      return lookup[index];
    }).setOutput([width * height * 4]);
  }
  process(real, imag) {
    const source = this.#mod2(real, imag);
    this.#max = Math.max(...source.toArray()); // Slow*er* but not too bad, I checked.
    const result = this.#render(source, this.#gain, this.#maxindex, this.#lookup); // MAKE THESE CONSTANTS
    source.delete();
    return new Uint8ClampedArray(result);
  }
  resetGain() {
    this.#gain = this.#maxindex / this.#max;
    this.#max = 0;
  }
}

// Partial port of class from QuantumData.java in original
class GameRender {
  qd; // QuantumData
  colourmap; // AmpColourMap (for now)
  #image; // ImageData
  width() {return Math.trunc(this.qd.width);}
  height() {return Math.trunc(this.qd.height);}
  constructor(source) {
    this.qd = source;
    this.#image = new ImageData(new Uint8ClampedArray(this.width() * this.height() * 4), this.width()); //ImageData replaced BufferedImage
    this.colourmap = new AmpColourMap(this.width(), this.height(), this.qd.getGPU());
  }
  update() {
    // "it may seem perverse to calculate 'data' only to copy it into 'image'"
    // "rather than just calculate image.  but i profiled and it's faster."
    this.#image = new ImageData(this.colourmap.process(this.qd.getReal(), this.qd.getImag()), this.width());
    this.colourmap.resetGain();
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
  mask; background = new Image();
  counters = [];
  #allGoalsSatisfiedThisRound; #allGoalsSatisfiedThreadSafe = false;
  #quantumframetimemillis; // long
  #scale; // float
  constructor(basename, filename) {
    this.controlstate = new ControlState();
    loadLevel(basename, filename, this);
  }
  init(scale, dt, maxtilt, thousanditertimesecs) {
    this.#scale = scale;
    this.qd = new QuantumData(Math.trunc(canvas.width/scale), Math.trunc(canvas.height/scale), this.controlstate);
    this.qd.setDeltaT(dt);
    this.qd.setMaxTilt(maxtilt);
    this.gr = new GameRender(this.qd);
    this.#quantumframetimemillis = thousanditertimesecs;
  }
  addGaussian(x, y, sigma, px, py, a) {
    this.qd.addGaussian(Math.trunc(x/this.#scale), Math.trunc(y/this.#scale), sigma, Math.trunc(px/this.#scale), Math.trunc(py/this.#scale), a);
  }
  addGaussianQUnits(x, y, sigma, px, py, a) {
    this.qd.addGaussian(x, y, sigma, px, py, a);
    this.qd.texturize();
  }
  setBackground(background) {
    this.background = background;
  }
  addPotentialPlane(tl, tr, bl, br, mask) {
    const m = this.#getSubMask(mask);
    // "find extremes"
    let top = this.qd.height; let left = this.qd.width; let bottom = 0; let right = 0;
    for(let x = 0; x < this.qd.width; x++) {
      for(let y = 0; y < this.qd.height; y++) {
        if(m[x + y * this.qd.width]) {
          if(x < left) {left = x;}
          if(x > right) {right = x;}
          if(y < top) {top = y;}
          if(y > bottom) {bottom = y;}
        }
      }
    }
    // "add potential"
    let potwidth = right-left;
    if(potwidth == 0) {potwidth = 1;}
    let potheight = bottom-top;
    if(potheight == 0) {potheight = 1;}
    for(let x = left; x <= right; x++) {
      for(let y = top; y <= bottom; y++) {
        if(m[x + y * this.qd.width]) {
          const rx = (x - left) / potwidth;
          const ry = (y - top) / potheight;
          // "potx0, potxh = potential at x,0 and x,height"
          const potx0 = tl + (tr - tl) * rx;
          const potxh = bl + (br - bl) * rx;
          const p = potx0 + (potxh - potx0) * ry;
          this.qd.addPotential(x, y, p);
        }
      }
    }
  }
  setWalls(mask) {
    this.qd.setWalls(this.#getSubMask(mask));
  }
  #getSubMask(desired_colour) {
    const wanted = [Number("0x" + desired_colour.slice(1, 3)), Number("0x" + desired_colour.slice(3, 5)), Number("0x" + desired_colour.slice(5, 7))];
    let best_dist = Math.POSITIVE_INFINITY;
    let chosen_colour = -1;
    const submask = new Array(Math.trunc(this.mask.width/this.#scale) * Math.trunc(this.mask.height/this.#scale)).fill(false);
    for(let x = 0; x < this.qd.width; x++) {
      for(let y = 0; y < this.qd.height; y++) {
        const imgX = Math.trunc(x * this.#scale);
        const imgY = Math.trunc(y * this.#scale);
        const imgPos = (imgX + this.mask.width * imgY);
        const thisRed = this.mask.data[imgPos*4 + 0];
        const thisGreen = this.mask.data[imgPos*4 + 1];
        const thisBlue = this.mask.data[imgPos*4 + 2];
        if(thisRed == wanted[0] && thisGreen == wanted[1] && thisBlue == wanted[2]) {
          submask[x + this.qd.width * y] = true;
        }
      } 
    }
    return submask;
  }
  #getCOG(mask) {
    let xtot = 0; let ytot = 0; let n = 0;
    for(let x = 0; x < this.qd.width; x++) {
      for(let y = 0; y < this.qd.height; y++) {
        if(mask[x + this.qd.width * y]) {
          xtot += x;
          ytot += y;
          n++;
        }
      }      
    }
    return [xtot/n, ytot/n];
  }
  getQD() {
    return this.qd;
  }
  updateGraphics() {
    if(this.controlstate.resetRequested()) {
      this.resetInitialState();
    }
    this.gr.update();
    this.#checkCounters();
    // Following two lines replace lc.repaint in the original: clear canvas and then draw data in gr.#image to it.
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(this.background, 0, 0, canvas.width, canvas.height);
    createImageBitmap(this.gr.getImage()).then(renderer => ctx.drawImage(renderer, 0, 0, canvas.width, canvas.height)); // Creates new bitmap image using imageData and then scales it. This code was taken from Stack Overflow post made by user "Kaiido" on 18-07-2018. Accessed 27-02-2025. https://stackoverflow.com/questions/51387989/change-image-size-with-ctx-putimagedata
    for(let i = 0; i < this.counters.length; i++) {
      ctx.fillStyle = this.counters[i].getColour();
      const pos = this.counters[i].getTextLocation();
      ctx.fillText(this.counters[i].getText(), (pos[0] * this.#scale) - (ctx.measureText(this.counters[i].getText()).width / 2), pos[1] * this.#scale)
    }
  }
  #checkCounters() {
    this.qd.getCounterScores(this.counters);
    this.#allGoalsSatisfiedThisRound = true;
    for(let i = 0; i < this.counters.length; i++) {
      this.counters[i].check();
    }
    this.#allGoalsSatisfiedThreadSafe = this.#allGoalsSatisfiedThisRound;
  }
  shouldTerminate() {
    return this.#allGoalsSatisfiedThreadSafe || this.controlstate.menuRequested();
  }
  resetInitialState() {
    this.qd.resetInitialState();
    for(let i = 0; i < this.counters.length; i++) {
      this.counters[i].reset();
    }
  }
  addGoal(mask, target) {
    const submask = this.#getSubMask(mask);
    this.qd.addCounter(submask);
    const gc = new PosGoalCounter(this, target, this.#getCOG(submask));
    this.counters.push(gc);
  }
  addCollapse(mask, target, sigma) {
    const submask = this.#getSubMask(mask);
    this.qd.addCounter(submask);
    const cc = new CollapseCounter(this, target, this.#getCOG(submask), sigma);
    this.counters.push(cc);
  }
  clearWaveFunction() {
    this.qd.clearWaveFunction();
  }
  reportUnsatisfiedGoal() {
    this.#allGoalsSatisfiedThisRound = false;
  }
  quantumFrameTimeNanos() {
    return this.#quantumframetimemillis;
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
  menuRequested() {
    // "called multiple times"
		// "but once menu is requested during class lifespan"
		// "it can't be cancelled, so no need to reset"
		return this.#menu;
  }
  #up = false; #down = false; #left = false; #right = false; #reset = false; #menu = false;
  #queueReset() {
    this.#reset = true;
  }
  #queueMenu() {
    this.#menu = true;
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
      case "KeyM":
        // Hacky key-up stops multiple
        if(!value) {
          this.#queueMenu();
        }
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
ctx.font = "18px Arial";

//Firefox unsupported
if (navigator.userAgent.indexOf("Firefox") != -1) {
  ctx.fillText("Unfortunately, Firefox is not supported by this application. Please switch to a different browser.")
}

let manager = 0;

function keyDownEvent(e) {
  manager.controlstate.set(e.code, true);
}

function keyUpEvent(e) {
  manager.controlstate.set(e.code, false);
}

document.addEventListener('keydown', keyDownEvent);
document.addEventListener('keyup', keyUpEvent);

// Begin implementation of UpdateTask.run()
const gfxframetime = 33; // "30 fps"
let quantum_frames_per_gfx_frame = 0;
let lastframetime = 0;
let quantumframes_this_frame = 0;
let totalQFrames = 0;
let totalGfxFrames = 0;

const levelNames = ["c-bounce", "c-shape", "easysnake", "smorgasbord", "toofast", "medtraps", "hardtraps_introcollapse", "caketray0", "diffract", "caketray1", "annularmake", "annularwait", "annularwaitmove", "hardshell", "caketray2", "tunneling3"]
getNextLevel();