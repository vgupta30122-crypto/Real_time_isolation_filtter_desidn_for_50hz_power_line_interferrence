/**
 * ============================================================================
 * Real-Time Isolation Filter Design for 50Hz Power Line Interference
 * 3D Virtual Reality & Biomedical Signal Processing Laboratory
 * ============================================================================
 */

(function () {
  'use strict';

  // --- Experiment State & Configuration ---
  const state = {
    // Bio-signal settings
    signalType: 'ecg', // 'ecg' | 'eeg' | 'sine'
    heartRate: 72,      // BPM (for ECG)
    
    // 50Hz Power Line Interference (PLI)
    noiseAmp: 1.20,     // mV
    mainsFreq: 50.0,    // Hz (45.0 - 55.0 Hz)
    harmonic3rd: true,  // 150 Hz harmonic
    emFieldVisible: true,
    
    // Galvanic Isolation Stage
    isolationActive: true,
    groundLoopVoltage: 1.8, // Volts CM when non-isolated
    cmrrIsolated: 122,      // dB
    cmrrNonIsolated: 60,    // dB
    
    // Real-Time 50Hz Notch Filter Stage
    filterActive: true,
    filterType: 'iir',      // 'iir' | 'twint' | 'lms'
    notchFreq: 50.0,        // Hz (45.0 - 55.0 Hz)
    qFactor: 16.0,          // Quality factor (Q = f0 / BW)
    
    // Oscilloscope Display
    timebase: 100,          // ms / division (10 divisions total = 1000 ms)
    scaleCH: 0.5,           // mV / division
    showCH1: true,          // Clean Ref (Cyan)
    showCH2: true,          // Noisy Input (Yellow)
    showCH3: true,          // Filtered Output (Green)
    isPaused: false,
    
    // Audio Immersion
    audioActive: false,

    // Guided Experiment Mode (Step 1 to 8)
    guideStep: 1,
    
    // Live Computed Metrics
    metrics: {
      noiseRMS: 848.5,      // microvolts
      errorRMS: 42.1,       // microvolts
      snrIn: -4.2,          // dB
      snrOut: 27.8,         // dB
      reductionPct: 95.0,   // %
      attenuation: -46.8,   // dB
      residual50Hz: 12.4,   // microvolts
      fps: 60
    }
  };

  // --- Mathematical & Digital Signal Processing Engine ---
  const SAMPLING_RATE = 1000; // 1000 Hz sample rate (standard biomedical acquisition)
  const BUFFER_SIZE = 1000;   // 1 second circular ring buffer
  
  const rawBuffer = new Float32Array(BUFFER_SIZE);
  const filteredBuffer = new Float32Array(BUFFER_SIZE);
  const cleanRefBuffer = new Float32Array(BUFFER_SIZE);
  let bufferIndex = 0;
  let globalSampleCount = 0;

  // IIR Biquad Direct Form II Transposed State
  let iirState = { w1: 0, w2: 0 };
  let biquadCoeffs = { b0: 1, b1: 0, b2: 0, a1: 0, a2: 0 };

  // Adaptive LMS Filter State
  let lmsState = { wCos: 0, wSin: 0, mu: 0.006 };

  // --- Real-Time Cooley-Tukey Radix-2 FFT Engine ---
  class FastFourierTransform {
    constructor(size = 1024) {
      this.size = size;
      this.cosTable = new Float32Array(size / 2);
      this.sinTable = new Float32Array(size / 2);
      for (let i = 0; i < size / 2; i++) {
        this.cosTable[i] = Math.cos((2 * Math.PI * i) / size);
        this.sinTable[i] = Math.sin((2 * Math.PI * i) / size);
      }
      this.bitRev = new Uint32Array(size);
      let j = 0;
      for (let i = 0; i < size; i++) {
        this.bitRev[i] = j;
        let k = size >> 1;
        while (k <= j) {
          j -= k;
          k >>= 1;
        }
        j += k;
      }
      // Hanning window to prevent spectral leakage
      this.window = new Float32Array(size);
      for (let i = 0; i < size; i++) {
        this.window[i] = 0.5 * (1 - Math.cos((2 * Math.PI * i) / (size - 1)));
      }
      this.real = new Float32Array(size);
      this.imag = new Float32Array(size);
      this.mag = new Float32Array(size / 2);
    }

    compute(inputArray, offset = 0) {
      const N = this.size;
      for (let i = 0; i < N; i++) {
        const idx = (offset + i) % inputArray.length;
        const rev = this.bitRev[i];
        this.real[rev] = inputArray[idx] * this.window[i];
        this.imag[rev] = 0;
      }
      for (let halfSize = 1; halfSize < N; halfSize <<= 1) {
        const step = halfSize << 1;
        const tableStep = N / step;
        for (let i = 0; i < N; i += step) {
          for (let k = 0; k < halfSize; k++) {
            const tIdx = k * tableStep;
            const cos = this.cosTable[tIdx];
            const sin = this.sinTable[tIdx];
            const uR = this.real[i + k];
            const uI = this.imag[i + k];
            const vR = this.real[i + k + halfSize] * cos + this.imag[i + k + halfSize] * sin;
            const vI = this.imag[i + k + halfSize] * cos - this.real[i + k + halfSize] * sin;
            this.real[i + k] = uR + vR;
            this.imag[i + k] = uI + vI;
            this.real[i + k + halfSize] = uR - vR;
            this.imag[i + k + halfSize] = uI - vI;
          }
        }
      }
      const norm = 2.0 / (N * 0.5); // normalized for Hanning coherent gain
      for (let i = 0; i < N / 2; i++) {
        const r = this.real[i] * norm;
        const im = this.imag[i] * norm;
        const magLin = Math.sqrt(r * r + im * im);
        // Magnitude in dBV (clamped to -90 dB)
        this.mag[i] = 20 * Math.log10(Math.max(1e-5, magLin));
      }
      return this.mag;
    }
  }

  const fftEngine = new FastFourierTransform(1024);
  let latestRawSpectrum = new Float32Array(512);
  let latestFltSpectrum = new Float32Array(512);

  /**
   * Recalculates 2nd-Order Digital IIR Notch Biquad Filter Coefficients
   * Transfer Function: H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)
   */
  function updateBiquadCoefficients() {
    let f0 = state.notchFreq;
    let Q = state.qFactor;
    const fs = SAMPLING_RATE;

    // Gracefully handle invalid parameter combinations & Nyquist limit
    if (f0 <= 0 || f0 >= fs / 2) {
      console.warn('Invalid notch frequency clamped to 50 Hz:', f0);
      f0 = 50.0;
      state.notchFreq = 50.0;
    }
    if (Q <= 0.1) {
      Q = 0.5;
      state.qFactor = 0.5;
    }
    
    const w0 = (2 * Math.PI * f0) / fs;
    const deltaW = w0 / Q;
    const beta = Math.tan(deltaW / 2);
    
    const norm = 1 / (1 + beta);
    biquadCoeffs.b0 = norm;
    biquadCoeffs.b1 = -2 * Math.cos(w0) * norm;
    biquadCoeffs.b2 = norm;
    biquadCoeffs.a1 = -2 * Math.cos(w0) * norm;
    biquadCoeffs.a2 = (1 - beta) * norm;

    // Update real-time frequency response plot in controls dock
    renderFilterFrequencyResponse();
  }

  /**
   * Process a single sample through the 2nd-Order IIR Notch Filter
   */
  function processIIRSample(x) {
    const { b0, b1, b2, a1, a2 } = biquadCoeffs;
    const y = b0 * x + iirState.w1;
    iirState.w1 = b1 * x - a1 * y + iirState.w2;
    iirState.w2 = b2 * x - a2 * y;
    return y;
  }

  /**
   * Active Twin-T Analog Filter Emulation
   */
  function processTwinTSample(x, t) {
    const deltaF = Math.abs(state.mainsFreq - state.notchFreq);
    const bw = state.notchFreq / state.qFactor;
    const atten = 1 / Math.sqrt(1 + Math.pow(bw / Math.max(0.01, 2 * deltaF), 2));
    
    const cleanBio = generateBioSignal(t);
    const noiseComponent = x - cleanBio;
    return cleanBio + noiseComponent * atten * 0.02;
  }

  /**
   * Adaptive LMS Filter (Least Mean Squares 50Hz cancellation)
   */
  function processLMSSample(x, t) {
    const w0 = 2 * Math.PI * state.notchFreq * t;
    const xCos = Math.cos(w0);
    const xSin = Math.sin(w0);
    
    const estimatedNoise = lmsState.wCos * xCos + lmsState.wSin * xSin;
    const error = x - estimatedNoise;
    
    // Weight update
    lmsState.wCos += 2 * lmsState.mu * error * xCos;
    lmsState.wSin += 2 * lmsState.mu * error * xSin;
    
    return error;
  }

  /**
   * Generates physiological clean biomedical signal at time t
   */
  function generateBioSignal(t) {
    if (state.signalType === 'ecg') {
      // Physiologically realistic cardiac P-Q-R-S-T waveform using Gaussian wavelets
      const period = 60 / state.heartRate;
      const tau = (t % period) / period; // Normalized cardiac phase [0, 1]
      
      const gaussian = (amp, center, width) => 
        amp * Math.exp(-Math.pow((tau - center) / width, 2));

      // P wave (atrial depolarization)
      const p = gaussian(0.18, 0.18, 0.035);
      // Q wave (septal depolarization)
      const q = gaussian(-0.16, 0.35, 0.012);
      // R wave (ventricular depolarization peak)
      const r = gaussian(1.35, 0.39, 0.015);
      // S wave (late ventricular depolarization)
      const s = gaussian(-0.32, 0.43, 0.018);
      // T wave (ventricular repolarization)
      const tWave = gaussian(0.28, 0.65, 0.065);
      // Gentle baseline respiration wander (0.25 Hz)
      const wander = 0.05 * Math.sin(2 * Math.PI * 0.25 * t);

      return p + q + r + s + tWave + wander;
    } 
    else if (state.signalType === 'eeg') {
      // EEG Alpha Rhythm (10 Hz) with spontaneous spindle burst modulation + Beta (20 Hz)
      const alphaBurst = (0.7 + 0.3 * Math.sin(2 * Math.PI * 0.8 * t)) * Math.sin(2 * Math.PI * 10.2 * t);
      const beta = 0.18 * Math.sin(2 * Math.PI * 21.0 * t);
      const theta = 0.12 * Math.sin(2 * Math.PI * 5.5 * t);
      return (alphaBurst * 0.35 + beta + theta);
    } 
    else {
      // Pure 10Hz calibration test sine wave (1.0 mV peak)
      return 1.0 * Math.sin(2 * Math.PI * 10.0 * t);
    }
  }

  /**
   * Generates 50Hz Power Line Interference and Ground Loop hum
   */
  function generateInterference(t) {
    const fMains = state.mainsFreq;
    // Fundamental 50Hz hum
    let pli = state.noiseAmp * Math.sin(2 * Math.PI * fMains * t);
    
    // 3rd harmonic (150 Hz)
    if (state.harmonic3rd) {
      pli += (state.noiseAmp * 0.3) * Math.sin(2 * Math.PI * 3 * fMains * t + 0.8);
    }
    
    // Galvanic Isolation Stage:
    // When Isolation is OFF, a ground loop exists. Chassis common-mode mains voltage (e.g. 1.8V)
    // converts to differential noise via finite CMRR (60 dB -> 1000:1 reduction)
    if (!state.isolationActive) {
      const cmNoise = (state.groundLoopVoltage / Math.pow(10, state.cmrrNonIsolated / 20)) * 1000; // in mV
      pli += cmNoise * Math.sin(2 * Math.PI * fMains * t + 1.2);
      pli += 0.08 * (Math.random() - 0.5); // ground loop thermal noise
    } else {
      // Optocoupler isolation barrier breaks ground loop; CMRR > 120 dB -> negligible leakage
      const cmLeakage = (state.groundLoopVoltage / Math.pow(10, state.cmrrIsolated / 20)) * 1000;
      pli += cmLeakage * Math.sin(2 * Math.PI * fMains * t + 1.2);
    }

    return pli;
  }

  /**
   * Generates real-time sample block and updates metrics
   */
  function stepSignalProcessing(samplesToGenerate) {
    if (state.isPaused) return;

    for (let i = 0; i < samplesToGenerate; i++) {
      const t = (globalSampleCount++) / SAMPLING_RATE;
      
      const cleanSignal = generateBioSignal(t);
      const interference = generateInterference(t);
      const rawSignal = cleanSignal + interference;
      
      let filteredSignal = rawSignal;
      if (state.filterActive) {
        if (state.filterType === 'iir') {
          filteredSignal = processIIRSample(rawSignal);
        } else if (state.filterType === 'twint') {
          filteredSignal = processTwinTSample(rawSignal, t);
        } else if (state.filterType === 'lms') {
          filteredSignal = processLMSSample(rawSignal, t);
        }
      }

      // Store in circular history buffers
      rawBuffer[bufferIndex] = rawSignal;
      filteredBuffer[bufferIndex] = filteredSignal;
      cleanRefBuffer[bufferIndex] = cleanSignal;
      
      bufferIndex = (bufferIndex + 1) % BUFFER_SIZE;
    }

    // Periodically update FFT & Performance Metrics (every ~100 ms)
    if (globalSampleCount % 100 < samplesToGenerate) {
      computeMetrics();
    }
  }

  /**
   * Compute Signal-to-Noise Ratio (SNR), RMS Noise, and True FFT Attenuation
   */
  function computeMetrics() {
    let powerSignal = 0;
    let powerNoiseIn = 0;
    let powerErrorOut = 0;
    const n = BUFFER_SIZE;

    for (let i = 0; i < n; i++) {
      const s = cleanRefBuffer[i];
      const r = rawBuffer[i];
      const f = filteredBuffer[i];
      
      powerSignal += s * s;
      powerNoiseIn += Math.pow(r - s, 2);
      powerErrorOut += Math.pow(f - s, 2);
    }

    powerSignal /= n;
    powerNoiseIn = Math.max(1e-12, powerNoiseIn / n);
    powerErrorOut = Math.max(1e-12, powerErrorOut / n);

    const noiseRMS_uV = Math.sqrt(powerNoiseIn) * 1000;
    const errorRMS_uV = Math.sqrt(powerErrorOut) * 1000;

    state.metrics.noiseRMS = noiseRMS_uV.toFixed(1);
    state.metrics.errorRMS = errorRMS_uV.toFixed(1);

    // Compute true FFT spectrums
    latestRawSpectrum = fftEngine.compute(rawBuffer, bufferIndex);
    latestFltSpectrum = fftEngine.compute(filteredBuffer, bufferIndex);

    // Bin corresponding to notch frequency (N=1024, fs=1000 => bin = f * 1024 / 1000)
    const binNotch = Math.min(511, Math.max(0, Math.round((state.notchFreq / SAMPLING_RATE) * 1024)));
    const rawPeakDB = latestRawSpectrum[binNotch];
    const fltPeakDB = latestFltSpectrum[binNotch];
    
    // Attenuation around notch frequency
    const attenDB = fltPeakDB - rawPeakDB;
    state.metrics.attenuation = attenDB <= 0 ? attenDB.toFixed(1) : `+${attenDB.toFixed(1)}`;

    // SNR In & SNR Out
    state.metrics.snrIn = (10 * Math.log10(powerSignal / powerNoiseIn)).toFixed(1);
    state.metrics.snrOut = (10 * Math.log10(powerSignal / powerErrorOut)).toFixed(1);

    // Noise Reduction %
    const redPct = Math.max(0, Math.min(100, (1 - errorRMS_uV / Math.max(1, noiseRMS_uV)) * 100));
    state.metrics.reductionPct = redPct.toFixed(1);

    // Residual 50Hz voltage
    state.metrics.residual50Hz = (Math.pow(10, fltPeakDB / 20) * 1000).toFixed(1);

    // Update DOM HUD metrics
    updateHUDDOM();
  }

  /**
   * Renders Filter Frequency Response |H(f)| on dock canvas
   */
  function renderFilterFrequencyResponse() {
    const canvas = document.getElementById('canvas-filter-response');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const w = canvas.width;
    const h = canvas.height;

    ctx.fillStyle = '#050a14';
    ctx.fillRect(0, 0, w, h);

    // Reticle lines: 0 dB, -3 dB, -40 dB
    const y0dB = 18;
    const yMinus3dB = 32;
    const yBottom = h - 16;

    ctx.strokeStyle = 'rgba(0, 229, 255, 0.12)';
    ctx.lineWidth = 1;

    ctx.beginPath();
    ctx.moveTo(35, y0dB); ctx.lineTo(w - 10, y0dB);
    ctx.moveTo(35, yMinus3dB); ctx.lineTo(w - 10, yMinus3dB);
    ctx.moveTo(35, yBottom); ctx.lineTo(w - 10, yBottom);
    ctx.stroke();

    ctx.fillStyle = '#8b9bb4';
    ctx.font = '10px "JetBrains Mono", monospace';
    ctx.fillText('0dB', 5, y0dB + 3);
    ctx.fillText('-3dB', 5, yMinus3dB + 3);
    ctx.fillText('-40dB', 5, yBottom);

    const fMax = 150;
    const f0 = state.notchFreq;
    const Q = state.qFactor;
    const fs = SAMPLING_RATE;
    const deltaF = f0 / Q;

    // Response curve
    ctx.strokeStyle = state.filterActive ? '#00ff88' : '#64748b';
    ctx.lineWidth = 2;
    ctx.beginPath();

    const numPoints = 120;
    for (let i = 0; i < numPoints; i++) {
      const f = (i / (numPoints - 1)) * fMax;
      const x = 35 + (f / fMax) * (w - 45);

      let magDB = 0;
      if (state.filterActive) {
        const wRad = (2 * Math.PI * f) / fs;
        const cosW = Math.cos(wRad);
        const sinW = Math.sin(wRad);
        const cos2W = Math.cos(2 * wRad);
        const sin2W = Math.sin(2 * wRad);

        const numR = biquadCoeffs.b0 + biquadCoeffs.b1 * cosW + biquadCoeffs.b2 * cos2W;
        const numI = -biquadCoeffs.b1 * sinW - biquadCoeffs.b2 * sin2W;

        const denR = 1 + biquadCoeffs.a1 * cosW + biquadCoeffs.a2 * cos2W;
        const denI = -biquadCoeffs.a1 * sinW - biquadCoeffs.a2 * sin2W;

        const numMag2 = numR * numR + numI * numI;
        const denMag2 = denR * denR + denI * denI;
        const magLin = Math.sqrt(numMag2 / Math.max(1e-12, denMag2));
        magDB = 20 * Math.log10(Math.max(1e-3, magLin));
      }

      const clampedDB = Math.max(-50, Math.min(2, magDB));
      const y = y0dB + (-clampedDB / 50) * (yBottom - y0dB);

      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // Notch Center Marker
    const xNotch = 35 + (f0 / fMax) * (w - 45);
    ctx.strokeStyle = '#00e5ff';
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(xNotch, 6);
    ctx.lineTo(xNotch, h - 6);
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.fillStyle = '#00e5ff';
    ctx.font = 'bold 10px "JetBrains Mono"';
    ctx.fillText(`f₀=${f0.toFixed(1)}Hz`, Math.max(35, xNotch - 24), h - 3);

    // Update label
    const bwLabel = document.getElementById('val-filter-bw');
    if (bwLabel) {
      bwLabel.innerText = `BW: ${deltaF.toFixed(1)} Hz (Q=${Q.toFixed(1)})`;
    }
  }

  // --- Web Audio API 50Hz Hum Immersion ---
  let audioCtx = null;
  let osc50 = null;
  let osc150 = null;
  let masterGain = null;

  function initAudio() {
    if (audioCtx) return;
    try {
      const AudioContext = window.AudioContext || window.webkitAudioContext;
      audioCtx = new AudioContext();

      masterGain = audioCtx.createGain();
      masterGain.gain.setValueAtTime(0, audioCtx.currentTime);
      masterGain.connect(audioCtx.destination);

      osc50 = audioCtx.createOscillator();
      osc50.type = 'sine';
      osc50.frequency.setValueAtTime(state.mainsFreq, audioCtx.currentTime);

      osc150 = audioCtx.createOscillator();
      osc150.type = 'sine';
      osc150.frequency.setValueAtTime(state.mainsFreq * 3, audioCtx.currentTime);

      const gain150 = audioCtx.createGain();
      gain150.gain.setValueAtTime(0.25, audioCtx.currentTime);
      osc150.connect(gain150);
      gain150.connect(masterGain);

      osc50.connect(masterGain);

      osc50.start();
      osc150.start();
    } catch (e) {
      console.warn('Web Audio API not supported or blocked:', e);
    }
  }

  function updateAudioGain() {
    if (!audioCtx || !masterGain) return;
    if (!state.audioActive || state.isPaused) {
      masterGain.gain.setTargetAtTime(0, audioCtx.currentTime, 0.05);
      return;
    }

    let targetVol = 0.08 * (state.noiseAmp / 2.0);
    if (!state.isolationActive) targetVol *= 2.0;

    if (state.filterActive) {
      targetVol *= Math.pow(10, parseFloat(state.metrics.attenuation) / 20);
    }

    masterGain.gain.setTargetAtTime(Math.min(0.25, Math.max(0, targetVol)), audioCtx.currentTime, 0.05);
    osc50.frequency.setTargetAtTime(state.mainsFreq, audioCtx.currentTime, 0.05);
    osc150.frequency.setTargetAtTime(state.mainsFreq * 3, audioCtx.currentTime, 0.05);
  }

  // --- Three.js 3D Virtual Reality Lab Scene ---
  let scene, camera, renderer, controls;
  let dsoTexture, dsoCanvas, dsoCtx;
  let fftTexture, fftCanvas, fftCtx;
  let emRings = [];
  let signalPulses = [];
  let pcbLedIso, pcbLedFilter, pcbLedClip;
  let vrButtonElement = null;

  // Preset Camera View Coordinates
  const cameraPresets = {
    overview: { pos: new THREE.Vector3(0, 3.2, 5.0), target: new THREE.Vector3(0, 1.2, 0) },
    dso:      { pos: new THREE.Vector3(-0.9, 1.7, 1.7), target: new THREE.Vector3(-0.9, 1.55, 0.2) },
    spectrum: { pos: new THREE.Vector3(0.9, 1.7, 1.7), target: new THREE.Vector3(0.9, 1.55, 0.2) },
    circuit:  { pos: new THREE.Vector3(0.0, 2.1, 1.2), target: new THREE.Vector3(0.0, 1.15, 0.4) },
    patient:  { pos: new THREE.Vector3(-2.2, 1.8, 2.2), target: new THREE.Vector3(-1.8, 1.1, 0.2) }
  };

  let targetCamPos = cameraPresets.overview.pos.clone();
  let targetCamLook = cameraPresets.overview.target.clone();
  let isTransitioningCamera = false;

  function transitionCameraTo(position, target) {
    targetCamPos.copy(position);
    targetCamLook.copy(target);
    isTransitioningCamera = true;
  }

  function initThreeJS() {
    const container = document.getElementById('canvas-container');
    const width = container.clientWidth || window.innerWidth;
    const height = container.clientHeight || window.innerHeight;

    // 1. Scene
    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x070a13);
    scene.fog = new THREE.FogExp2(0x070a13, 0.06);

    // 2. Camera
    camera = new THREE.PerspectiveCamera(50, width / height, 0.1, 50);
    camera.position.copy(cameraPresets.overview.pos);

    // 3. Renderer with WebGL 2 & WebXR VR enabled
    renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false, powerPreference: 'high-performance' });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.shadowMap.enabled = true;
    renderer.shadowMap.type = THREE.PCFSoftShadowMap;
    renderer.xr.enabled = true; // WebXR support
    container.appendChild(renderer.domElement);

    // Prevent context menu on right click so right-click pan works without interruption
    const preventContextMenu = (e) => e.preventDefault();
    renderer.domElement.addEventListener('contextmenu', preventContextMenu);
    container.addEventListener('contextmenu', preventContextMenu);

    // Append WebXR VR Button
    const vrContainer = document.getElementById('vr-container');
    if (typeof VRButton !== 'undefined') {
      vrButtonElement = VRButton.createButton(renderer);
      vrButtonElement.className = 'vr-cta-btn';
      vrContainer.appendChild(vrButtonElement);
    }

    // 4. Orbit Controls (Rotation, Zoom, Pan)
    if (typeof THREE.OrbitControls !== 'undefined') {
      controls = new THREE.OrbitControls(camera, renderer.domElement);
      controls.enableDamping = true;
      controls.dampingFactor = 0.05;
      controls.enableRotate = true;
      controls.rotateSpeed = 0.9;
      controls.enableZoom = true;
      controls.zoomSpeed = 1.2;
      controls.enablePan = true;
      controls.panSpeed = 0.8;
      controls.screenSpacePanning = true; // Standard screen-space panning on right click
      controls.maxPolarAngle = Math.PI / 2 - 0.02; // Don't clip below floor
      controls.minDistance = 0.5;
      controls.maxDistance = 14.0;
      controls.target.copy(cameraPresets.overview.target);

      // Explicit mouse buttons configuration
      controls.mouseButtons = {
        LEFT: THREE.MOUSE.ROTATE,
        MIDDLE: THREE.MOUSE.DOLLY,
        RIGHT: THREE.MOUSE.PAN
      };
      controls.touches = {
        ONE: THREE.TOUCH.ROTATE,
        TWO: THREE.TOUCH.DOLLY_PAN
      };

      // Cancel camera preset animation immediately when user begins manual mouse drag / scroll
      controls.addEventListener('start', () => {
        isTransitioningCamera = false;
      });
    }

    // 5. Lighting
    setupLighting();

    // 6. Build Virtual Environment
    buildLabRoom();
    buildWorkbench();
    buildOscilloscope();
    buildSpectrumAnalyzer();
    buildFilterCircuitBoard();
    buildPatientPhantom();
    buildPowerLineConduit();

    // 7. Window Resize
    window.addEventListener('resize', onWindowResize);
  }

  function setupLighting() {
    const ambientLight = new THREE.AmbientLight(0x1a263d, 1.2);
    scene.add(ambientLight);

    const benchLight = new THREE.SpotLight(0xe0f2fe, 2.5);
    benchLight.position.set(0, 4.5, 1.2);
    benchLight.target.position.set(0, 1.2, 0.3);
    benchLight.angle = Math.PI / 3.5;
    benchLight.penumbra = 0.4;
    benchLight.castShadow = true;
    benchLight.shadow.mapSize.width = 2048;
    benchLight.shadow.mapSize.height = 2048;
    benchLight.shadow.camera.near = 0.5;
    benchLight.shadow.camera.far = 8.0;
    benchLight.shadow.bias = -0.0008;
    scene.add(benchLight);
    scene.add(benchLight.target);

    const cyanAccent = new THREE.PointLight(0x00e5ff, 1.5, 4.5);
    cyanAccent.position.set(-1.2, 2.0, 0.8);
    scene.add(cyanAccent);

    const greenAccent = new THREE.PointLight(0x00ff88, 1.2, 4.5);
    greenAccent.position.set(1.2, 2.0, 0.8);
    scene.add(greenAccent);
  }

  function buildLabRoom() {
    const floorGeo = new THREE.PlaneGeometry(16, 16);
    const floorMat = new THREE.MeshStandardMaterial({
      color: 0x0a101f,
      roughness: 0.35,
      metalness: 0.4
    });
    const floor = new THREE.Mesh(floorGeo, floorMat);
    floor.rotation.x = -Math.PI / 2;
    floor.position.y = 0;
    floor.receiveShadow = true;
    scene.add(floor);

    const gridHelper = new THREE.GridHelper(16, 32, 0x00e5ff, 0x182642);
    gridHelper.position.y = 0.005;
    scene.add(gridHelper);

    const wallGeo = new THREE.PlaneGeometry(16, 7);
    const wallMat = new THREE.MeshStandardMaterial({
      color: 0x0d1527,
      roughness: 0.6,
      metalness: 0.2
    });
    const backWall = new THREE.Mesh(wallGeo, wallMat);
    backWall.position.set(0, 3.5, -3.5);
    scene.add(backWall);

    const bannerCanvas = document.createElement('canvas');
    bannerCanvas.width = 1024;
    bannerCanvas.height = 256;
    const bctx = bannerCanvas.getContext('2d');
    bctx.fillStyle = '#101a30';
    bctx.fillRect(0, 0, 1024, 256);
    bctx.strokeStyle = '#00e5ff';
    bctx.lineWidth = 6;
    bctx.strokeRect(10, 10, 1004, 236);
    
    bctx.fillStyle = '#00e5ff';
    bctx.font = 'bold 36px "JetBrains Mono", sans-serif';
    bctx.fillText('BIOMEDICAL INSTRUMENTATION & CALIBRATION LAB', 50, 75);
    
    bctx.fillStyle = '#8b9bb4';
    bctx.font = '22px "Inter", sans-serif';
    bctx.fillText('STANDARDS: IEC 60601-1 (MEDICAL ELECTRICAL SAFETY) • GALVANIC ISOLATION > 120 dB', 50, 125);
    bctx.fillText('POWER GRID: 230V RMS @ 50.0 Hz • STRAY CAPACITIVE COUPLING NOISE MITIGATION', 50, 165);

    const bannerTex = new THREE.CanvasTexture(bannerCanvas);
    const bannerMesh = new THREE.Mesh(
      new THREE.PlaneGeometry(6.5, 1.6),
      new THREE.MeshBasicMaterial({ map: bannerTex })
    );
    bannerMesh.position.set(0, 3.8, -3.48);
    scene.add(bannerMesh);
  }

  function buildWorkbench() {
    const tableGroup = new THREE.Group();

    const topGeo = new THREE.BoxGeometry(4.2, 0.1, 1.8);
    const topMat = new THREE.MeshStandardMaterial({
      color: 0x162238,
      roughness: 0.3,
      metalness: 0.6
    });
    const tabletop = new THREE.Mesh(topGeo, topMat);
    tabletop.position.set(0, 1.15, 0.3);
    tabletop.castShadow = true;
    tabletop.receiveShadow = true;
    tableGroup.add(tabletop);

    const legGeo = new THREE.CylinderGeometry(0.045, 0.045, 1.15, 16);
    const legMat = new THREE.MeshStandardMaterial({
      color: 0x475569,
      roughness: 0.2,
      metalness: 0.8
    });

    const legPositions = [
      [-1.95, 0.575, -0.45],
      [ 1.95, 0.575, -0.45],
      [-1.95, 0.575,  1.05],
      [ 1.95, 0.575,  1.05]
    ];

    legPositions.forEach(([x, y, z]) => {
      const leg = new THREE.Mesh(legGeo, legMat);
      leg.position.set(x, y, z);
      leg.castShadow = true;
      tableGroup.add(leg);
    });

    const shelfGeo = new THREE.BoxGeometry(4.0, 0.06, 0.65);
    const shelfMat = new THREE.MeshStandardMaterial({
      color: 0x0f172a,
      roughness: 0.4,
      metalness: 0.5
    });
    const shelf = new THREE.Mesh(shelfGeo, shelfMat);
    shelf.position.set(0, 1.55, -0.2);
    shelf.castShadow = true;
    tableGroup.add(shelf);

    scene.add(tableGroup);
  }

  function buildOscilloscope() {
    const dsoGroup = new THREE.Group();
    dsoGroup.position.set(-1.0, 1.58, -0.15);

    const chassisGeo = new THREE.BoxGeometry(1.35, 0.75, 0.4);
    const chassisMat = new THREE.MeshStandardMaterial({
      color: 0x1e293b,
      roughness: 0.5,
      metalness: 0.3
    });
    const chassis = new THREE.Mesh(chassisGeo, chassisMat);
    chassis.castShadow = true;
    dsoGroup.add(chassis);

    const frontGeo = new THREE.PlaneGeometry(1.28, 0.68);
    const frontMat = new THREE.MeshStandardMaterial({
      color: 0x0f172a,
      roughness: 0.7
    });
    const front = new THREE.Mesh(frontGeo, frontMat);
    front.position.z = 0.201;
    dsoGroup.add(front);

    dsoCanvas = document.createElement('canvas');
    dsoCanvas.width = 1024;
    dsoCanvas.height = 640;
    dsoCtx = dsoCanvas.getContext('2d');

    dsoTexture = new THREE.CanvasTexture(dsoCanvas);
    dsoTexture.minFilter = THREE.LinearFilter;
    dsoTexture.magFilter = THREE.LinearFilter;

    const screenGeo = new THREE.PlaneGeometry(0.85, 0.58);
    const screenMat = new THREE.MeshBasicMaterial({ map: dsoTexture });
    const screenMesh = new THREE.Mesh(screenGeo, screenMat);
    screenMesh.position.set(-0.18, 0, 0.205);
    dsoGroup.add(screenMesh);

    const knobGeo = new THREE.CylinderGeometry(0.035, 0.035, 0.03, 16);
    const knobMat = new THREE.MeshStandardMaterial({ color: 0x64748b, metalness: 0.7, roughness: 0.3 });

    const knob1 = new THREE.Mesh(knobGeo, knobMat);
    knob1.rotation.x = Math.PI / 2;
    knob1.position.set(0.42, 0.15, 0.215);
    dsoGroup.add(knob1);

    const knob2 = new THREE.Mesh(knobGeo, knobMat);
    knob2.rotation.x = Math.PI / 2;
    knob2.position.set(0.42, -0.05, 0.215);
    dsoGroup.add(knob2);

    const bncGeo = new THREE.CylinderGeometry(0.02, 0.02, 0.035, 12);
    const bnc1 = new THREE.Mesh(bncGeo, new THREE.MeshStandardMaterial({ color: 0x00e5ff, metalness: 0.8 }));
    bnc1.rotation.x = Math.PI / 2;
    bnc1.position.set(0.35, -0.22, 0.215);
    dsoGroup.add(bnc1);

    const bnc2 = new THREE.Mesh(bncGeo, new THREE.MeshStandardMaterial({ color: 0xffb703, metalness: 0.8 }));
    bnc2.rotation.x = Math.PI / 2;
    bnc2.position.set(0.44, -0.22, 0.215);
    dsoGroup.add(bnc2);

    const bnc3 = new THREE.Mesh(bncGeo, new THREE.MeshStandardMaterial({ color: 0x00ff88, metalness: 0.8 }));
    bnc3.rotation.x = Math.PI / 2;
    bnc3.position.set(0.53, -0.22, 0.215);
    dsoGroup.add(bnc3);

    scene.add(dsoGroup);
  }

  function buildSpectrumAnalyzer() {
    const fftGroup = new THREE.Group();
    fftGroup.position.set(1.0, 1.58, -0.15);

    const chassisGeo = new THREE.BoxGeometry(1.35, 0.75, 0.4);
    const chassisMat = new THREE.MeshStandardMaterial({
      color: 0x1e293b,
      roughness: 0.5,
      metalness: 0.3
    });
    const chassis = new THREE.Mesh(chassisGeo, chassisMat);
    chassis.castShadow = true;
    fftGroup.add(chassis);

    fftCanvas = document.createElement('canvas');
    fftCanvas.width = 1024;
    fftCanvas.height = 640;
    fftCtx = fftCanvas.getContext('2d');

    fftTexture = new THREE.CanvasTexture(fftCanvas);
    fftTexture.minFilter = THREE.LinearFilter;
    fftTexture.magFilter = THREE.LinearFilter;

    const screenGeo = new THREE.PlaneGeometry(1.15, 0.58);
    const screenMat = new THREE.MeshBasicMaterial({ map: fftTexture });
    const screenMesh = new THREE.Mesh(screenGeo, screenMat);
    screenMesh.position.set(0, 0, 0.205);
    fftGroup.add(screenMesh);

    scene.add(fftGroup);
  }

  function buildFilterCircuitBoard() {
    const pcbGroup = new THREE.Group();
    pcbGroup.position.set(0.0, 1.21, 0.45);

    const pcbGeo = new THREE.BoxGeometry(1.4, 0.02, 0.75);
    const pcbMat = new THREE.MeshStandardMaterial({
      color: 0x064e3b,
      roughness: 0.4,
      metalness: 0.2
    });
    const pcb = new THREE.Mesh(pcbGeo, pcbMat);
    pcb.receiveShadow = true;
    pcb.castShadow = true;
    pcbGroup.add(pcb);

    const trenchGeo = new THREE.BoxGeometry(0.025, 0.025, 0.72);
    const trenchMat = new THREE.MeshBasicMaterial({ color: 0x00e5ff });
    const trench = new THREE.Mesh(trenchGeo, trenchMat);
    trench.position.set(-0.15, 0.005, 0);
    pcbGroup.add(trench);

    const barrierGeo = new THREE.PlaneGeometry(0.72, 0.22);
    const barrierMat = new THREE.MeshBasicMaterial({
      color: 0x00e5ff,
      transparent: true,
      opacity: 0.35,
      side: THREE.DoubleSide
    });
    const barrier = new THREE.Mesh(barrierGeo, barrierMat);
    barrier.rotation.y = Math.PI / 2;
    barrier.position.set(-0.15, 0.11, 0);
    pcbGroup.add(barrier);

    const optoGeo = new THREE.BoxGeometry(0.12, 0.04, 0.08);
    const icMat = new THREE.MeshStandardMaterial({ color: 0x0f172a, roughness: 0.3 });
    const optoIC = new THREE.Mesh(optoGeo, icMat);
    optoIC.position.set(-0.15, 0.03, 0);
    pcbGroup.add(optoIC);

    const opampGeo = new THREE.BoxGeometry(0.14, 0.04, 0.09);
    const opamp = new THREE.Mesh(opampGeo, icMat);
    opamp.position.set(0.25, 0.03, -0.15);
    pcbGroup.add(opamp);

    const dspGeo = new THREE.BoxGeometry(0.18, 0.035, 0.18);
    const dsp = new THREE.Mesh(dspGeo, icMat);
    dsp.position.set(0.25, 0.03, 0.15);
    pcbGroup.add(dsp);

    const capGeo = new THREE.BoxGeometry(0.05, 0.06, 0.03);
    const capMat = new THREE.MeshStandardMaterial({ color: 0xef4444, roughness: 0.2 });
    for (let c = 0; c < 3; c++) {
      const cap = new THREE.Mesh(capGeo, capMat);
      cap.position.set(0.05 + c * 0.06, 0.04, -0.22);
      pcbGroup.add(cap);
    }

    const termGeo = new THREE.BoxGeometry(0.08, 0.06, 0.18);
    const termMat = new THREE.MeshStandardMaterial({ color: 0x15803d, roughness: 0.5 });

    const inTerm = new THREE.Mesh(termGeo, termMat);
    inTerm.position.set(-0.62, 0.04, 0);
    pcbGroup.add(inTerm);

    const outTerm = new THREE.Mesh(termGeo, termMat);
    outTerm.position.set(0.62, 0.04, 0);
    pcbGroup.add(outTerm);

    const ledGeo = new THREE.SphereGeometry(0.018, 12, 12);
    
    const ledPwr = new THREE.Mesh(ledGeo, new THREE.MeshBasicMaterial({ color: 0x00ff88 }));
    ledPwr.position.set(-0.55, 0.03, 0.28);
    pcbGroup.add(ledPwr);

    pcbLedIso = new THREE.Mesh(ledGeo, new THREE.MeshBasicMaterial({ color: 0x00e5ff }));
    pcbLedIso.position.set(-0.15, 0.03, 0.28);
    pcbGroup.add(pcbLedIso);

    pcbLedFilter = new THREE.Mesh(ledGeo, new THREE.MeshBasicMaterial({ color: 0x00ff88 }));
    pcbLedFilter.position.set(0.35, 0.03, 0.28);
    pcbGroup.add(pcbLedFilter);

    pcbLedClip = new THREE.Mesh(ledGeo, new THREE.MeshBasicMaterial({ color: 0x330000 }));
    pcbLedClip.position.set(-0.4, 0.03, 0.28);
    pcbGroup.add(pcbLedClip);

    const pulseGeo = new THREE.SphereGeometry(0.012, 8, 8);
    for (let p = 0; p < 12; p++) {
      const pulseMat = new THREE.MeshBasicMaterial({ color: 0xffb703 });
      const pulse = new THREE.Mesh(pulseGeo, pulseMat);
      pulse.position.set(-0.6 + p * 0.1, 0.025, 0);
      pcbGroup.add(pulse);
      signalPulses.push({
        mesh: pulse,
        progress: p / 12,
        speed: 0.4
      });
    }

    scene.add(pcbGroup);
  }

  function buildPatientPhantom() {
    const patientGroup = new THREE.Group();
    patientGroup.position.set(-1.8, 1.25, 0.5);

    const torsoGeo = new THREE.CylinderGeometry(0.24, 0.28, 0.65, 20);
    const torsoMat = new THREE.MeshStandardMaterial({
      color: 0x334155,
      roughness: 0.6,
      metalness: 0.1
    });
    const torso = new THREE.Mesh(torsoGeo, torsoMat);
    torso.position.y = 0.35;
    torso.castShadow = true;
    patientGroup.add(torso);

    const headGeo = new THREE.SphereGeometry(0.14, 16, 16);
    const head = new THREE.Mesh(headGeo, torsoMat);
    head.position.y = 0.78;
    patientGroup.add(head);

    const elecGeo = new THREE.CylinderGeometry(0.025, 0.025, 0.015, 12);
    
    const ra = new THREE.Mesh(elecGeo, new THREE.MeshStandardMaterial({ color: 0xff3366 }));
    ra.rotation.x = Math.PI / 2;
    ra.position.set(-0.16, 0.45, 0.22);
    patientGroup.add(ra);

    const la = new THREE.Mesh(elecGeo, new THREE.MeshStandardMaterial({ color: 0xffb703 }));
    la.rotation.x = Math.PI / 2;
    la.position.set(0.16, 0.45, 0.22);
    patientGroup.add(la);

    const ll = new THREE.Mesh(elecGeo, new THREE.MeshStandardMaterial({ color: 0x00ff88 }));
    ll.rotation.x = Math.PI / 2;
    ll.position.set(0.08, 0.20, 0.25);
    patientGroup.add(ll);

    const createCable = (from, to, color) => {
      const curve = new THREE.QuadraticBezierCurve3(
        from,
        new THREE.Vector3((from.x + to.x) / 2, Math.min(from.y, to.y) - 0.15, (from.z + to.z) / 2),
        to
      );
      const tubeGeo = new THREE.TubeGeometry(curve, 16, 0.008, 8, false);
      const tubeMat = new THREE.MeshStandardMaterial({ color, roughness: 0.5 });
      return new THREE.Mesh(tubeGeo, tubeMat);
    };

    const cable1 = createCable(
      new THREE.Vector3(-1.8 - 0.16, 1.70, 0.72),
      new THREE.Vector3(-0.62, 1.25, 0.45),
      0xff3366
    );
    scene.add(cable1);

    const cable2 = createCable(
      new THREE.Vector3(-1.8 + 0.16, 1.70, 0.72),
      new THREE.Vector3(-0.62, 1.25, 0.45),
      0xffb703
    );
    scene.add(cable2);

    scene.add(patientGroup);
  }

  function buildPowerLineConduit() {
    const powerGroup = new THREE.Group();
    powerGroup.position.set(0, 4.2, 0.5);

    const pipeGeo = new THREE.CylinderGeometry(0.04, 0.04, 6.0, 16);
    const pipeMat = new THREE.MeshStandardMaterial({
      color: 0x64748b,
      metalness: 0.8,
      roughness: 0.2
    });
    const pipe = new THREE.Mesh(pipeGeo, pipeMat);
    pipe.rotation.z = Math.PI / 2;
    powerGroup.add(pipe);

    const markerCanvas = document.createElement('canvas');
    markerCanvas.width = 256; markerCanvas.height = 64;
    const mctx = markerCanvas.getContext('2d');
    mctx.fillStyle = '#ffb703';
    mctx.fillRect(0, 0, 256, 64);
    mctx.fillStyle = '#000000';
    mctx.font = 'bold 24px monospace';
    mctx.fillText('230V AC 50Hz', 35, 40);

    const markerTex = new THREE.CanvasTexture(markerCanvas);
    const markerMesh = new THREE.Mesh(
      new THREE.PlaneGeometry(0.8, 0.2),
      new THREE.MeshBasicMaterial({ map: markerTex })
    );
    markerMesh.position.set(0, -0.1, 0.05);
    powerGroup.add(markerMesh);

    for (let r = 0; r < 4; r++) {
      const ringGeo = new THREE.TorusGeometry(0.4 + r * 0.45, 0.015, 8, 32);
      const ringMat = new THREE.MeshBasicMaterial({
        color: 0xffb703,
        transparent: true,
        opacity: 0.35
      });
      const ring = new THREE.Mesh(ringGeo, ringMat);
      ring.rotation.x = Math.PI / 2;
      powerGroup.add(ring);
      emRings.push({ mesh: ring, baseScale: 1.0 + r * 0.5, phase: r * (Math.PI / 2) });
    }

    scene.add(powerGroup);
  }

  // --- Dynamic 2D Canvas Renders for 3D & 2D Instruments ---

  /**
   * Universal 3-Signal Oscilloscope Renderer
   */
  function drawOscilloscope(ctx, w, h) {
    if (!ctx) return;

    // Background
    ctx.fillStyle = '#050c18';
    ctx.fillRect(0, 0, w, h);

    // Reticle Grid (10 divisions X, 8 divisions Y)
    const cols = 10;
    const rows = 8;
    const dx = w / cols;
    const dy = h / rows;

    ctx.strokeStyle = 'rgba(0, 229, 255, 0.12)';
    ctx.lineWidth = 1;

    for (let x = 0; x <= cols; x++) {
      ctx.beginPath();
      ctx.moveTo(x * dx, 0);
      ctx.lineTo(x * dx, h);
      ctx.stroke();
    }

    for (let y = 0; y <= rows; y++) {
      ctx.beginPath();
      ctx.moveTo(0, y * dy);
      ctx.lineTo(w, y * dy);
      ctx.stroke();
    }

    // Center axes tick marks
    ctx.strokeStyle = 'rgba(0, 229, 255, 0.3)';
    ctx.beginPath();
    ctx.moveTo(w / 2, 0); ctx.lineTo(w / 2, h);
    ctx.moveTo(0, h / 2); ctx.lineTo(w, h / 2);
    ctx.stroke();

    // Top Header Readout
    ctx.fillStyle = state.isPaused ? '#ffb703' : '#00e5ff';
    ctx.font = 'bold 20px "JetBrains Mono", monospace';
    ctx.fillText(state.isPaused ? 'DSO PAUSED' : 'TRIG\'D AUTO', 25, 34);

    ctx.fillStyle = '#8b9bb4';
    ctx.font = '16px "JetBrains Mono", monospace';
    ctx.fillText(`TB: ${state.timebase} ms/div`, 220, 34);
    ctx.fillText(`fs: 1000 Hz`, 420, 34);

    // Channel Legends
    ctx.font = '14px "JetBrains Mono", monospace';
    if (state.showCH1) {
      ctx.fillStyle = '#00e5ff';
      ctx.fillText(`CH1 (REF): ${state.scaleCH} mV/div`, 25, h - 18);
    }
    if (state.showCH2) {
      ctx.fillStyle = '#ffb703';
      ctx.fillText(`CH2 (IN): ${state.scaleCH} mV/div`, 260, h - 18);
    }
    if (state.showCH3) {
      ctx.fillStyle = '#00ff88';
      ctx.fillText(`CH3 (FLT): ${state.scaleCH} mV/div`, 500, h - 18);
    }

    // Plot Traces from circular buffer
    const samplesToShow = Math.min(BUFFER_SIZE, Math.round((state.timebase * cols)));
    const startIdx = (bufferIndex - samplesToShow + BUFFER_SIZE) % BUFFER_SIZE;

    const yCenter = h * 0.5;
    const vScale = (dy / state.scaleCH);

    // Trace 1: Clean Reference (Cyan)
    if (state.showCH1) {
      ctx.strokeStyle = '#00e5ff';
      ctx.lineWidth = 2.0;
      ctx.shadowColor = 'rgba(0, 229, 255, 0.7)';
      ctx.shadowBlur = 4;
      ctx.beginPath();
      for (let i = 0; i < samplesToShow; i++) {
        const idx = (startIdx + i) % BUFFER_SIZE;
        const x = (i / samplesToShow) * w;
        const y = yCenter - cleanRefBuffer[idx] * vScale;
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }

    // Trace 2: Raw Noisy Input (Yellow)
    if (state.showCH2) {
      ctx.strokeStyle = '#ffb703';
      ctx.lineWidth = 2.4;
      ctx.shadowColor = 'rgba(255, 183, 3, 0.8)';
      ctx.shadowBlur = 6;
      ctx.beginPath();
      for (let i = 0; i < samplesToShow; i++) {
        const idx = (startIdx + i) % BUFFER_SIZE;
        const x = (i / samplesToShow) * w;
        const y = yCenter - rawBuffer[idx] * vScale;
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }

    // Trace 3: Filtered Output (Neon Green)
    if (state.showCH3) {
      ctx.strokeStyle = '#00ff88';
      ctx.lineWidth = 2.8;
      ctx.shadowColor = 'rgba(0, 255, 136, 0.9)';
      ctx.shadowBlur = 8;
      ctx.beginPath();
      for (let i = 0; i < samplesToShow; i++) {
        const idx = (startIdx + i) % BUFFER_SIZE;
        const x = (i / samplesToShow) * w;
        const y = yCenter - filteredBuffer[idx] * vScale;
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }

    ctx.shadowBlur = 0;
  }

  function renderOscilloscopeScreen() {
    drawOscilloscope(dsoCtx, dsoCanvas.width, dsoCanvas.height);
    dsoTexture.needsUpdate = true;
  }

  /**
   * Universal FFT Spectrum Analyzer Renderer
   */
  function drawSpectrum(ctx, w, h) {
    if (!ctx) return;

    ctx.fillStyle = '#060a14';
    ctx.fillRect(0, 0, w, h);

    const dbMin = -80;
    const dbMax = 10;
    const fMax = 150; // 0 to 150 Hz

    ctx.strokeStyle = 'rgba(0, 229, 255, 0.12)';
    ctx.lineWidth = 1;

    // Horizontal dB lines (-60, -45, -30, -15, 0 dB)
    for (let db = dbMin; db <= dbMax; db += 15) {
      const y = h - ((db - dbMin) / (dbMax - dbMin)) * (h - 85) - 45;
      ctx.beginPath();
      ctx.moveTo(55, y);
      ctx.lineTo(w - 20, y);
      ctx.stroke();

      ctx.fillStyle = '#8b9bb4';
      ctx.font = '13px "JetBrains Mono", monospace';
      ctx.fillText(`${db}dB`, 10, y + 4);
    }

    // Vertical Frequency lines (25, 50, 75, 100, 125, 150 Hz)
    for (let f = 25; f <= fMax; f += 25) {
      const x = 55 + (f / fMax) * (w - 75);
      ctx.beginPath();
      ctx.moveTo(x, 40);
      ctx.lineTo(x, h - 45);
      ctx.stroke();

      ctx.fillStyle = f === 50 ? '#00e5ff' : '#8b9bb4';
      ctx.font = f === 50 ? 'bold 15px "JetBrains Mono"' : '13px "JetBrains Mono"';
      ctx.fillText(`${f}Hz`, x - 16, h - 24);
    }

    // Title
    ctx.fillStyle = '#00e5ff';
    ctx.font = 'bold 20px "JetBrains Mono", monospace';
    ctx.fillText('REAL-TIME FFT POWER SPECTRUM (Hanning Window, 1024-pt)', 55, 30);

    // Plot Spectrum curves from actual FFT calculation
    const binsToPlot = Math.round((fMax / (SAMPLING_RATE / 2)) * (fftEngine.size / 2));

    // Raw Spectrum (Amber / Yellow)
    ctx.strokeStyle = 'rgba(255, 183, 3, 0.8)';
    ctx.lineWidth = 2.0;
    ctx.beginPath();
    for (let i = 0; i < binsToPlot; i++) {
      const freq = (i / (fftEngine.size / 2)) * (SAMPLING_RATE / 2);
      const x = 55 + (freq / fMax) * (w - 75);
      const db = Math.max(dbMin, Math.min(dbMax, latestRawSpectrum[i]));
      const y = h - ((db - dbMin) / (dbMax - dbMin)) * (h - 85) - 45;
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // Filtered Spectrum (Neon Green)
    ctx.strokeStyle = '#00ff88';
    ctx.lineWidth = 2.6;
    ctx.shadowColor = 'rgba(0, 255, 136, 0.8)';
    ctx.shadowBlur = 6;
    ctx.beginPath();
    for (let i = 0; i < binsToPlot; i++) {
      const freq = (i / (fftEngine.size / 2)) * (SAMPLING_RATE / 2);
      const x = 55 + (freq / fMax) * (w - 75);
      const db = Math.max(dbMin, Math.min(dbMax, latestFltSpectrum[i]));
      const y = h - ((db - dbMin) / (dbMax - dbMin)) * (h - 85) - 45;
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    }
    ctx.stroke();
    ctx.shadowBlur = 0;

    // Highlight Notch Center Frequency Marker
    const markerX = 55 + (state.notchFreq / fMax) * (w - 75);
    ctx.strokeStyle = '#00e5ff';
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(markerX, 40);
    ctx.lineTo(markerX, h - 45);
    ctx.stroke();
    ctx.setLineDash([]);

    // Notch Tag Readout
    ctx.fillStyle = '#00e5ff';
    ctx.font = 'bold 14px "JetBrains Mono"';
    ctx.fillText(`f₀: ${state.notchFreq}Hz | Attenuation: ${state.metrics.attenuation} dB`, Math.min(w - 320, markerX + 8), 75);
  }

  function renderSpectrumScreen() {
    drawSpectrum(fftCtx, fftCanvas.width, fftCanvas.height);
    fftTexture.needsUpdate = true;
  }

  /**
   * Renders Expanded 2D Modal Canvases
   */
  function renderExpandedInstruments() {
    const dsoExp = document.getElementById('canvas-expanded-dso');
    if (dsoExp && dsoExp.offsetParent !== null) {
      drawOscilloscope(dsoExp.getContext('2d'), dsoExp.width, dsoExp.height);
    }
    const fftExp = document.getElementById('canvas-expanded-fft');
    if (fftExp && fftExp.offsetParent !== null) {
      drawSpectrum(fftExp.getContext('2d'), fftExp.width, fftExp.height);
    }
  }

  // --- Animation & Render Loop ---
  let lastTime = performance.now();
  let frameCount = 0;
  let fpsTimer = 0;

  function animate(now) {
    const dt = Math.min(0.1, (now - lastTime) / 1000);
    lastTime = now;

    // Compute FPS
    frameCount++;
    fpsTimer += dt;
    if (fpsTimer >= 0.5) {
      state.metrics.fps = Math.round(frameCount / fpsTimer);
      document.getElementById('fps-counter').innerText = `${state.metrics.fps} FPS`;
      frameCount = 0;
      fpsTimer = 0;
    }

    // Step physiological DSP simulation (approx 16 samples per 60fps frame = 1000 Hz)
    const samplesThisFrame = Math.round(SAMPLING_RATE * dt);
    stepSignalProcessing(samplesThisFrame);

    // Update 3D Canvas Textures
    renderOscilloscopeScreen();
    renderSpectrumScreen();

    // Update Expanded 2D Instrument View if open
    renderExpandedInstruments();

    // Animate Electromagnetic Field Waves
    if (state.emFieldVisible) {
      emRings.forEach(item => {
        item.mesh.visible = true;
        item.phase += 2 * Math.PI * state.mainsFreq * dt * 0.05;
        const scale = item.baseScale + 0.15 * Math.sin(item.phase);
        item.mesh.scale.set(scale, scale, scale);
        item.mesh.material.opacity = (state.noiseAmp / 3.0) * 0.45;
      });
    } else {
      emRings.forEach(item => { item.mesh.visible = false; });
    }

    // Animate PCB Signal Pulses
    signalPulses.forEach(item => {
      item.progress = (item.progress + item.speed * dt) % 1.0;
      const x = -0.6 + item.progress * 1.2;
      item.mesh.position.x = x;
      
      if (x < -0.15) {
        item.mesh.material.color.setHex(0xffb703); // Corrupted
      } else if (x < 0.2) {
        item.mesh.material.color.setHex(0x00e5ff); // Galvanically Isolated
      } else {
        item.mesh.material.color.setHex(state.filterActive ? 0x00ff88 : 0xffb703); // Filtered Clean
      }
    });

    // Update PCB LEDs
    if (pcbLedIso) {
      pcbLedIso.material.color.setHex(state.isolationActive ? 0x00e5ff : 0x072738);
    }
    if (pcbLedFilter) {
      pcbLedFilter.material.color.setHex(state.filterActive ? 0x00ff88 : 0x06331a);
    }
    if (pcbLedClip) {
      pcbLedClip.material.color.setHex(state.noiseAmp > 1.8 ? 0xff3366 : 0x330000);
    }

    // Smooth Camera Transition when a preset or lab step is clicked
    if (isTransitioningCamera && controls) {
      camera.position.lerp(targetCamPos, 0.06);
      controls.target.lerp(targetCamLook, 0.06);

      if (camera.position.distanceTo(targetCamPos) < 0.02 && controls.target.distanceTo(targetCamLook) < 0.02) {
        camera.position.copy(targetCamPos);
        controls.target.copy(targetCamLook);
        isTransitioningCamera = false;
      }
    }
    if (controls) {
      controls.update();
    }

    // Update audio hum gain if enabled
    if (state.audioActive) {
      updateAudioGain();
    }

    // Render WebGL / WebXR frame
    renderer.render(scene, camera);
  }

  function onWindowResize() {
    const container = document.getElementById('canvas-container');
    const width = container.clientWidth || window.innerWidth;
    const height = container.clientHeight || window.innerHeight;
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
    renderer.setSize(width, height);
  }

  // --- UI DOM Synchronization & Event Handlers ---

  function updateHUDDOM() {
    const snrInEl = document.getElementById('metric-snr-in');
    const snrOutEl = document.getElementById('metric-snr-out');
    const noiseRMSEl = document.getElementById('metric-noise-rms');
    const errorRMSEl = document.getElementById('metric-error-rms');
    const redEl = document.getElementById('metric-reduction');
    const attenEl = document.getElementById('metric-attenuation');
    const resEl = document.getElementById('metric-residual');
    const isoBadge = document.getElementById('badge-iso-status');
    const fltBadge = document.getElementById('badge-filter-status');
    const bwBadge = document.getElementById('badge-bw-status');

    if (snrInEl) snrInEl.innerHTML = `${state.metrics.snrIn} <span class="unit">dB</span>`;
    if (snrOutEl) {
      snrOutEl.innerHTML = `${state.metrics.snrOut} <span class="unit">dB</span>`;
      snrOutEl.className = parseFloat(state.metrics.snrOut) > 15 ? 'metric-value snr-good' : 'metric-value snr-bad';
    }
    if (noiseRMSEl) noiseRMSEl.innerHTML = `${state.metrics.noiseRMS} <span class="unit">µV</span>`;
    if (errorRMSEl) errorRMSEl.innerHTML = `${state.metrics.errorRMS} <span class="unit">µV</span>`;
    if (redEl) redEl.innerHTML = `${state.metrics.reductionPct} <span class="unit">%</span>`;
    if (attenEl) attenEl.innerHTML = `${state.metrics.attenuation} <span class="unit">dB</span>`;
    if (resEl) resEl.innerHTML = `${state.metrics.residual50Hz} <span class="unit">µV</span>`;

    if (isoBadge) {
      isoBadge.innerText = state.isolationActive ? 'ACTIVE (ISOLATED)' : 'BYPASS (GROUND LOOP)';
      isoBadge.className = state.isolationActive ? 'status-tag active' : 'status-tag inactive';
    }

    if (fltBadge) {
      fltBadge.innerText = state.filterActive ? 'FILTERING' : 'BYPASSED';
      fltBadge.className = state.filterActive ? 'status-tag active' : 'status-tag inactive';
    }

    if (bwBadge) {
      const bw = (state.notchFreq / state.qFactor).toFixed(1);
      bwBadge.innerText = `${bw} Hz (Q=${state.qFactor.toFixed(1)})`;
    }

    // Stream status indicator
    const simBadge = document.getElementById('sim-status-badge');
    const simText = document.getElementById('sim-state-text');
    if (simBadge && simText) {
      if (state.isPaused) {
        simBadge.innerHTML = '<span class="dot amber"></span><span>PAUSED</span>';
      } else {
        simBadge.innerHTML = '<span class="dot green"></span><span>LIVE</span>';
      }
    }
  }

  function setupUIEventHandlers() {
    // 1. Camera View Navigation Pills
    const camButtons = document.querySelectorAll('#camera-nav .cam-btn');
    camButtons.forEach(btn => {
      btn.addEventListener('click', () => {
        camButtons.forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        const view = btn.dataset.cam;
        if (cameraPresets[view]) {
          transitionCameraTo(cameraPresets[view].pos, cameraPresets[view].target);
        }
      });
    });

    // 2. Bio-Signal Selection Pills
    const sigPills = document.querySelectorAll('#signal-type-pills .btn-pill');
    sigPills.forEach(pill => {
      pill.addEventListener('click', () => {
        sigPills.forEach(p => p.classList.remove('active'));
        pill.classList.add('active');
        state.signalType = pill.dataset.signal;
        
        const hrGroup = document.getElementById('group-heart-rate');
        if (hrGroup) hrGroup.style.display = state.signalType === 'ecg' ? 'block' : 'none';

        const tag = document.getElementById('sig-type-tag');
        if (tag) tag.innerText = state.signalType.toUpperCase();
      });
    });

    // 3. Heart Rate Slider
    const sliderHR = document.getElementById('slider-hr');
    if (sliderHR) {
      sliderHR.addEventListener('input', (e) => {
        state.heartRate = parseInt(e.target.value);
        document.getElementById('val-hr').innerText = state.heartRate;
      });
    }

    // 4. Noise Amplitude Slider (PLI)
    const sliderAmp = document.getElementById('slider-pli-amp');
    if (sliderAmp) {
      sliderAmp.addEventListener('input', (e) => {
        state.noiseAmp = parseFloat(e.target.value);
        document.getElementById('val-pli-amp').innerText = `${state.noiseAmp.toFixed(2)} mV`;
        updateAudioGain();
      });
    }

    // 5. Mains Grid Frequency Slider
    const sliderMains = document.getElementById('slider-mains-freq');
    if (sliderMains) {
      sliderMains.addEventListener('input', (e) => {
        state.mainsFreq = parseFloat(e.target.value);
        document.getElementById('val-mains-freq').innerText = `${state.mainsFreq.toFixed(1)} Hz`;
        updateBiquadCoefficients();
        updateAudioGain();
      });
    }

    // 6. 3rd Harmonic Toggle
    const checkHarmonic = document.getElementById('check-harmonic');
    if (checkHarmonic) {
      checkHarmonic.addEventListener('change', (e) => {
        state.harmonic3rd = e.target.checked;
        updateAudioGain();
      });
    }

    // 7. EM Radiation Waves Toggle
    const checkEM = document.getElementById('check-em-field');
    if (checkEM) {
      checkEM.addEventListener('change', (e) => {
        state.emFieldVisible = e.target.checked;
      });
    }

    // 8. Galvanic Isolation Barrier Toggle
    const checkIso = document.getElementById('check-isolation');
    if (checkIso) {
      checkIso.addEventListener('change', (e) => {
        state.isolationActive = e.target.checked;
        const cmGroup = document.getElementById('group-ground-loop');
        if (cmGroup) cmGroup.style.display = state.isolationActive ? 'none' : 'block';
        updateHUDDOM();
        updateAudioGain();
      });
    }

    const sliderCM = document.getElementById('slider-cm-noise');
    if (sliderCM) {
      sliderCM.addEventListener('input', (e) => {
        state.groundLoopVoltage = parseFloat(e.target.value);
        document.getElementById('val-cm-noise').innerText = `${state.groundLoopVoltage.toFixed(2)} V`;
      });
    }

    // 9. Real-Time Notch Filter State Toggle
    const checkFilter = document.getElementById('check-filter');
    if (checkFilter) {
      checkFilter.addEventListener('change', (e) => {
        state.filterActive = e.target.checked;
        renderFilterFrequencyResponse();
        updateHUDDOM();
        updateAudioGain();
      });
    }

    // 10. Filter Architecture Pills (IIR, Twin-T, LMS)
    const filterPills = document.querySelectorAll('#filter-type-pills .btn-pill');
    filterPills.forEach(pill => {
      pill.addEventListener('click', () => {
        filterPills.forEach(p => p.classList.remove('active'));
        pill.classList.add('active');
        state.filterType = pill.dataset.filter;
        iirState.w1 = 0; iirState.w2 = 0;
        updateBiquadCoefficients();
      });
    });

    // 11. Notch Center Frequency Slider (f0)
    const sliderNotchFreq = document.getElementById('slider-notch-freq');
    if (sliderNotchFreq) {
      sliderNotchFreq.addEventListener('input', (e) => {
        state.notchFreq = parseFloat(e.target.value);
        document.getElementById('val-notch-freq').innerText = `${state.notchFreq.toFixed(1)} Hz`;
        updateBiquadCoefficients();
      });
    }

    // 12. Quality Factor Slider (Q)
    const sliderQ = document.getElementById('slider-q-factor');
    if (sliderQ) {
      sliderQ.addEventListener('input', (e) => {
        state.qFactor = parseFloat(e.target.value);
        document.getElementById('val-q-factor').innerText = `${state.qFactor.toFixed(1)}`;
        updateBiquadCoefficients();
        updateHUDDOM();
      });
    }

    // 13. Channel Visibility Checkboxes
    const ch1Box = document.getElementById('check-ch1');
    if (ch1Box) {
      ch1Box.addEventListener('change', (e) => { state.showCH1 = e.target.checked; });
    }
    const ch2Box = document.getElementById('check-ch2');
    if (ch2Box) {
      ch2Box.addEventListener('change', (e) => { state.showCH2 = e.target.checked; });
    }
    const ch3Box = document.getElementById('check-ch3');
    if (ch3Box) {
      ch3Box.addEventListener('change', (e) => { state.showCH3 = e.target.checked; });
    }

    // 14. Oscilloscope Timebase & Scale
    const sliderTB = document.getElementById('slider-timebase');
    if (sliderTB) {
      sliderTB.addEventListener('input', (e) => {
        state.timebase = parseInt(e.target.value);
        document.getElementById('val-timebase').innerText = `${state.timebase} ms/div`;
        const instTB = document.getElementById('inst-dso-tb');
        if (instTB) instTB.innerText = `TB: ${state.timebase} ms/div`;
      });
    }

    const sliderCHScale = document.getElementById('slider-scale-ch');
    if (sliderCHScale) {
      sliderCHScale.addEventListener('input', (e) => {
        state.scaleCH = parseFloat(e.target.value);
        document.getElementById('val-scale-ch').innerText = `${state.scaleCH.toFixed(1)} mV/div`;
      });
    }

    // 15. Simulation State Controls (Play, Pause, Reset)
    const btnSimPlay = document.getElementById('btn-sim-play');
    if (btnSimPlay) {
      btnSimPlay.addEventListener('click', () => {
        state.isPaused = false;
        updateHUDDOM();
        updateAudioGain();
      });
    }

    const btnSimPause = document.getElementById('btn-sim-pause');
    if (btnSimPause) {
      btnSimPause.addEventListener('click', () => {
        state.isPaused = true;
        updateHUDDOM();
        updateAudioGain();
      });
    }

    const btnSimReset = document.getElementById('btn-sim-reset');
    if (btnSimReset) {
      btnSimReset.addEventListener('click', () => {
        globalSampleCount = 0;
        bufferIndex = 0;
        rawBuffer.fill(0);
        filteredBuffer.fill(0);
        cleanRefBuffer.fill(0);
        iirState.w1 = 0; iirState.w2 = 0;
      });
    }

    // 16. Audio Hum Button
    const btnAudio = document.getElementById('btn-toggle-audio');
    if (btnAudio) {
      btnAudio.addEventListener('click', () => {
        initAudio();
        state.audioActive = !state.audioActive;
        btnAudio.classList.toggle('active', state.audioActive);
        document.getElementById('audio-btn-text').innerText = state.audioActive ? 'Mute 50Hz Hum' : 'Listen to 50Hz Hum';
        updateAudioGain();
      });
    }

    // 17. Controls Dock Collapse Button
    const btnToggleDock = document.getElementById('btn-toggle-dock');
    const dock = document.getElementById('controls-dock');
    if (btnToggleDock && dock) {
      btnToggleDock.addEventListener('click', () => {
        dock.classList.toggle('collapsed');
        const icon = document.getElementById('dock-toggle-icon');
        icon.className = dock.classList.contains('collapsed') ? 'fa-solid fa-chevron-left' : 'fa-solid fa-chevron-right';
      });
    }

    // 18. Modals Management
    // Theory Modal
    const modalTheory = document.getElementById('modal-theory');
    const btnOpenTheory = document.getElementById('btn-open-theory');
    const btnCloseTheory = document.getElementById('btn-close-theory');

    if (btnOpenTheory && modalTheory) {
      btnOpenTheory.addEventListener('click', () => modalTheory.classList.add('open'));
    }
    if (btnCloseTheory && modalTheory) {
      btnCloseTheory.addEventListener('click', () => modalTheory.classList.remove('open'));
    }

    // Expanded Instruments Modal
    const modalInstruments = document.getElementById('modal-instruments');
    const btnOpenInstruments = document.getElementById('btn-open-instruments');
    const btnCloseInstruments = document.getElementById('btn-close-instruments');

    if (btnOpenInstruments && modalInstruments) {
      btnOpenInstruments.addEventListener('click', () => modalInstruments.classList.add('open'));
    }
    if (btnCloseInstruments && modalInstruments) {
      btnCloseInstruments.addEventListener('click', () => modalInstruments.classList.remove('open'));
    }

    // Channel toggles inside expanded instrument modal
    const btnCh1Toggle = document.getElementById('btn-dso-ch1-toggle');
    if (btnCh1Toggle) {
      btnCh1Toggle.addEventListener('click', () => {
        state.showCH1 = !state.showCH1;
        btnCh1Toggle.classList.toggle('active', state.showCH1);
        if (ch1Box) ch1Box.checked = state.showCH1;
      });
    }
    const btnCh2Toggle = document.getElementById('btn-dso-ch2-toggle');
    if (btnCh2Toggle) {
      btnCh2Toggle.addEventListener('click', () => {
        state.showCH2 = !state.showCH2;
        btnCh2Toggle.classList.toggle('active', state.showCH2);
        if (ch2Box) ch2Box.checked = state.showCH2;
      });
    }
    const btnCh3Toggle = document.getElementById('btn-dso-ch3-toggle');
    if (btnCh3Toggle) {
      btnCh3Toggle.addEventListener('click', () => {
        state.showCH3 = !state.showCH3;
        btnCh3Toggle.classList.toggle('active', state.showCH3);
        if (ch3Box) ch3Box.checked = state.showCH3;
      });
    }

    // Guided Lab Steps Modal (8 Steps)
    const modalGuide = document.getElementById('modal-guide');
    const btnOpenGuide = document.getElementById('btn-open-guide');
    const btnCloseGuide = document.getElementById('btn-close-guide');

    if (btnOpenGuide && modalGuide) {
      btnOpenGuide.addEventListener('click', () => modalGuide.classList.add('open'));
    }
    if (btnCloseGuide && modalGuide) {
      btnCloseGuide.addEventListener('click', () => modalGuide.classList.remove('open'));
    }

    // Guide Step Navigation
    const btnGuideNext = document.getElementById('btn-guide-next');
    const btnGuidePrev = document.getElementById('btn-guide-prev');

    if (btnGuideNext) {
      btnGuideNext.addEventListener('click', () => {
        if (state.guideStep < 8) {
          goToGuideStep(state.guideStep + 1);
        }
      });
    }
    if (btnGuidePrev) {
      btnGuidePrev.addEventListener('click', () => {
        if (state.guideStep > 1) {
          goToGuideStep(state.guideStep - 1);
        }
      });
    }

    // Quick Apply Buttons on each guide step
    const stepButtons = document.querySelectorAll('.btn-apply-step');
    stepButtons.forEach(btn => {
      btn.addEventListener('click', () => {
        const step = parseInt(btn.dataset.step);
        applyLabStep(step);
      });
    });

    // 19. Reset Experiment to Defaults
    const btnReset = document.getElementById('btn-reset-exp');
    if (btnReset) {
      btnReset.addEventListener('click', resetDefaults);
    }

    // 20. Export CSV Data
    const btnCSV = document.getElementById('btn-export-csv');
    if (btnCSV) {
      btnCSV.addEventListener('click', exportCSV);
    }

    // 21. Take PNG Snapshot
    const btnPNG = document.getElementById('btn-export-png');
    if (btnPNG) {
      btnPNG.addEventListener('click', takeSnapshot);
    }

    // 22. Print Lab Report
    const btnPrint1 = document.getElementById('btn-print-report');
    const btnPrint2 = document.getElementById('btn-dock-print-report');
    const btnPrint3 = document.getElementById('btn-guide-print-report');

    [btnPrint1, btnPrint2, btnPrint3].forEach(btn => {
      if (btn) btn.addEventListener('click', generatePrintReport);
    });
  }

  /**
   * Navigation inside Guided 8-Step Mode
   */
  function goToGuideStep(step) {
    state.guideStep = step;
    
    // Update progress bar
    const bar = document.getElementById('guide-progress-bar');
    if (bar) bar.style.width = `${(step / 8) * 100}%`;

    const counter = document.getElementById('guide-step-counter');
    if (counter) counter.innerText = `Step ${step} of 8`;

    // Show active step page
    const pages = document.querySelectorAll('.guide-step-page');
    pages.forEach(p => {
      p.classList.toggle('active', parseInt(p.dataset.step) === step);
    });

    // Button states
    const btnPrev = document.getElementById('btn-guide-prev');
    const btnNext = document.getElementById('btn-guide-next');
    if (btnPrev) btnPrev.disabled = (step === 1);
    if (btnNext) {
      if (step === 8) {
        btnNext.innerHTML = 'Completed <i class="fa-solid fa-check"></i>';
        btnNext.disabled = true;
      } else {
        btnNext.innerHTML = 'Next Step <i class="fa-solid fa-chevron-right"></i>';
        btnNext.disabled = false;
      }
    }

    // If reaching step 8, populate summary table
    if (step === 8) {
      updateGuideSummaryTable();
    }
  }

  function updateGuideSummaryTable() {
    const sSig = document.getElementById('sum-sig');
    const sMains = document.getElementById('sum-mains');
    const sNotch = document.getElementById('sum-notch');
    const sQ = document.getElementById('sum-q');
    const sSnrIn = document.getElementById('sum-snr-in');
    const sSnrOut = document.getElementById('sum-snr-out');
    const sRed = document.getElementById('sum-reduction');
    const sAtten = document.getElementById('sum-atten');

    if (sSig) sSig.innerText = `${state.signalType.toUpperCase()} (${state.heartRate} BPM)`;
    if (sMains) sMains.innerText = `${state.mainsFreq.toFixed(1)} Hz`;
    if (sNotch) sNotch.innerText = `${state.notchFreq.toFixed(1)} Hz`;
    if (sQ) sQ.innerText = `${state.qFactor.toFixed(1)}`;
    if (sSnrIn) sSnrIn.innerText = `${state.metrics.snrIn} dB`;
    if (sSnrOut) sSnrOut.innerText = `${state.metrics.snrOut} dB`;
    if (sRed) sRed.innerText = `${state.metrics.reductionPct}%`;
    if (sAtten) sAtten.innerText = `${state.metrics.attenuation} dB`;
  }

  /**
   * Applies guided experiment step configurations (Steps 1 through 8)
   */
  function applyLabStep(step) {
    if (step === 1) {
      // Step 1: Configure Bio-Signal
      state.signalType = 'ecg';
      state.heartRate = 72;
      state.noiseAmp = 0.0;
      state.isolationActive = true;
      state.filterActive = false;
      transitionCameraTo(cameraPresets.overview.pos, cameraPresets.overview.target);
    } 
    else if (step === 2) {
      // Step 2: Observe Clean Signal
      state.noiseAmp = 0.0;
      state.isolationActive = true;
      state.filterActive = false;
      state.showCH1 = true;
      state.showCH2 = false;
      state.showCH3 = false;
      transitionCameraTo(cameraPresets.dso.pos, cameraPresets.dso.target);
    } 
    else if (step === 3) {
      // Step 3: Introduce 50Hz PLI
      state.noiseAmp = 1.20;
      state.mainsFreq = 50.0;
      state.harmonic3rd = true;
      state.isolationActive = false; // Ground loop active
      state.filterActive = false;
      state.showCH1 = true;
      state.showCH2 = true;
      state.showCH3 = false;
      transitionCameraTo(cameraPresets.dso.pos, cameraPresets.dso.target);
    } 
    else if (step === 4) {
      // Step 4: Inspect FFT Spectrum
      transitionCameraTo(cameraPresets.spectrum.pos, cameraPresets.spectrum.target);
    } 
    else if (step === 5) {
      // Step 5: Enable Notch Filter
      state.isolationActive = true;
      state.filterActive = true;
      state.notchFreq = 50.0;
      state.qFactor = 16.0;
      transitionCameraTo(cameraPresets.circuit.pos, cameraPresets.circuit.target);
    } 
    else if (step === 6) {
      // Step 6: Compare Input and Filtered
      state.showCH1 = true;
      state.showCH2 = true;
      state.showCH3 = true;
      transitionCameraTo(cameraPresets.dso.pos, cameraPresets.dso.target);
    } 
    else if (step === 7) {
      // Step 7: Evaluate Performance
      transitionCameraTo(cameraPresets.overview.pos, cameraPresets.overview.target);
    } 
    else if (step === 8) {
      // Step 8: Summary
      updateGuideSummaryTable();
    }

    // Synchronize UI Controls
    const checkIso = document.getElementById('check-isolation');
    if (checkIso) checkIso.checked = state.isolationActive;
    const checkFlt = document.getElementById('check-filter');
    if (checkFlt) checkFlt.checked = state.filterActive;
    const sliderAmp = document.getElementById('slider-pli-amp');
    if (sliderAmp) {
      sliderAmp.value = state.noiseAmp;
      document.getElementById('val-pli-amp').innerText = `${state.noiseAmp.toFixed(2)} mV`;
    }

    updateBiquadCoefficients();
    updateHUDDOM();
    updateAudioGain();
  }

  function resetDefaults() {
    state.signalType = 'ecg';
    state.heartRate = 72;
    state.noiseAmp = 1.20;
    state.mainsFreq = 50.0;
    state.harmonic3rd = true;
    state.emFieldVisible = true;
    state.isolationActive = true;
    state.filterActive = true;
    state.filterType = 'iir';
    state.notchFreq = 50.0;
    state.qFactor = 16.0;
    state.timebase = 100;
    state.scaleCH = 0.5;
    state.showCH1 = true;
    state.showCH2 = true;
    state.showCH3 = true;
    state.isPaused = false;
    state.audioActive = false;
    state.guideStep = 1;

    // Reset controls UI
    const checkIso = document.getElementById('check-isolation');
    if (checkIso) checkIso.checked = true;
    const checkFlt = document.getElementById('check-filter');
    if (checkFlt) checkFlt.checked = true;
    const sliderAmp = document.getElementById('slider-pli-amp');
    if (sliderAmp) {
      sliderAmp.value = 1.20;
      document.getElementById('val-pli-amp').innerText = '1.20 mV';
    }
    const sliderNotch = document.getElementById('slider-notch-freq');
    if (sliderNotch) {
      sliderNotch.value = 50.0;
      document.getElementById('val-notch-freq').innerText = '50.0 Hz';
    }
    const sliderQ = document.getElementById('slider-q-factor');
    if (sliderQ) {
      sliderQ.value = 16.0;
      document.getElementById('val-q-factor').innerText = '16.0';
    }

    transitionCameraTo(cameraPresets.overview.pos, cameraPresets.overview.target);
    goToGuideStep(1);

    updateBiquadCoefficients();
    updateHUDDOM();
    updateAudioGain();
  }

  /**
   * Export Numerical CSV Data
   */
  function exportCSV() {
    let csv = `# Experiment: 50Hz Power Line Interference Isolation & Filter Virtual Lab\n`;
    csv += `# Medical Standard: IEC 60601-1\n`;
    csv += `# Export Date: ${new Date().toISOString()}\n`;
    csv += `# Bio-Signal Source: ${state.signalType.toUpperCase()} (${state.heartRate} BPM)\n`;
    csv += `# Sampling Frequency: ${SAMPLING_RATE} Hz\n`;
    csv += `# Mains Power Frequency: ${state.mainsFreq.toFixed(1)} Hz\n`;
    csv += `# Digital Notch Center (f0): ${state.notchFreq.toFixed(1)} Hz | Q: ${state.qFactor.toFixed(1)}\n`;
    csv += `# Input SNR: ${state.metrics.snrIn} dB | Output SNR: ${state.metrics.snrOut} dB\n`;
    csv += `# Noise Reduction: ${state.metrics.reductionPct}%\n`;
    csv += `# 50Hz Attenuation: ${state.metrics.attenuation} dB\n`;
    csv += `# Data Disclaimer: Simulated Educational Biomedical Data\n`;
    csv += `Time_ms,Clean_Reference_mV,Noisy_Input_mV,Filtered_Output_mV,Error_mV\n`;

    const n = BUFFER_SIZE;
    for (let i = 0; i < n; i++) {
      const idx = (bufferIndex + i) % BUFFER_SIZE;
      const clean = cleanRefBuffer[idx];
      const raw = rawBuffer[idx];
      const flt = filteredBuffer[idx];
      const err = flt - clean;
      csv += `${i},${clean.toFixed(4)},${raw.toFixed(4)},${flt.toFixed(4)},${err.toFixed(4)}\n`;
    }

    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `50Hz_Isolation_Filter_Data_${Date.now()}.csv`;
    a.click();
    URL.revokeObjectURL(url);
  }

  function takeSnapshot() {
    renderer.render(scene, camera);
    const dataURL = renderer.domElement.toDataURL('image/png');
    const a = document.createElement('a');
    a.href = dataURL;
    a.download = `VirtualLab_Snapshot_${Date.now()}.png`;
    a.click();
  }

  /**
   * Generates Print-Friendly Laboratory PDF Report
   */
  function generatePrintReport() {
    // Fill parameters
    const pSig = document.getElementById('p-sig-type');
    const pHr = document.getElementById('p-hr');
    const pAmp = document.getElementById('p-noise-amp');
    const pMains = document.getElementById('p-mains-freq');
    const pHarm = document.getElementById('p-harmonic');
    const pIso = document.getElementById('p-iso');
    const pState = document.getElementById('p-flt-state');
    const pNotch = document.getElementById('p-notch-freq');
    const pQ = document.getElementById('p-q');
    const pDate = document.getElementById('print-date');

    if (pSig) pSig.innerText = state.signalType.toUpperCase();
    if (pHr) pHr.innerText = state.heartRate;
    if (pAmp) pAmp.innerText = `${state.noiseAmp.toFixed(2)} mV`;
    if (pMains) pMains.innerText = `${state.mainsFreq.toFixed(1)} Hz`;
    if (pHarm) pHarm.innerText = state.harmonic3rd ? 'Enabled (150 Hz)' : 'Disabled';
    if (pIso) pIso.innerText = state.isolationActive ? 'Active (> 120 dB CMRR)' : 'Bypassed (Ground Loop)';
    if (pState) pState.innerText = state.filterActive ? 'Active (Filtering)' : 'Bypassed';
    if (pNotch) pNotch.innerText = `${state.notchFreq.toFixed(1)} Hz`;
    if (pQ) pQ.innerText = state.qFactor.toFixed(1);
    if (pDate) pDate.innerText = `Date: ${new Date().toLocaleDateString()} ${new Date().toLocaleTimeString()}`;

    // Fill metrics
    const mNoise = document.getElementById('m-noise-rms');
    const mErr = document.getElementById('m-error-rms');
    const mSnrIn = document.getElementById('m-snr-in');
    const mSnrOut = document.getElementById('m-snr-out');
    const mRed = document.getElementById('m-reduction');
    const mAtt = document.getElementById('m-atten');

    if (mNoise) mNoise.innerText = `${state.metrics.noiseRMS} µV`;
    if (mErr) mErr.innerText = `${state.metrics.errorRMS} µV`;
    if (mSnrIn) mSnrIn.innerText = `${state.metrics.snrIn} dB`;
    if (mSnrOut) mSnrOut.innerText = `${state.metrics.snrOut} dB`;
    if (mRed) mRed.innerText = `${state.metrics.reductionPct}%`;
    if (mAtt) mAtt.innerText = `${state.metrics.attenuation} dB`;

    // Capture DSO Snapshot image into report
    const dsoDataURL = dsoCanvas.toDataURL('image/png');
    const imgEl = document.getElementById('print-snapshot-img');
    if (imgEl) imgEl.src = dsoDataURL;

    // Student notes
    const notesInput = document.getElementById('student-observations');
    const notesPrint = document.getElementById('print-student-notes');
    if (notesInput && notesPrint && notesInput.value.trim() !== '') {
      notesPrint.innerText = notesInput.value;
    }

    // Trigger Print Dialog (allowing user to Save as PDF)
    window.print();
  }

  // --- Bootstrap Initialization ---
  window.addEventListener('DOMContentLoaded', () => {
    updateBiquadCoefficients();
    initThreeJS();
    setupUIEventHandlers();
    updateHUDDOM();
    renderFilterFrequencyResponse();

    // Start WebGL / WebXR Animation Loop
    renderer.setAnimationLoop(animate);
  });

})();
